#include "ExternalProblem.h"
#include "FispactProblem.h"

#include "../../include/userobjects/FispactNuclearDataPaths.C"
/// Custom user object includes
#include "FispactSchedule.h"

#include "HDF5Utils.h"
#include "MooseError.h"
#include "MooseTypes.h"

/// Fispact includes

//// PugiXML include
#include "PhotonSharingData.h"
#include "fispactutil.hpp"
#include "libmesh/id_types.h"
#include "pugixml.hpp"

/// Cpp includes
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <mpi.h>

// Avogadro's number
#define AVOGADRO 6.0221408e+23

#define MOLAR_MASS_DATASET_DIMS 1

registerMooseObject("FizzyApp", FispactProblem);

InputParameters FispactProblem::validParams() {
  InputParameters params = ExternalProblem::validParams();

  /// Suppress parameters we are not using
  params.suppressParameter<bool>("allow_invalid_solution");
  params.suppressParameter<bool>("boundary_restricted_elem_integrity_check");
  params.suppressParameter<bool>("boundary_restricted_node_integrity_check");
  params.suppressParameter<bool>("check_uo_aux_state");
  params.suppressParameter<bool>("error_on_jacobian_nonzero_reallocation");
  params.suppressParameter<std::vector<std::vector<TagName>>>(
      "extra_tag_matrices");
  params.suppressParameter<std::vector<TagName>>("extra_tag_solutions");
  params.suppressParameter<std::vector<std::vector<TagName>>>(
      "extra_tag_vectors");
  params.suppressParameter<bool>("force_restart");
  params.suppressParameter<bool>("fv_bcs_integrity_check");
  params.suppressParameter<bool>("material_dependency_check");
  params.suppressParameter<unsigned int>("near_null_space_dimension");
  params.suppressParameter<unsigned int>("null_space_dimension");
  params.suppressParameter<unsigned int>("transpose_null_space_dimension");
  params.suppressParameter<bool>("immediately_print_invalid_solution");
  params.suppressParameter<bool>("identify_variable_groups_in_nl");
  params.suppressParameter<std::vector<LinearSystemName>>("linear_sys_names");

  /// New parameters we need for FISPACT
  params.addRequiredParam<FileName>("neutron_flux_file",
                                    "HDF5 file storing neutron flux");
  params.addRequiredParam<int>(
      "neutron_flux_tally_id",
      "Path within HDF5 file for the vector storing neutron flux");

  params.addRequiredParam<std::string>(
      "fispact_nuclear_data_uo", "Name of the FispactNuclearDataPaths objects "
                                 "to use for setting nuclear data");

  params.addRequiredParam<std::string>("neutron_bin_type",
                                       "neutron binning scheme for FISPACT");

  params.addRequiredParam<UserObjectName>(
      "fispact_schedule_uo", "Name of the FispactSchedule user object defining "
                             "the FISPACT flux schedule");

  params.addParam<bool>(
      "read_materials_from_xml", false,
      "Parameter determining whether user wishes to read materaial nuclide "
      "compositions from openMC XML file");

  params.addParam<FileName>(
      "materials_xml_file", "materials.xml",
      "Path and name of material xml file user wishes to use.");

  params.addParam<bool>(
      "write_photon_flux", false,
      "Boolean value used to determine whether to wite photon spectra to hdf5 "
      "after Fizzy has finished running.");

  params.addParam<FileName>(
      "photon_flux_filename", "photon_flux.h5",
      "Filename for the h5 file containing the output photon spectra");

  params.addParam<int>(
      "num_photon_bins", 24,
      "The number of bins to sort the output photon flux into");

  params.addParam<bool>(
      "comm_photon_flux", false,
      "Boolean value used to indicate whether to use boost::interprocess to "
      "communicate photon spectra through IPC.");

  params.addRequiredParam<FileName>(
      "molar_mass_data",
      "Filename for HDF5 file containing molar mass data for all isotopes");

  return params;
}

FispactProblem::FispactProblem(const InputParameters &params)
    : ExternalProblem(params), _fp_monitor(fispactLogName()),
      _fp_nuclear_data(_fp_monitor),
      _fp_nuclear_data_uo(getParam<std::string>("fispact_nuclear_data_uo")),
      _neutron_flux_filename(getParam<FileName>("neutron_flux_file")),
      _neutron_flux_tally_id(getParam<int>("neutron_flux_tally_id")),
      _photon_flux_filename(getParam<FileName>("photon_flux_filename")),
      _materials_from_xml(getParam<bool>("read_materials_from_xml")),
      _materials_xml_file(getParam<FileName>("materials_xml_file")),
      _neutron_bin_type(getParam<std::string>("neutron_bin_type")),
      _num_photon_bins(getParam<int>("num_photon_bins")),
      _schedule_uo_name(getParam<UserObjectName>("fispact_schedule_uo")),
      _local_domain_strength(0), _total_domain_strength(0),
      _write_photon_flux(getParam<bool>("write_photon_flux")),
      _comm_photon_flux(getParam<bool>("comm_photon_flux")),
      _interprocess_segment_name(generateInterprocessName()),
      _molar_mass_data_filename(getParam<FileName>("molar_mass_data")) {

  /**
   * Load materials from xml file if read_materials_from_xml is set to true,
   * and a materials xml filename has been passed
   */
  if (_materials_from_xml) {
    if (!isParamSetByUser("materials_xml_file")) {
      mooseWarning("read_materials_from_xml is set to true, but "
                   "materials_xml_file is not set! Defaulting to " +
                   _materials_xml_file);
    }
    /// Populate _mat_definitions with materials from openmc xml
    read_material_xml_data();
  }

  /**
   * If write_photon_flux was set to true then check that user input a
   * filename, if not use default
   */
  if (_write_photon_flux && !isParamSetByUser("photon_flux_filename")) {
    _console << "write_photon_flux is set to true but photon_flux_filename is "
                "not set! Photon flux filename defaulting to " +
                    _photon_flux_filename;
  }

  /// Set neutron bin type
  setNeutronBins();

  /// Check a corresponding material exists for all mesh blocks
  // checkMaterialsExist();

  /// Initialise FISPACT
  fp::GlobalInitialise(_fp_monitor);

  /// Load molar mass data
  loadMolarMasses();

  /// Read neutron flux from h5 file
  if (_comm_photon_flux) {

#ifdef LIBMESH_HAVE_BOOST

    /// Remove any shared memory region with a similar name just in case
    bi::shared_memory_object::remove(_interprocess_segment_name.c_str());

    /// Calculate space needed for shared mem region
    unsigned long shared_memory_size = calculateMemorySize();

    /// Create shared memory region
    _segment = bi::managed_shared_memory(bi::create_only,
                                         _interprocess_segment_name.c_str(),
                                         shared_memory_size);
#else
    mooseError("_comm_photon_flux is set to true but libmesh was not built "
               "with BOOST. No communication occuring.");
#endif
  }
}

void FispactProblem::initialSetup() {
  ExternalProblem::initialSetup();

  /// Set nuclear data paths
  setNuclearData(_fp_nuclear_data_uo);

  /// Get n FISPACT Schedule Times
  FispactSchedule &schedule = getUserObject<FispactSchedule>(_schedule_uo_name);

  /// Set _n_inventories
  _n_inventories = schedule.getTimes().size();

  /// Reserve space in our solution vector
  _photon_energy_spectra.resize(
      _mesh.nActiveLocalElem() * _n_inventories * _num_photon_bins, 0);

  /// Reserve space for element strengths vector
  _element_strengths.resize(_mesh.nActiveLocalElem() * _n_inventories, 0);

  _local_domain_strength.reserve(_n_inventories);

  _total_domain_strength.reserve(_n_inventories);

  /// Set up local element index map
  int local_elem_idx = 0;
  for (const libMesh::Elem *element : *_mesh.getActiveLocalElementRange()) {
    _local_elem_index.insert(
        std::pair<int, int>(element->id(), local_elem_idx++));
  }
}

void FispactProblem::timestepSetup() {

  ExternalProblem::timestepSetup();

  if (_comm_photon_flux) {
    if (timeStep() > 1) {

      _segment.destroy<PhotonSharingData>("photon_sharing_instance");
    }
  }
}

void FispactProblem::externalSolve() {

  if (!_solved) {

    /// Set up fispact input data
    fp::InputData fispact_input(_fp_monitor);
    fp::OutputData fispact_output(_fp_monitor);

    int counter = 1;
    for (const libMesh::Elem *element : *_mesh.getActiveLocalElementRange()) {

      /// Get element id
      dof_id_type elem_id = element->id();

      _console << std::endl
               << "Elem ID: " << std::to_string(elem_id) << std::endl;

      std::vector<double> neutron_flux = readElementNeutronFlux(
          _neutron_flux_filename, elem_id, _neutron_flux_tally_id);

      /// Output how many elements have been checked
      _console << counter++ << "/" << _mesh.getMesh().n_active_local_elem()
               << std::endl;

      /// If there is flux in the element, run FISPACT
      if (isFlux(neutron_flux)) {
        /// Here we are assuming the input mesh is in centimeters
        double element_volume = element->volume();

        /// Get element material definition
        const FISPACTMaterial &mat = getElementMaterial(elem_id);

        setFispactInputData(_fp_monitor, mat, neutron_flux, element_volume,
                            fispact_input);
        /// Run FISPACT
        fp::Process(fispact_input, _fp_nuclear_data, fispact_output,
                    _fp_monitor, process_callback);

        /// If photon bins aren't set yet, set them
        if (_photon_bins.empty()) {
          setPhotonBins(fispact_output);
        }

        /// Vector to store photon energy spectra in photons/cc-s
        std::vector<double> element_photon_energy_spectrum(_num_photon_bins, 0);

        /// Loop over number of inventories to get all data for current element
        for (int inv_index = 0; inv_index < _n_inventories; inv_index++) {

          convertGammaEvToCount(
              fispact_input, fispact_output.getGammaSpectrumBins(inv_index + 1),
              element_photon_energy_spectrum);

          _photon_energy_spectra.insert(
              _photon_energy_spectra.begin() +
                  photonEnergySpectraIdx(inv_index, elem_id),
              element_photon_energy_spectrum.begin(),
              element_photon_energy_spectrum.end());

          /**
           * Calculate photon source strength of this element and insert it into
           * _element_strengths
           */
          insertElementStrength(inv_index, element,
                                element_photon_energy_spectrum);
        }
      }
    }
    /// Calculate _local_domain_strength
    calculateLocalDomainStrength();

    if (_write_photon_flux) {
      const std::vector<double> &inv_times =
          getUserObject<FispactSchedule>(_schedule_uo_name)
              .getCumulativeTimes();
      writePhotonFlux(_photon_flux_filename, inv_times);
    }

    fp::GlobalFinalise(_fp_monitor);

    _solved = true;
  }
}

void FispactProblem::syncSolutions(ExternalProblem::Direction direction) {
  if (direction == ExternalProblem::Direction::FROM_EXTERNAL_APP) {

    if (_comm_photon_flux) {
#ifdef LIBMESH_HAVE_BOOST

      /// Find the the inventory index associated to time @"time()"
      const std::vector<double> &schedule_times =
          getUserObject<FispactSchedule>(_schedule_uo_name)
              .getCumulativeTimes();

      auto schedule_iterator =
          std::find(schedule_times.begin(), schedule_times.end(), time());
      if (schedule_iterator == schedule_times.end()) {
        mooseError("Current time " + std::to_string(time()) +
                   " does not match any entry in the FISPACT schedule. Cannot "
                   "do interprocess communication");
      }
      size_t inventory_idx =
          std::distance(schedule_times.begin(), schedule_iterator);

      /// Calculate TotalDomainStrength
      getTotalDomainStrength();

      std::vector<double> inv_photon_spectra = std::vector<double>(
          _photon_energy_spectra.begin() +
              (inventory_idx * _num_photon_bins * _mesh.nActiveLocalElem()),
          _photon_energy_spectra.begin() +
              (inventory_idx * _num_photon_bins * _mesh.nActiveLocalElem()) +
              (_mesh.nActiveLocalElem() * _num_photon_bins));

      std::vector<double> inv_element_strengths =
          std::vector<double>(_element_strengths.begin() +
                                  (inventory_idx * _mesh.nActiveLocalElem()),
                              _element_strengths.begin() +
                                  (inventory_idx * _mesh.nActiveLocalElem()) +
                                  _mesh.nActiveLocalElem());

      /// Create a shared instantiation of the photon sharing class
      PhotonSharingData *photon_sharing_instance =
          _segment.construct<PhotonSharingData>("photon_sharing_instance")(
              _segment, inv_photon_spectra, inv_element_strengths,
              _num_photon_bins, _mesh.getMesh().n_active_local_elem(),
              _total_domain_strength[inventory_idx],
              _local_domain_strength[inventory_idx], _photon_bins,
              _local_elem_index);
#else
      mooseError("_comm_photon_flux is set to true but libmesh was not
                 built with BOOST. No communication occuring.");
#endif
    }
  }
  if (direction == ExternalProblem::Direction::TO_EXTERNAL_APP) {
    if (_comm_photon_flux) {
#ifdef LIBMESH_HAVE_BOOST

#else

      mooseError("_comm_photon_flux is set to true but libmesh was not
                 built with BOOST. No communication occuring.");
#endif
    }
  }
}

void FispactProblem::writePhotonFlux(
    const std::string &filename, const std::vector<double> &inventory_times) {
/**
 * Do use MPI HDF5 driver, as we are calling H5Dcreate on each rank the same
 * number of times
 */
#ifdef H5_HAVE_PARALLEL
  bool parallel = true;
#else
  mooseWarning(
      "Writing photon flux requires HDF5 compiled with MPI. "
      "Simulation continuing, but photon spectra will not be written.");
  return;
#endif
  /// Open OpenMC statepoint file with Neutron Flux
  hid_t file_id =
      hdf5_utils::file_open(filename.c_str(), 'w', parallel, comm().get());

  for (int inv_id = 0; inv_id < _n_inventories; inv_id++) {

    std::string dataset_name =
        "photon_flux_" + std::to_string(inventory_times[inv_id]);

    int ndim = 2;
    unsigned long x_dim = _mesh.getMesh().n_active_elem();
    int y_dim = _num_photon_bins;
    /// Set up dimensions for photon flux dataspace
    hsize_t dataspace_dims[ndim];
    dataspace_dims[0] = x_dim;
    dataspace_dims[1] = y_dim;

    // Create filespace
    hid_t dataspace = H5Screate_simple(ndim, dataspace_dims, NULL);
    hid_t h5_dataset =
        H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_DOUBLE, dataspace,
                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(dataspace);

    for (const libMesh::Elem *element : *_mesh.getActiveLocalElementRange()) {

      /// Create hdf5 dataspace to be our memory space for writing
      hsize_t hyperslab_dims[2] = {1, static_cast<hsize_t>(y_dim)};

      /**
       * offset and count are used to select our hyperslab.
       * Given the dimensions of the dataspace are num_elems * 24, our offset
       * selection should be the element_id who's flux we wish to write
       */
      hsize_t start[2] = {element->id(), 0};

      /// count specifies the number of entries we wish to write in each
      /// dimension
      hsize_t count[2] = {1, _num_photon_bins};

      int flux_idx = photonEnergySpectraIdx(inv_id, element->id());

      hdf5_utils::write_double_hyperslab(
          h5_dataset, nullptr, ndim, hyperslab_dims, start, count,
          &_photon_energy_spectra[flux_idx], parallel);
    }
    /// Close all the HDF5 bits and pieces
    H5Dclose(h5_dataset);
  }

  /// Write bins to hdf5 file as well
  // writePhotonFluxBins(file_id, parallel);

  H5Fclose(file_id);
}

void FispactProblem::writePhotonFluxBins(const hid_t &file_id,
                                         const bool parallel) {
  /// Check _photon_bins are actually set
  std::vector<double> &photon_bins = getPhotonBins();

  /// Set up hsize_t object to hold dataset dimensions
  int ndim = 1;
  hsize_t bin_dataset_dims[ndim];
  bin_dataset_dims[0] = photon_bins.size();

  std::string dataset_name = "photon_bins";

  /// Create a dataspace for the photon bins
  hid_t filespace_photon_bins = H5Screate_simple(1, bin_dataset_dims, NULL);

  /// Create a dataset for the photon bins using the dataspace
  hid_t h5_dataset_photon_bins =
      H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_DOUBLE,
                filespace_photon_bins, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dclose(h5_dataset_photon_bins);

  hdf5_utils::write_dataset_lowlevel(file_id, dataset_name.c_str(), ndim,
                                     bin_dataset_dims, H5T_NATIVE_DOUBLE,
                                     photon_bins.data(), parallel);
}

void FispactProblem::setPhotonBins(const fp::OutputData &fispact_output) {

  // Retrieve gamma spectrum boundaries from first fispact inv step
  size_t fispact_step = 1;
  std::vector<double> bounds =
      fispact_output.getGammaSpectrumBoundaries(fispact_step);

  _photon_bins = bounds;

  // Scale photon bin entries to get them into eV
  for (int i = 0; i < _photon_bins.size(); i++) {
    if (i == 0) {
      continue;
    }
    _photon_bins[i] *= 1e6;
  }

  // Broadcast to all processors. Prevents issues resulting from all the
  // elements on one processor having zero flux, and hence _photon_bins never
  // getting set
  comm().broadcast(_photon_bins);
}

std::vector<double>
FispactProblem::readElementNeutronFlux(const std::string &filename,
                                       const dof_id_type &elem_id,
                                       const int &tally_id) const {
  /**
   * Don't do this with HDF5 MPI driver, as we are going to call H5Dopen a
   * different number of times on each rank!
   */
  bool parallel = false;
  /// Open OpenMC statepoint file with Neutron Flux
  hid_t file_id =
      hdf5_utils::file_open(filename.c_str(), 'r', parallel, comm().get());

  /// Open up tallies group in HDF5 file
  hid_t tallies_group_id = hdf5_utils::open_group(file_id, "tallies");

  /// Open up the neutron flux tally group within the tallies group
  std::string tally_group_name = "tally " + std::to_string(tally_id);
  hid_t neutron_tally_group_id =
      hdf5_utils::open_group(tallies_group_id, tally_group_name);

  /// Get number of tally realizations
  int n_realizations;
  hdf5_utils::read_int(neutron_tally_group_id, "n_realizations",
                       &n_realizations, parallel, true);

  /// Define size of selection of neutron flux
  hsize_t dataset_dims[3];
  hsize_t hyperslab_dims[3]{static_cast<hsize_t>(_num_neutron_bins), 1, 1};
  hsize_t start[]{elem_id * _num_neutron_bins, 0, 0};
  hsize_t count[]{static_cast<hsize_t>(_num_neutron_bins), 1, 1};

  // Check this statepoint file is of
  hid_t flux_dataset =
      hdf5_utils::open_dataset(neutron_tally_group_id, "results");
  hdf5_utils::get_shape(flux_dataset, dataset_dims);

  if (dataset_dims[0] != _num_neutron_bins * _mesh.getMesh().n_elem()) {
    mooseError("Provided tally dimensions do not conform to this mesh and "
               "neutron bin structure. Tally dimensions are (" +
               std::to_string(dataset_dims[0]) + "," +
               std::to_string(dataset_dims[1]) + "," +
               std::to_string(dataset_dims[2]) + ")");
  }

  std::vector<double> neutron_flux_results(_num_neutron_bins, 0);

  /// Number of dimensions of hyperslab to read
  hsize_t rank = 3;
  hdf5_utils::read_double_hyperslab(neutron_tally_group_id, "results", rank,
                                    hyperslab_dims, start, count,
                                    neutron_flux_results.data(), parallel);

  hdf5_utils::close_group(tallies_group_id);
  hdf5_utils::close_group(neutron_tally_group_id);

  hdf5_utils::file_close(file_id);

  /**
   * Divide every value in neutron flux results vector by n_realizations to get
   * the mean value
   */
  for (auto &bin : neutron_flux_results) {
    bin /= n_realizations;
  }

  return neutron_flux_results;
}

/**
 * TODO: Break up this function into setFispactInputFlux and
 * setFispactInputMaterial
 */
void FispactProblem::setFispactInputData(
    const fp::FispactMonitor &monitor, const FISPACTMaterial &material,
    const std::vector<double> &neutron_flux, const double &volume,
    fp::InputData &input) const {

  /// Set neutron flux
  input.setFlux(_neutron_bins, neutron_flux);
  input.setFluxWallLoading(1.0);
  input.setFluxName("neutrons");

  /// Get density from mat density, in g/cm^3!
  double density = material.getDensity();
  input.setDensity(density);

  /// Set atoms threshold
  input.setAtomsThreshold(1.0e3);

  /// Volume read in is in cm^3, so we need to scale by 1e-6, as mass is in kg
  double total_mass = density * volume * 1e-6;

  //
  if (material.getMaterialType() == "MASS") {

    std::vector<int> atomic_numbers;
    std::vector<double> percent;

    /// Set total mass
    input.setMassTotal(total_mass);

    const std::unordered_map<std::string, double> &nuclideFractionMap =
        material.getNuclideFractionMap();

    atomic_numbers.reserve(nuclideFractionMap.size());

    for (auto &[element_name, mass_fraction] : nuclideFractionMap) {

      atomic_numbers.push_back(
          fp::util::GetAtomicNumberFromElementName(monitor, element_name));
    }

    input.setMass(atomic_numbers, material.getNuclideFractions());

  } else if (material.getMaterialType() == "FUEL") {

    /// Get material map, that maps from map[nuclide_name] -> mass_fraction
    const std::unordered_map<std::string, double> &nuclideFractionMap =
        material.getNuclideFractionMap();

    std::vector<int> zais;
    zais.reserve(nuclideFractionMap.size());
    std::vector<double> atoms;
    atoms.reserve(nuclideFractionMap.size());

    /// For all key (isotope name) value (mass_fraction) pairs in map, calculate
    /// the number of atoms pertaining to each isotope and append to input fuel
    for (const auto &[isotope_name, mass_fraction] : nuclideFractionMap) {

      double zai_mass = total_mass * mass_fraction;

      zais.push_back(fp::util::GetZai(monitor, isotope_name));

      atoms.push_back(
          getNumAtoms(zai_mass, _molar_mass_map.at(zais.back()), AVOGADRO));
    }
    input.setFuel(zais, atoms);
  }

  /// Get the fispact input schdule from the user object and set it in the
  /// FISPACT input
  setFispactSchedule(input);
}

void FispactProblem::setFispactSchedule(fp::InputData &input) const {

  FispactSchedule &schedule = getUserObject<FispactSchedule>(_schedule_uo_name);

  const std::vector<double> &flux_schedule = schedule.getFluxSchedule();
  const std::vector<double> &times = schedule.getTimes();
  input.setSchedule(times, flux_schedule);
}

const FISPACTMaterial &
FispactProblem::getElementMaterial(dof_id_type &elem_id) {

  libMesh::Elem *elem = _mesh.elemPtr(elem_id);

  // Query the warehouse to see if a FISPACTMaterial exists on the block this
  // element is assigned to
  std::vector<GeneralUserObject *> objs;
  theWarehouse()
      .query()
      .condition<AttribSystem>("UserObject")
      .condition<AttribSubdomains>(elem->subdomain_id())
      .queryInto(objs);

  // Remove extraneous user objects that are not FISPACTMaterials. Having done
  // this, only one object should remain in the vector, and it should be the
  for (auto it = objs.begin(); it != objs.end();) {
    if ((*it)->type() != "FISPACTMaterial") {
      it = objs.erase(it);
    } else {
      ++it;
    }
  }

  // FISPACT material pertaining to this block. If there ismore than one object,
  // then two FISPACTMaterials are assigned to this block, and that makes no
  // blimmin sense does it
  if (objs.empty()) {
    mooseError("Unable to find FISPACTMaterial object on block " +
               std::to_string(elem->subdomain_id()));
  } else if (objs.size() > 1) {
    mooseError("More than 1 FISPACTMaterial definition exists on block " +
               std::to_string(elem->subdomain_id()));
  }

  /// Return the FISPACTMaterial
  return *(static_cast<FISPACTMaterial *>(objs[0]));
}

bool FispactProblem::isFlux(const std::vector<double> &flux) const {
  /// Check if there is any neutron flux in current element
  bool is_zero_flux =
      std::all_of(flux.begin(), flux.end(), [](double j) { return j == 0; });

  /**
   * return true if there is flux, false if there isn't, as the function name
   * implies
   */
  return !is_zero_flux;
}

void FispactProblem::read_material_xml_data() {
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(_materials_xml_file.c_str());

  if (!result) {
    mooseError("No file called " + _materials_xml_file +
               " could be found, exiting.");
  }
  for (pugi::xml_node material : doc.child("materials").children()) {

    /// Get material name
    std::string material_name = material.attribute("name").value();
    /// Get material density
    double density =
        std::stod(material.child("density").attribute("value").value());

    /// Get atomic composition of material
    std::vector<std::pair<std::string, double>> atomic_comp;

    for (pugi::xml_node nuclide : material.children("nuclide")) {
      std::pair<std::string, double> nuclide_symbol_and_percentage =
          std::make_pair(std::string(nuclide.attribute("name").value()),
                         std::stod(nuclide.attribute("ao").value()));
      atomic_comp.push_back(nuclide_symbol_and_percentage);
    }

    /// Create material definition
    MaterialDefinition material_def{material_name, atomic_comp, density};
    /// Insert material definition into material map
    _mat_definitions.insert(std::make_pair(material_name, material_def));
  }
}

void FispactProblem::setNeutronBins() {
  if (_neutron_bin_type == "G1102") {
    _neutron_bins = fp::groups::G1102();
    _num_neutron_bins = 1102;
  }
  /// TODO: Add more bins types
}

void FispactProblem::convertGammaEvToCount(
    const fp::InputData &input,
    const std::vector<double> &photon_energy_spectra_ev,
    std::vector<double> &photon_energy_spectra_per_cc_s) {

  /// Reserve memory for photons per cc per s vector
  photon_energy_spectra_per_cc_s.resize(photon_energy_spectra_ev.size());

  /// Get inventory density and mass
  double inv_density = input.getDensity();
  double inv_mass = input.getMassTotal();

  const std::vector<double> &photon_flux_bins = getPhotonBins();
  /// Convert from MeV/s to per cc per s for each bin
  for (int i = 0; i < photon_energy_spectra_ev.size(); i++) {
    double bin_energy = photon_flux_bins[i] +
                        ((photon_flux_bins[i + 1] - photon_flux_bins[i]) / 2);

    /**
     * per_cc_per_s = MeV/s * (inventory_density/(inventory_mass *
     * energy_bin_midpoint))
     */
    double per_cc_per_s =
        photon_energy_spectra_ev[i] * (inv_density / (inv_mass * bin_energy));

    photon_energy_spectra_per_cc_s[i] = per_cc_per_s;
  }
}

double FispactProblem::calculateElementStrength(
    const libMesh::Elem *element, const std::vector<double> element_flux) {
  double element_strength = 0;

  for (double flux_bin : element_flux) {
    element_strength += flux_bin * element->volume();
  }

  return element_strength;
}

void FispactProblem::insertElementStrength(
    const int inv_index, const libMesh::Elem *element,
    const std::vector<double> element_flux) {
  double elem_strength = calculateElementStrength(element, element_flux);

  int idx =
      (inv_index * _mesh.nActiveLocalElem()) + _local_elem_index[element->id()];
  _element_strengths.insert(_element_strengths.begin() + idx, elem_strength);
}

void FispactProblem::calculateLocalDomainStrength() {

  for (int inv_index = 0; inv_index < _n_inventories; inv_index++) {
    for (const libMesh::Elem *element : *_mesh.getActiveLocalElementRange()) {

      double element_strength;

      /**
       * Attempt to retrieve element strength from map, otherwise catch
       * exception and give a useful error
       */
      int idx = (inv_index * _mesh.nActiveLocalElem()) +
                _local_elem_index[element->id()];

      element_strength = _element_strengths.at(idx);

      _local_domain_strength[inv_index] += element_strength;
    }
  }
}

void FispactProblem::getTotalDomainStrength() {

  _total_domain_strength = _local_domain_strength;
  comm().sum(_total_domain_strength);
}

int FispactProblem::calculateMemorySize() {
  /// Get number of active local elements
  dof_id_type n_local_elem = _mesh.getMesh().n_active_local_elem();

  size_t photon_flux_map_size =
      ((_num_photon_bins * sizeof(double)) + sizeof(int)) * n_local_elem;

  size_t element_strengths_map_size =
      (sizeof(int) + sizeof(double)) * n_local_elem;

  size_t memory_size = photon_flux_map_size + element_strengths_map_size +
                       (sizeof(int) * 2) + (sizeof(double) * 2);
  /**
   * Really naive way of doing this, but currently giving a 20% buffer to
   * account for the memory space required by Boost allocators and such
   */
  return memory_size * 2;
}

std::string FispactProblem::fispactLogName() {
  std::string log_name = "FISPACT_app_" +
                         this->getMooseApp().getInputFileNames()[0] +
                         std::to_string(processor_id()) + ".log";
  return log_name;
}

void FispactProblem::setNuclearData(const std::string &fp_nuclear_data_uo) {

  FispactNuclearDataPaths &nuclear_data_paths =
      getUserObject<FispactNuclearDataPaths>(fp_nuclear_data_uo);

  nuclear_data_paths.loadNuclearData(_fp_nuclear_data, _fp_monitor,
                                     load_callback);
  // nd_reader.setPath(FISPACT_ND_IND_NUC_KEY,
  //                   nd_base_path + "decay2020/decay_2020_index.txt");
  //
  // nd_reader.setPath(FISPACT_ND_XS_ENDF_KEY,
  //                   nd_base_path + "TENDL2021data/gendf-1102");
  // nd_reader.setPath(FISPACT_ND_PROB_TAB_KEY,
  //                   nd_base_path + "TENDL2021data/tp-1102-294");
  //
  // nd_reader.setPath(FISPACT_ND_FY_ENDF_KEY,
  //                   nd_base_path + "GEFY61data/gefy61_nfy");
  // nd_reader.setPath(FISPACT_ND_SF_ENDF_KEY,
  //                   nd_base_path + "GEFY61data/gefy61_sfy");
  //
  // nd_reader.setPath(FISPACT_ND_DK_ENDF_KEY,
  //                   nd_base_path + "decay2020/decay_2020");
  //
  // nd_reader.setPath(FISPACT_ND_ABSORP_KEY, nd_base_path +
  // "decay2012/abs_2012");
  //
  // nd_reader.setPath(FISPACT_ND_HAZARDS_KEY,
  //                   nd_base_path + "/decay2012/hazards_2012");
  // nd_reader.setPath(FISPACT_ND_CLEAR_KEY,
  //                   nd_base_path + "/decay2012/clear_2012");
  // nd_reader.setPath(FISPACT_ND_A2DATA_KEY, nd_base_path +
  // "/decay2012/a2_2012");
  //
  // nd_reader.load(_fp_nuclear_data, &FispactProblem::load_callback);
}

const std::string FispactProblem::generateInterprocessName() {
  char mpi_proc_name[MPI_MAX_PROCESSOR_NAME];
  int len = 0;
  int err = MPI_Get_processor_name(mpi_proc_name, &len);

  std::string ipc_name = std::string(mpi_proc_name);
  ipc_name += "_" + std::to_string(comm().rank());

  return ipc_name;
}

double FispactProblem::getNumAtoms(const double &mass, const double &molar_mass,
                                   const double &avogadro) const {
  return (mass / molar_mass) * avogadro;
}

void FispactProblem::loadMolarMasses() {

  hid_t molar_mass_file = hdf5_utils::file_open(
      _molar_mass_data_filename.c_str(), 'r', false, comm().get());

  hid_t molar_mass_dataset =
      hdf5_utils::open_dataset(molar_mass_file, "MolarMass");

  hid_t element_symbols_dataset =
      hdf5_utils::open_dataset(molar_mass_file, "symbol");

  hsize_t molar_mass_dims[MOLAR_MASS_DATASET_DIMS];

  hdf5_utils::get_shape(molar_mass_dataset, molar_mass_dims);

  // Vectors to store molar masses and respective element symbols
  std::vector<double> molar_masses(molar_mass_dims[0]);
  std::vector<std::string> element_symbols(molar_mass_dims[0]);

  hdf5_utils::read_double(molar_mass_dataset, nullptr, molar_masses.data(),
                          false);

  hdf5_utils::read_string(element_symbols_dataset, nullptr, element_symbols, 8,
                          false);

  for (int i = 0; i < molar_masses.size(); i++) {

    int zai = fp::util::GetZai(_fp_monitor, element_symbols[i]);
    std::pair<int, double> key_value =
        std::pair<int, double>(zai, molar_masses[i]);
    _molar_mass_map.insert(key_value);
  }
}
