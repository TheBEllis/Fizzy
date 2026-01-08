#include "ExternalProblem.h"
#include "FispactProblem.h"

/// Custom user object includes
#include "FispactSchedule.h"

#include "HDF5Utils.h"
#include "MooseError.h"
#include "MooseTypes.h"

/// Fispact includes

//// PugiXML include
#include "fispactutil.hpp"
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

  params.addRequiredParam<FileName>("fispact_nuclear_data_path",
                                    "Path to FISPACT nuclear data");

  params.addRequiredParam<std::string>("neutron_bin_type",
                                       "neutron binning scheme for FISPACT");

  params.addRequiredParam<std::string>(
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
      _fp_nuclear_data_path(getParam<FileName>("fispact_nuclear_data_path")),
      _neutron_flux_filename(getParam<FileName>("neutron_flux_file")),
      _neutron_flux_tally_id(getParam<int>("neutron_flux_tally_id")),
      _photon_flux_filename(getParam<FileName>("photon_flux_filename")),
      _materials_from_xml(getParam<bool>("read_materials_from_xml")),
      _materials_xml_file(getParam<FileName>("materials_xml_file")),
      _neutron_bin_type(getParam<std::string>("neutron_bin_type")),
      _num_photon_bins(getParam<int>("num_photon_bins")),
      _schedule_uo_name(getParam<std::string>("fispact_schedule_uo")),
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

  /// Set nuclear data paths
  setNuclearData(_fp_nuclear_data_path);

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

void FispactProblem::syncSolutions(ExternalProblem::Direction direction) {
  if (direction == ExternalProblem::Direction::FROM_EXTERNAL_APP) {

    if (_comm_photon_flux) {
#ifdef LIBMESH_HAVE_BOOST

      /// Create vector to store all local domain strengths
      std::vector<double> local_domain_strengths;

      /// Gather local domain strengths from all processors
      comm().allgather(_local_domain_strength, local_domain_strengths);

      /// Calculate TotalDomainStrength
      getTotalDomainStrength();

      /// Create a shared instantiation of the photon sharing class
      PhotonSharingData *photon_sharing_instance =
          _segment.construct<PhotonSharingData>(
              "PhotonSharingData photon_sharing_instance")(
              _segment, _photon_fluxes, _element_strengths, _num_photon_bins,
              (int)_mesh.getMesh().n_active_local_elem(),
              _total_domain_strength, _local_domain_strength,
              local_domain_strengths, _photon_bins);
#else
      mooseError("_comm_photon_flux is set to true but libmesh was not built "
                 "with BOOST. No communication occuring.");
#endif
    }
  }
}

void FispactProblem::externalSolve() {

  /// Set up fispact input data
  fp::InputData fispact_input(_fp_monitor);
  fp::OutputData fispact_output(_fp_monitor);

  int counter = 1;
  for (const libMesh::Elem *element : *_mesh.getActiveLocalElementRange()) {

    /// Get element id
    dof_id_type elem_id = element->id();

    _console << "Elem ID: " << std::to_string(elem_id) << std::endl;

    std::vector<double> neutron_flux = readElementNeutronFlux(
        _neutron_flux_filename, elem_id, _neutron_flux_tally_id);

    /// Output how many elements have been checked
    _console << counter++ << "/" << _mesh.getMesh().n_active_local_elem()
             << std::endl;

    std::vector<double> photon_spectra(_num_photon_bins, 0);

    /// Check if there is any neutron flux in current element
    bool is_flux = isFlux(neutron_flux);

    /// If there is flux in the element, run FISPACT
    if (is_flux) {
      /// Here we are assuming the input mesh is in centimeters
      double element_volume = element->volume();

      /// Get element material definition
      const FISPACTMaterial &mat = getElementMaterial(elem_id);

      setFispactInputData(_fp_monitor, mat, neutron_flux, element_volume,
                          fispact_input);
      /// Run FISPACT
      fp::Process(fispact_input, _fp_nuclear_data, fispact_output, _fp_monitor,
                  process_callback);

      /// If photon bins aren't set yet, set them
      if (_photon_bins.empty()) {
        setPhotonBins(fispact_output);
      }

      convertGammaEvToCount(fispact_input,
                            fispact_output.getGammaSpectrumBins(1),
                            photon_spectra);
    }

    _photon_fluxes.insert(
        std::pair<int, std::vector<double>>(elem_id, photon_spectra));
    /**
     * Calculate photon source strength of this element and insert it into
     * _element_strengths map
     */
    insertElementStrength(element, _photon_fluxes.at(elem_id));

    /// Add calculated element strength to _local_domain_strength
    updateLocalDomainStrength(element);
  }

  if (_write_photon_flux) {
    writePhotonFlux(_photon_flux_filename);
  }

  fp::GlobalFinalise(_fp_monitor);
}

void FispactProblem::writePhotonFlux(const std::string &filename) {
/**
 * Do use MPI HDF5 driver, as we are calling H5Dcreate on each rank the same
 * number of times
 */
#ifdef H5_HAVE_PARALLEL
  bool parallel = true;
#elif
  parallel = false;
#endif
  /// Open OpenMC statepoint file with Neutron Flux
  hid_t file_id =
      hdf5_utils::file_open(filename.c_str(), 'w', parallel, comm().get());

  std::string dataset_name = "photon_flux";
  int ndim = 2;
  int x_dim = _mesh.getMesh().n_active_elem();

  int y_dim = _num_photon_bins;
  hsize_t dataspace_dims[ndim];

  /// Set up dimensions for photon flux dataspace
  dataspace_dims[0] = x_dim;
  dataspace_dims[1] = y_dim;

  hid_t filespace = H5Screate_simple(ndim, dataspace_dims, NULL);
  hid_t h5_dataset =
      H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_DOUBLE, filespace,
                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Sclose(filespace);

  for (auto &element_flux_pair : _photon_fluxes) {
    /**
     * CONTEXT: A hyperslab is a section of a hdf5 dataspace, don't be put off
     * by the fancy name
     */

    /// Create hdf5 dataspace to be our memory space for writing
    hsize_t hyperslab_dims[2] = {1, y_dim};

    /// Retrieve element id from _photon_fluxes map
    hsize_t element_id = (hsize_t)element_flux_pair.first;

    /**
     * offset and count are used to select our hyperslab.
     * Given the dimensions of the dataspace are num_elems * 24, our offset
     * selection should be the element_id who's flux we wish to write
     */
    hsize_t start[2] = {element_id, 0};

    /// count specifies the number of entries we wish to write in each dimension
    hsize_t count[2] = {1, _num_photon_bins};

    hdf5_utils::write_double_hyperslab(h5_dataset, nullptr, ndim,
                                       hyperslab_dims, start, count,
                                       element_flux_pair.second.data(), true);
  }
  /// Close all the HDF5 bits and pieces
  H5Dclose(h5_dataset);

  /// Write bins to hdf5 file as well
  writePhotonFluxBins(file_id, parallel);

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
  int fispact_step = 1;
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

  /// Define size of selection of neutron flux
  hsize_t dims[3]{static_cast<hsize_t>(_num_neutron_bins), 1, 1};
  hsize_t start[]{elem_id * _num_neutron_bins, 0, 0};
  hsize_t count[]{static_cast<hsize_t>(_num_neutron_bins), 1, 1};

  /// Get number of tally realizations
  int n_realizations;
  hdf5_utils::read_int(neutron_tally_group_id, "n_realizations",
                       &n_realizations, parallel, true);

  std::vector<double> neutron_flux_results(_num_neutron_bins, 0);

  /// Number of dimensions of hyperslab to read
  hsize_t rank = 3;
  hdf5_utils::read_double_hyperslab(neutron_tally_group_id, "results", rank,
                                    dims, start, count,
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

    std::vector<int> atomic_numbers;
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

    /// For all key (isotope name) value (mass_fraction) pairs in map, calculate
    /// the number of atoms pertaining to each isotope and append to input fuel
    for (const auto &[isotope_name, mass_fraction] : nuclideFractionMap) {

      double zai_mass = total_mass * mass_fraction;

      double zai = fp::util::GetZai(monitor, isotope_name);

      double atoms = getNumAtoms(zai_mass, _molar_mass_map.at(zai), AVOGADRO);

      input.appendFuel(zai, atoms);
    }
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
    const fp::InputData &input, const std::vector<double> &photon_spectra,
    std::vector<double> &photons_per_cc_per_s) {

  /// Reserve memory for photons per cc per s vector
  photons_per_cc_per_s.resize(photon_spectra.size());

  /// Get inventory density and mass
  double inv_density = input.getDensity();
  double inv_mass = input.getMassTotal();

  const std::vector<double> &photon_flux_bins = getPhotonBins();
  /// Convert from MeV/s to per cc per s for each bin
  for (int i = 0; i < photon_spectra.size(); i++) {
    double bin_energy = photon_flux_bins[i] +
                        ((photon_flux_bins[i + 1] - photon_flux_bins[i]) / 2);

    /**
     * per_cc_per_s = MeV/s * (inventory_density/(inventory_mass *
     * energy_bin_midpoint))
     */
    double per_cc_per_s =
        photon_spectra[i] * (inv_density / (inv_mass * bin_energy));

    photons_per_cc_per_s[i] = per_cc_per_s;
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
    const libMesh::Elem *element, const std::vector<double> element_flux) {
  double elem_strength = calculateElementStrength(element, element_flux);

  _element_strengths.insert(
      std::pair<int, double>(element->id(), elem_strength));
}

void FispactProblem::updateLocalDomainStrength(const libMesh::Elem *element) {
  double element_strength;

  /**
   * Attempt to retrieve element strength from map, otherwise catch exception
   * and give a useful error
   */
  try {
    element_strength = _element_strengths.at(element->id());
  } catch (std::out_of_range) {
    mooseError("Attempted to access element strength for an element ID that "
               "has not had strength calcalculated");
  }

  _local_domain_strength += element_strength;
}

void FispactProblem::getTotalDomainStrength() {

  _total_domain_strength = _local_domain_strength;
  comm().sum(_total_domain_strength);
}

int FispactProblem::calculateMemorySize() {
  /// Get number of active local elements
  int n_local_elem = _mesh.getMesh().n_active_local_elem();

  unsigned long photon_flux_map_size =
      ((_num_photon_bins * sizeof(double)) + sizeof(int)) * n_local_elem;

  unsigned long element_strengths_map_size =
      (sizeof(int) + sizeof(double)) * n_local_elem;

  unsigned long memory_size = photon_flux_map_size +
                              element_strengths_map_size + (sizeof(int) * 2) +
                              (sizeof(double) * 2);
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

void FispactProblem::setNuclearData(const std::string &nd_base_path) {
  fp::io::NuclearDataReader nd_reader(_fp_monitor);

  nd_reader.setPath(FISPACT_ND_IND_NUC_KEY,
                    nd_base_path + "decay2020/decay_2020_index.txt");

  nd_reader.setPath(FISPACT_ND_XS_ENDF_KEY,
                    nd_base_path + "TENDL2021data/gendf-1102");
  nd_reader.setPath(FISPACT_ND_PROB_TAB_KEY,
                    nd_base_path + "TENDL2021data/tp-1102-294");

  nd_reader.setPath(FISPACT_ND_FY_ENDF_KEY,
                    nd_base_path + "GEFY61data/gefy61_nfy");
  nd_reader.setPath(FISPACT_ND_SF_ENDF_KEY,
                    nd_base_path + "GEFY61data/gefy61_sfy");

  nd_reader.setPath(FISPACT_ND_DK_ENDF_KEY,
                    nd_base_path + "decay2020/decay_2020");
  nd_reader.setPath(FISPACT_ND_ABSORP_KEY, nd_base_path + "decay/abs_2012");

  nd_reader.setPath(FISPACT_ND_HAZARDS_KEY,
                    nd_base_path + "/decay/hazards_2012");
  nd_reader.setPath(FISPACT_ND_CLEAR_KEY, nd_base_path + "/decay/clear_2012");
  nd_reader.setPath(FISPACT_ND_A2DATA_KEY, nd_base_path + "/decay/a2_2012");

  nd_reader.load(_fp_nuclear_data, &FispactProblem::load_callback);
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
