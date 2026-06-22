#include "EnergyGroups.h"
#include "FispactInventoryManager.h"
#include "FispactProblem.h"
//// PugiXML include
#include "FizzyEnums.h"
#include "pugixml.hpp"

/// Cpp includes
#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <memory>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <mpi.h>

// Avogadro's number
constexpr double AVOGADRO = 6.02214076e+23;

#define MOLAR_MASS_DATASET_DIMS 1

registerMooseObject("FizzyApp", FispactProblem);

InputParameters FispactProblem::validParams() {
  InputParameters params = ExternalProblem::validParams();

  params.addClassDescription(
      "Run a FISPACT-II solve using MOOSE input syntax, and obtain inventory "
      "quantities as aux variables.");

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

  params.addRequiredParam<UserObjectName>(
      "fispact_nuclear_data_uo", "Name of the FispactNuclearDataPaths objects "
                                 "to use for setting nuclear data");

  params.addRequiredParam<UserObjectName>(
      "fispact_schedule_uo", "Name of the FispactSchedule user object defining "
                             "the FISPACT flux schedule");

  params.addRequiredParam<UserObjectName>(
      "fispact_input_flux_uo",
      "Name of the FispactFluxInput user object defining "
      "the input flux spectra");

  params.addParam<bool>("mock_fispact", false,
                        "Parameter exclusively for unit testing when we don't "
                        "wish to use the actual FISPACT library");
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

  params.addParam<bool>("exclude_xrays", false,
                        "Boolean value used to determine whether to exclude "
                        "xrays from output gamma spectrum.");

  params.addParam<FileName>(
      "photon_flux_filename", "photon_flux.h5",
      "Filename for the h5 file containing the output photon spectra");

  params.addParam<uint64_t>(
      "num_photon_bins", 24,
      "The number of bins to sort the output photon flux into");

  params.addParam<bool>(
      "comm_photon_flux", false,
      "Boolean value used to indicate whether to use boost::interprocess to "
      "communicate photon spectra through IPC.");

  params.addRequiredParam<FileName>("molar_mass_data",
                                    "Filename for HDF5 file containing molar "
                                    "mass data in g/mol for all isotopes");

  params.addParam<MooseEnum>(
      "conversion_type", MooseEnum("LETHARGY ENERGY", "LETHARGY"),
      "Setting to determine whether to convert by lethargy or energy");

  params.addParam<double>("output_inventory_time",
                          "When using the Steady executioner, which inventory "
                          "step should be placed in interprocess data");

  params.addParam<double>("rtol", 2e-3, "Relative FISPACT solver tolerance.");
  params.addParam<double>("atol", 1e4, "Absolute FISPACT solver tolerance.");

  params.addParam<bool>(
      "uniform_sampling", false,
      "When using distributed sampling, this forced all elements to be sampled "
      "uniformly, with their weighting changed to account for their emission "
      "strengths.");
  return params;
}

FispactProblem::FispactProblem(const InputParameters &params)
    : ExternalProblem(params),
      _photon_flux_filename(getParam<FileName>("photon_flux_filename")),
      _materials_from_xml(getParam<bool>("read_materials_from_xml")),
      _materials_xml_file(getParam<FileName>("materials_xml_file")),
      _n_photon_bins(getParam<uint64_t>("num_photon_bins")),
      _fp_schedule_uo_name(getParam<UserObjectName>("fispact_schedule_uo")),
      _fp_flux_uo_name(getParam<UserObjectName>("fispact_input_flux_uo")),
      _fp_nuclear_data_uo_name(
          getParam<UserObjectName>("fispact_nuclear_data_uo")),
      _local_domain_strength(0), _total_domain_strength(0),
      _write_photon_flux(getParam<bool>("write_photon_flux")),
      _comm_photon_flux(getParam<bool>("comm_photon_flux")), _solved(false),
      _uniform(getParam<bool>("uniform_sampling")),
      _interprocess_segment_name(generateInterprocessName()),
      _molar_mass_data_filename(getParam<FileName>("molar_mass_data")),
      _atol(getParam<double>("atol")), _rtol(getParam<double>("rtol")),
      _exclude_xrays(getParam<bool>("exclude_xrays")) {

  /**
   * If write_photon_flux was set to true then check that user input a
   * filename, if not use default
   */
  if (_write_photon_flux && !isParamSetByUser("photon_flux_filename")) {
    paramWarning("write_photon_flux",
                 "write_photon_flux is set to true but photon_flux_filename is "
                 "not set! Photon flux filename defaulting to " +
                     _photon_flux_filename);
  }

  /// Set up local element index map
  int local_elem_idx = 0;
  for (const libMesh::Elem *element : *_mesh.getActiveLocalElementRange()) {
    _local_elem_index.insert(
        std::pair<int, int>(element->id(), local_elem_idx++));
  }

  if (_comm_photon_flux) {

#ifdef LIBMESH_HAVE_BOOST

    /// Remove any potentially left over shared memory segments
    bi::shared_memory_object::remove(_interprocess_segment_name.c_str());
    comm().barrier();

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

  /// Need to load molar masses here, as FispactMaterial will need molar mass
  /// data in its constructor
  loadMolarMasses();
}

FispactProblem::~FispactProblem() {
  /// Remove any shared memory region with a similar name just in case
  // bi::shared_memory_object::remove(_interprocess_segment_name.c_str());
}

void FispactProblem::initialSetup() {
  ExternalProblem::initialSetup();

  callFispactFactory();

  if (!_fp_ctxt) {
    mooseError("Fispact context could not be initialised! Maybe user has asked "
               "to use mock fispact outside of unit tests?");
  }
  // Init FISPACT
  _fp_ctxt->globalInitialise();

  resolveFispactUserObjects();

  setNuclearData();

  checkForEnergyGroupConsistency();

  setPhotonBins(_fp_ctxt->getUtils().getPhotonEnergyBounds(_n_photon_bins));

  _n_inventories = &(_fp_schedule_uo->getNumInventories());

  _n_solution_inventories = &(_fp_schedule_uo->getNumSolutionInventories());
  ///

  /// Reserve space in our solution vector
  // _photon_energy_spectra.resize(
  //     _mesh.nActiveLocalElem() * *_n_solution_inventories * _n_photon_bins,
  //     0);

  _photon_energy_spectra = std::make_unique<PhotonSpectra>(
      *_n_solution_inventories, _mesh.nActiveLocalElem(), _n_photon_bins);

  /// Reserve space for element strengths vector
  _element_strengths.resize(
      _mesh.nActiveLocalElem() * (*_n_solution_inventories), 0);

  _local_domain_strength.resize(*_n_solution_inventories, 0);

  _total_domain_strength.reserve(*_n_solution_inventories);

  /// Check user has passed output_inventory_time, if problem is Steady and
  /// they wish to use distributed sampling
  if (_comm_photon_flux) {
    if (!isParamSetByUser("output_inventory_time") && !isTransient()) {
      paramError("output_inventory_time",
                 "Parameter not set! When using a Steady executioner and "
                 "comm_photon_flux, user "
                 "must provide the inventory time to be communicated.");
    } else if (isParamSetByUser("output_inventory_time") && isTransient()) {
      paramWarning("output_inventory_time",
                   "Parameter is set, but executioner is Transient. Ignoring "
                   "parameter.");
    }
  }

  /**
   * Load materials from xml file if read_materials_from_xml is set to true,
   * and a materials xml filename has been passed
   */
  if (_materials_from_xml) {
    if (!isParamSetByUser("materials_xml_file")) {
      paramWarning("materials_xml_file",
                   "read_materials_from_xml is set to true, but "
                   "materials_xml_file is not set! Defaulting to " +
                       _materials_xml_file);
    }
    /// Populate _mat_definitions with materials from openmc xml
    read_material_xml_data();
  }
  /// Check a corresponding material exists for all mesh blocks
  // checkMaterialsExist();
  //
  if (isParamSetByUser("output_inventory_time")) {
    _output_inventory_time = (getParam<double>("output_inventory_time"));
  } else {
    _output_inventory_time = _fp_schedule_uo->getCumulativeTimes().back();
  }
}

void FispactProblem::resolveFispactUserObjects() {

  _fp_nuclear_data_uo =
      &getUserObject<FispactNuclearDataPaths>(_fp_nuclear_data_uo_name);
  _fp_schedule_uo = &getUserObject<FispactSchedule>(_fp_schedule_uo_name);
  _fp_flux_input_uo = &getUserObject<FispactFluxInput>(_fp_flux_uo_name);

  TheWarehouse::Query uo_query =
      theWarehouse().query().condition<AttribSystem>("UserObject");
  std::vector<UserObject *> userobjs;
  uo_query.queryInto(userobjs);

  for (const auto u : userobjs) {
    if (u->type() == "FispactMaterial") {

      _fp_fispact_materials.push_back(dynamic_cast<FispactMaterial *>(u));
    }
  }
}

void FispactProblem::timestepSetup() {

  ExternalProblem::timestepSetup();

  // if (_comm_photon_flux) {
  //   if (timeStep() > 1) {
  //
  //     _segment.destroy<PhotonSharingData>("photon_sharing_instance");
  //   }
  // }
}

void FispactProblem::externalSolve() {

  if (!_solved) {

    const auto &subdomain_names = getParam<std::vector<SubdomainName>>("block");
    auto mesh_subdomains_vec =
        MooseMeshUtils::getSubdomainIDs(_mesh, subdomain_names);
    std::set<SubdomainID> mesh_subdomains(mesh_subdomains_vec.begin(),
                                          mesh_subdomains_vec.end());
    int counter = 1;
    for (const libMesh::Elem *element : *_mesh.getActiveLocalElementRange()) {

      /// Get element id
      dof_id_type global_elem_id = element->id();

      std::vector<double> input_flux =
          _fp_flux_input_uo->getElemFlux(global_elem_id);

      /// Output how many elements have been checked
      _console << std::endl
               << "Elem ID: " << std::to_string(global_elem_id) << std::endl;
      _console << counter++ << "/" << _mesh.getMesh().n_active_local_elem()
               << std::endl;

      // If the energy groups of our nuclear data and input flux do not match,
      // convert input flux to energy grouping of loaded nuclear data
      if (_convert_energy_groups) {
        convertFluxEnergyGroups(input_flux,
                                _fp_flux_input_uo->getFluxEnergyGroups());
      }

      const FispactMaterial &input_material =
          getElementMaterial(global_elem_id);

      setFispactInputData(input_material, input_flux, element->volume(),
                          _fp_ctxt->getInput());

      if (mesh_subdomains.count(element->subdomain_id())) {

        /// Run FISPACT
        _fp_ctxt->process();

        /// Loop over number of inventories to get all data for current
        /// element
        for (int inv_index = 0; inv_index < *_n_solution_inventories;
             inv_index++) {

          /// Vector to store photon energy spectra in photons/cc-s
          std::vector<double> element_photon_energy_spectrum;

          /// Calculated Fispact inventories start at index 1, 0 is reserved
          /// for initial concentrations
          convertGammaEvToCount(
              _fp_ctxt->getOutput().getGammaSpectrumBins(inv_index + 1),
              element_photon_energy_spectrum);

          std::copy(element_photon_energy_spectrum.begin(),
                    element_photon_energy_spectrum.end(),
                    _photon_energy_spectra->spectrum_begin(
                        inv_index, _local_elem_index[global_elem_id]));

          insertElementStrength(inv_index, element,
                                element_photon_energy_spectrum);
        }

        // if inventory manager exists, store requested nuclide metrics
        if (hasUserObject("inv_manager")) {
          FispactInventoryManager &inv_manager =
              getUserObject<FispactInventoryManager>("inv_manager");
          inv_manager.extractInventoryData(*_fp_ctxt, global_elem_id);
        }
      }
    }
    calculateLocalDomainStrength();

    if (_write_photon_flux) {

      const std::vector<double> &inv_times =
          _fp_schedule_uo->getCumulativeTimes();
      writePhotonFlux(_photon_flux_filename, inv_times);
    }

    _fp_ctxt->globalFinalise();

    _solved = true;
  }
}

size_t FispactProblem::getFispactInventoryIndexFromTime() {
  const std::vector<double> &schedule_times =
      _fp_schedule_uo->getCumulativeTimes();

  double inventory_time =
      isTransient() ? time() : getParam<double>("output_inventory_time");

  auto schedule_iterator =
      std::find(schedule_times.begin(), schedule_times.end(), inventory_time);

  if (schedule_iterator == schedule_times.end()) {
    mooseError("Current time " + std::to_string(inventory_time));
  }

  size_t inventory_idx =
      std::distance(schedule_times.begin(), schedule_iterator);

  // Fispact Inventory index "0" is the initial inventory
  return inventory_idx + 1;
}

void FispactProblem::syncSolutions(ExternalProblem::Direction direction) {
  if (direction == ExternalProblem::Direction::FROM_EXTERNAL_APP) {

    if (_comm_photon_flux) {
#ifdef LIBMESH_HAVE_BOOST

      /// Find the the inventory index associated with current time
      /// -1 required due to indexing differences between the stored photon
      /// spectra, and FISPACTs own indexing scheme
      size_t photon_spectra_idx = getFispactInventoryIndexFromTime() - 1;

      /// Calculate TotalDomainStrength
      getTotalDomainStrength();

      /// Create a shared instantiation of the photon sharing class, the init
      /// data that are common across all timesteps
      if (timeStep() == 1) {
        _photon_sharing_instance = _segment.construct<PhotonSharingData>(
            "photon_sharing_instance")(_segment);
        _photon_sharing_instance->setPhotonBins(_photon_bins);
        _photon_sharing_instance->setNumPhotonBins(_n_photon_bins);
        _photon_sharing_instance->setNumLocalElems(_mesh.nActiveLocalElem());
        _photon_sharing_instance->setLocalElemIdMap(_local_elem_index);
      }

      // Setup data that varies per timestep
      _photon_sharing_instance->_is_setup = false;

      _photon_sharing_instance->setPhotonSpectra(
          _photon_energy_spectra->time_begin(photon_spectra_idx),
          _photon_energy_spectra->time_end(photon_spectra_idx));

      _photon_sharing_instance->setElementStrengths(
          _element_strengths.begin() +
              (photon_spectra_idx * _mesh.nActiveLocalElem()),
          _element_strengths.begin() +
              (photon_spectra_idx * _mesh.nActiveLocalElem()) +
              _mesh.nActiveLocalElem());

      _photon_sharing_instance->setLocalDomainStrength(
          _local_domain_strength[photon_spectra_idx]);
      _photon_sharing_instance->setTotalDomainStrength(
          _total_domain_strength[photon_spectra_idx]);

      _photon_sharing_instance->setSamplingMethod(_uniform);

#else
      mooseError("_comm_photon_flux is set to true but libmesh was not
                 built with BOOST. No communication occuring.");
#endif
    }
  }

  if (direction == ExternalProblem::Direction::TO_EXTERNAL_APP) {
  }
}

void FispactProblem::convertFluxEnergyGroups(
    std::vector<double> &flux, const std::vector<double> &input_energy_groups) {
  switch (getParam<MooseEnum>("conversion_type")) {
  case 0: // LETHARGY
    flux = _fp_ctxt->getUtils().GroupConvertByLethargy(
        input_energy_groups, flux, _flux_energy_groups);

  case 1: // ENERGY
    flux = _fp_ctxt->getUtils().GroupConvertByEnergy(input_energy_groups, flux,
                                                     _flux_energy_groups);
  }
}

/**
 * TODO: Break up this function into setFispactInputFlux and
 * setFispactInputMaterial
 */
void FispactProblem::setFispactInputData(const FispactMaterial &material,
                                         const std::vector<double> &flux,
                                         const double &volume,
                                         IFispactInputDataBase &input) const {
  input.setFlux(_flux_energy_groups, flux);
  input.setFluxWallLoading(1.0);
  input.setFluxName("neutrons");

  input.setSolverTolerance(_rtol, _atol);

  input.setExcludeXrays(_exclude_xrays);

  /// Get density from mat density, in g/cm^3!
  double density = material.getDensity();
  input.setDensity(density);

  /// Set atoms threshold
  input.setAtomsThreshold(1.0e1);

  double total_mass_grams = density * volume;

  if (material.getMaterialType() == "MASS") {

    std::vector<int> atomic_numbers;
    std::vector<double> percent;

    /// Set total mass, in kg!
    input.setMassTotal(total_mass_grams * 1e-3);

    const std::unordered_map<std::string, double> &nuclideFractionMap =
        material.getNuclideFractionMap();

    atomic_numbers.reserve(nuclideFractionMap.size());

    for (auto &[element_name, mass_fraction] : nuclideFractionMap) {

      atomic_numbers.push_back(
          _fp_ctxt->getUtils().GetAtomicNumberFromElementName(element_name));
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

    /// For all key (isotope name) value (mass_fraction) pairs in map,
    /// calculate the number of atoms pertaining to each isotope and append to
    /// input fuel
    for (const auto &[isotope_name, mass_fraction] : nuclideFractionMap) {

      double zai_mass = total_mass_grams * (mass_fraction);

      zais.push_back(_fp_ctxt->getUtils().GetZai(isotope_name));

      atoms.push_back(
          getNumAtoms(zai_mass, _molar_mass_map.at(isotope_name), AVOGADRO));
    }
    input.setFuel(zais, atoms);
  }

  double flux_sum = std::accumulate(flux.begin(), flux.end(), 0.0);

  setFispactSchedule(input, volume, flux_sum);
}

void FispactProblem::setFispactSchedule(IFispactInputDataBase &input,
                                        const double &element_volume,
                                        const double &neutron_flux_sum) const {

  std::vector<double> flux_schedule = _fp_schedule_uo->getFluxAmplitude();
  const std::vector<double> &times = _fp_schedule_uo->getTimes();

  /// Need to scale input flux amplitude by the total neutron flux in this
  /// element/ volume
  for (auto &flux_amplitude : flux_schedule) {
    flux_amplitude *= neutron_flux_sum / element_volume;
  }
  input.setSchedule(times, flux_schedule);
}

const FispactMaterial &
FispactProblem::getElementMaterial(const dof_id_type &elem_id) {

  libMesh::Elem *elem = _mesh.elemPtr(elem_id);

  // Query the warehouse to see if a FispactMaterial exists on the block this
  // element is assigned to
  std::vector<UserObject *> objs;
  theWarehouse()
      .query()
      .condition<AttribSystem>("UserObject")
      .condition<AttribSubdomains>(elem->subdomain_id())
      .queryInto(objs);

  // Remove extraneous user objects that are not FispactMaterials. Having done
  // this, only one object should remain in the vector, and it should be the
  for (auto it = objs.begin(); it != objs.end();) {
    if ((*it)->type() != "FispactMaterial") {
      it = objs.erase(it);
    } else {
      ++it;
    }
  }

  // FISPACT material pertaining to this block. If there ismore than one
  // object, then two FispactMaterials are assigned to this block, and that
  // makes no blimmin sense does it
  if (objs.empty()) {
    mooseError("Unable to find FispactMaterial object on block " +
               std::to_string(elem->subdomain_id()));
  } else if (objs.size() > 1) {
    mooseError("More than 1 FispactMaterial definition exists on block " +
               std::to_string(elem->subdomain_id()));
  }

  /// Return the FispactMaterial
  return *(static_cast<FispactMaterial *>(objs[0]));
}

bool FispactProblem::isFlux(const std::vector<double> &flux) const {
  /**
   * return true if there is flux, false if there isn't, as the function name
   * implies
   */
  bool is_zero_flux =
      std::all_of(flux.begin(), flux.end(), [](double j) { return j == 0; });

  return !is_zero_flux;
}

void FispactProblem::convertGammaEvToCount(
    const std::vector<double> &photon_energy_spectra_ev,
    std::vector<double> &photon_energy_spectra_per_s) {

  /// Reserve memory for photons per cc per s vector
  photon_energy_spectra_per_s.resize(photon_energy_spectra_ev.size(), 0);

  /// Assert in case, somehow, the spectra output by FISPACT has different bin
  /// structure to _photon_bins??
  mooseAssert(
      _photon_bins.size() == photon_energy_spectra_ev.size() + 1,
      "_photon_bins size should be one less that input photon spectra!");

  /// Convert from MeV/s to per cc per s for each bin
  for (int i = 0; i < photon_energy_spectra_ev.size(); i++) {
    double bin_energy =
        _photon_bins[i] + ((_photon_bins[i + 1] - _photon_bins[i]) / 2);

    /**
     * per_cc_per_s = MeV/s * (inventory_density/(inventory_mass *
     * energy_bin_midpoint))
     */
    /// _photon_bins are in eV, but Photon Spectra is returned in by FISPACT in
    /// MeV, so need 1e6
    photon_energy_spectra_per_s[i] =
        photon_energy_spectra_ev[i] * (1e6 / bin_energy);
  }
}

void FispactProblem::setPhotonBins(const std::vector<double> &photon_bins) {
  _photon_bins = photon_bins;
}

double FispactProblem::calculateElementStrength(
    const std::vector<double> element_flux) {
  double element_strength = 0;

  for (double flux_bin : element_flux) {
    element_strength += flux_bin;
  }

  return element_strength;
}

void FispactProblem::insertElementStrength(
    const int inv_index, const libMesh::Elem *element,
    const std::vector<double> element_flux) {
  double elem_strength = calculateElementStrength(element_flux);

  int idx =
      (inv_index * _mesh.nActiveLocalElem()) + _local_elem_index[element->id()];

  _element_strengths[idx] = elem_strength;
}

void FispactProblem::calculateLocalDomainStrength() {

  for (int inv_index = 0; inv_index < *_n_solution_inventories; inv_index++) {
    for (const libMesh::Elem *element : *_mesh.getActiveLocalElementRange()) {

      double element_strength;

      // Attempt to retrieve element strength from vector
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
      ((_n_photon_bins * sizeof(double)) + sizeof(int)) * n_local_elem;

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
  // std::string log_name = "FISPACT_app_" +
  //                        this->getMooseApp().getInputFileNames()[0] +
  //                        std::to_string(processor_id()) + ".log";

  std::string log_name =
      "FISPACT_app_" + std::to_string(processor_id()) + ".log";
  return log_name;
}

void FispactProblem::setNuclearData() {
  _fp_ctxt->setNuclearData(_fp_nuclear_data_uo->getNuclearDataPathMap());
}

const std::string FispactProblem::generateInterprocessName() {
  char mpi_proc_name[MPI_MAX_PROCESSOR_NAME];
  int len = 0;
  int err = MPI_Get_processor_name(mpi_proc_name, &len);

  MPI_Comm node_comm;
  int local_rank;
  MPI_Comm_split_type(comm().get(), MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                      &node_comm);
  MPI_Comm_rank(node_comm, &local_rank);

  std::string ipc_name = std::string(mpi_proc_name);
  ipc_name += "_" + std::to_string(local_rank);

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

    std::pair<std::string, double> key_value =
        std::pair<std::string, double>(element_symbols[i], molar_masses[i]);
    _molar_mass_map.insert(key_value);
  }
}

void FispactProblem::checkForEnergyGroupConsistency() {
  size_t n_nd_energy_groups = _fp_ctxt->getNuclearDataCrossSections();

  if (_fp_flux_input_uo->getNumEnergyGroups() != n_nd_energy_groups) {
    _flux_energy_groups =
        _fp_ctxt->getUtils().getNeutronEnergyBounds(n_nd_energy_groups);
    _convert_energy_groups = true;

    std::string conversion_type = getParam<MooseEnum>("conversion_type");
    mooseWarning("Input neutron flux group structure does not match that of "
                 "input nuclear data, converting using " +
                 conversion_type);
  } else {
    _flux_energy_groups = _fp_flux_input_uo->getFluxEnergyGroups();
  }
}

// TO DO: Move most of this functionality to FispactMaterial, allow
// FispactMaterial to take in a wider range of inputs eg density units, and do
// the hard work there
void FispactProblem::read_material_xml_data() {
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(_materials_xml_file.c_str());

  if (!result) {
    mooseError("No file called " + _materials_xml_file +
               " could be found, exiting.");
  }

  for (pugi::xml_node material : doc.child("materials").children()) {

    InputParameters params =
        _app.getFactory().getValidParams("FispactMaterial");

    /// Get material name
    std::string material_name = material.attribute("name").value();

    /// Get material density
    double density =
        std::stod(material.child("density").attribute("value").value());

    std::string density_units =
        material.child("density").attribute("units").value();

    std::vector<std::string> nuclides;
    std::vector<double> nuclide_fractions;

    for (pugi::xml_node nuclide : material.children("nuclide")) {
      nuclides.push_back(std::string(nuclide.attribute("name").value()));

      if (nuclide.attribute("wo")) {
        nuclide_fractions.push_back(std::stod(nuclide.attribute("wo").value()) /
                                    100);
      } else if (nuclide.attribute("ao")) {

        nuclide_fractions.push_back(
            -std::stod(nuclide.attribute("ao").value()));
      }
    }

    bool all_wo =
        std::all_of(nuclide_fractions.begin(), nuclide_fractions.end(),
                    [](double x) { return x >= 0.0; });
    bool all_ao =
        std::all_of(nuclide_fractions.begin(), nuclide_fractions.end(),
                    [](double x) { return x <= 0.0; });

    if (!(all_wo || all_ao)) {
      mooseError("Cannot mix atom and weight percents in material. Error when "
                 "parsing material xml");
    }

    std::string fraction_type;
    if (all_ao) {
      fraction_type = "ao";
      for (auto &fraction : nuclide_fractions) {
        fraction = -fraction;
      }
    }
    if (all_wo) {
      fraction_type = "wo";
    }

    params.set<double>("density") = density;
    params.set<MooseEnum>("material_type") = "FUEL";
    params.set<MooseEnum>("fraction_type") = fraction_type;
    params.set<MooseEnum>("density_units") = density_units;
    params.set<std::vector<SubdomainName>>("block") = {material_name};
    params.set<std::vector<std::string>>("nuclides") = nuclides;
    params.set<std::vector<double>>("nuclide_fraction") = nuclide_fractions;

    addUserObject("FispactMaterial", material_name, params);
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
  bool parallel = false;
  mooseWarning(
      "Writing photon flux requires HDF5 compiled with MPI. "
      "Simulation continuing, but photon spectra will not be written.");
  return;
#endif
  /// Open OpenMC statepoint file with Neutron Flux
  hid_t file_id =
      hdf5_utils::file_open(filename.c_str(), 'w', parallel, comm().get());

  for (int inv_id = 0; inv_id < *_n_solution_inventories; inv_id++) {

    std::string dataset_name =
        "photon_flux_" + std::to_string(inventory_times[inv_id]);

    int ndim = 2;
    unsigned long x_dim = _mesh.getMesh().n_active_elem();
    int y_dim = _n_photon_bins;
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
      hsize_t count[2] = {1, _n_photon_bins};

      hdf5_utils::write_double_hyperslab(
          h5_dataset, nullptr, ndim, hyperslab_dims, start, count,
          &*_photon_energy_spectra->spectrum_begin(
              inv_id, _local_elem_index[element->id()]),
          parallel);
    }
    /// Close all the HDF5 bits and pieces
    H5Dclose(h5_dataset);
  }

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

const double FispactProblem::avogadroNumber() const { return AVOGADRO; }

void FispactProblem::callFispactFactory() { _fp_ctxt = createFispactContext(); }
