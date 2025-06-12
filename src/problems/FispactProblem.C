#include "FispactProblem.h"

#include "FispactSchedule.h"
/// HDF5 include
#include "H5Cpp.h"

/// PugiXML include
#include "MooseError.h"
#include "MooseTypes.h"
#include "fispactinputdata.hpp"
#include "pugixml.hpp"
#include <filesystem>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>

registerMooseObject("FizzyApp", FispactProblem);

InputParameters FispactProblem::validParams() {
  InputParameters params = ExternalProblem::validParams();

  // Suppress parameters we are not using
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

  // New parameters we need for FISPACT
  params.addRequiredParam<std::string>("neutron_flux_file",
                                       "HDF5 file storing neutron flux");
  params.addRequiredParam<std::string>(
      "neutron_flux_hdf5_path",
      "Path within HDF5 file for the vector storing neutron flux");

  params.addRequiredParam<std::string>("fispact_nuclear_data_path",
                                       "Path to FISPACT nuclear data");

  params.addRequiredParam<std::string>("fispact_nuclear_data_path",
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
  params.addParam<std::string>(
      "materials_xml_file", "materials.xml",
      "Path and name of material xml file user wishes to use.");

  return params;
}

FispactProblem::FispactProblem(const InputParameters &params)
    : ExternalProblem(params), _fp_monitor(fispactLogName()),
      _fp_nuclear_data(_fp_monitor),
      _neutron_flux_filename(getParam<std::string>("neutron_flux_file")),
      _neutron_flux_hdf5_path(getParam<std::string>("neutron_flux_hdf5_path")),
      _materials_from_xml(getParam<bool>("read_materials_from_xml")),
      _neutron_bin_type(getParam<std::string>("neutron_bin_type")),
      _schedule_uo_name(getParam<std::string>("fispact_schedule_uo")) {

  // Initialise FISPACT
  fp::GlobalInitialise(_fp_monitor);

  // Load materials from xml file
  if (_materials_from_xml) {
    // Set file to get materials from
    _materials_xml_file = getParam<std::string>("materials_xml_file");
    // Populate _mat_definitions with materials from openmc xml
    read_material_xml_data();
  } else {
    if (isParamSetByUser("materials_xml_file")) {
      mooseError("materials_xml_file is set by user, but "
                 "read_materials_from_xml is false!");
    }
  }

  checkMaterialsExist();

  // Get the path to our nuclear data
  std::string fp_nuclear_data_path =
      getParam<std::string>("fispact_nuclear_data_path");

  // Set nuclear data paths
  setNuclearData(fp_nuclear_data_path);

  // Set neutron bin type
  setNeutronBins();

  // Read neutron flux from HDF5
  readNeutronFluxFromHDF5(_neutron_flux_filename, _neutron_flux_hdf5_path);
}

void FispactProblem::externalSolve() {
  // Set up fispact input data
  fp::InputData fispact_input(_fp_monitor);
  fp::OutputData fispact_output(_fp_monitor);

  int counter = 0;
  for (MeshBase::element_iterator element_iter =
           _mesh.activeLocalElementsBegin();
       element_iter != _mesh.activeLocalElementsEnd(); element_iter++) {

    std::cout << counter << "/" << _mesh.getMesh().n_active_local_elem()
              << std::endl;
    // Get element id
    int elem_id = (*element_iter)->id();

    // Check if there is any neutron flux in current element
    bool is_zero_flux = std::all_of(_neutron_fluxes[elem_id].begin(),
                                    _neutron_fluxes[elem_id].end(),
                                    [](double j) { return j == 0; });

    // If there is flux in the element, run FISPACT
    if (!is_zero_flux) {
      // Here we are assuming the input mesh is in centremeters
      double element_volume = (*element_iter)->volume();

      MaterialDefinition el_mat = getElementMaterial(elem_id);

      setFispactInputData(_fp_monitor, fispact_input, el_mat,
                          _neutron_fluxes[elem_id], _neutron_bins,
                          element_volume);
      // Run FISPACT!
      fp::Process(fispact_input, _fp_nuclear_data, fispact_output, _fp_monitor,
                  process_callback);
    }

    counter++;
  }
}

std::string FispactProblem::fispactLogName() {
  std::string log_name =
      "FISPACT_app_" + std::to_string(processor_id()) + ".log";
  return log_name;
}

void FispactProblem::setNuclearData(std::string nd_base_path) {
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

void FispactProblem::readNeutronFluxFromHDF5(std::string filename,
                                             std::string tally_dir) {
  // Open statepoint file as H5File
  H5::H5File h5f(filename.c_str(), H5F_ACC_RDONLY);

  // Create dataset and dataspace for the tally results
  H5::DataSet dataset_tally = h5f.openDataSet(tally_dir.c_str());
  H5::DataSpace dspace_tally = dataset_tally.getSpace();

  // Read in dimensions of tally results array
  hsize_t tally_array_dims[3];
  dspace_tally.getSimpleExtentDims(tally_array_dims, NULL);

  // Read in number of realizations to calculate mean
  // Set up dataset and dataspace for reading realization count
  int n_realizations;
  H5::DataSet dataset_realizations = h5f.openDataSet("n_realizations");
  H5::DataSpace dspace_realizations = dataset_realizations.getSpace();

  // Set up memory space for reading realization (batch) count
  hsize_t memspace_dimensions_realizations[3] = {1, 1, 1};
  H5::DataSpace memspace_realizations(1, memspace_dimensions_realizations);
  dataset_realizations.read(&n_realizations, H5::PredType::STD_I32LE,
                            memspace_realizations, dspace_realizations);

  // Set up vector of vectors to store neutron fluxes
  std::vector<std::vector<double>> neutron_fluxes(
      tally_array_dims[0] / _num_neutron_bins,
      std::vector<double>(_num_neutron_bins, 0));
  for (int i = 0; i < neutron_fluxes.size(); i++) {
    // Set up std::vector to store neutron flux data
    std::vector<double> neutron_flux_data(_num_neutron_bins, 0.0);

    // Set up counts and offsets for selecting hyperslab of tally array
    hsize_t dataCount[3] = {static_cast<hsize_t>(_num_neutron_bins), 1, 1};
    hsize_t dataOffset[3] = {static_cast<hsize_t>((_num_neutron_bins * i)), 0,
                             0};

    // Set up memory space for reading tally results
    hsize_t arr_len[3] = {static_cast<hsize_t>(_num_neutron_bins), 1, 1};
    H5::DataSpace memspace_tally(1, arr_len);
    dspace_tally.selectHyperslab(H5S_SELECT_SET, dataCount, dataOffset);

    // Read in neutron flux tally results
    dataset_tally.read(neutron_flux_data.data(), H5::PredType::IEEE_F64LE,
                       memspace_tally, dspace_tally);

    for (auto &neutron_flux : neutron_flux_data) {
      neutron_flux = neutron_flux / n_realizations;
    }
    neutron_fluxes[i] = neutron_flux_data;
  }

  // for now, only keep neutron fluxes for local elements
  for (MeshBase::element_iterator it = _mesh.activeLocalElementsBegin();
       it != _mesh.activeLocalElementsEnd(); it++) {
    int elem_id = (*it)->id();
    _neutron_fluxes.insert(std::make_pair(elem_id, neutron_fluxes[elem_id]));
  }
}

// To do: Break up this function into setFispactInputFlux and
// setFispactInputMaterial
void FispactProblem::setFispactInputData(fp::FispactMonitor &monitor,
                                         fp::InputData &input,
                                         MaterialDefinition &material,
                                         std::vector<double> &neutron_flux,
                                         const std::vector<double> &bins,
                                         double volume) {
  // Set neutron flux
  input.setFlux(bins, neutron_flux);

  input.setFluxWallLoading(1.0);

  input.setFluxName("neutrons");

  // get density from mat density
  input.setDensity(material._mat_density);
  input.setAtomsThreshold(1.0e3);

  // Volume read in is in cm ^ 3, so we need to scale by 1e-6, as mass is in kg
  double total_mass = material._mat_density * volume * 1e-6;
  input.setMassTotal(total_mass);

  std::vector<int> atomic_numbers;
  std::vector<double> percent;
  for (auto &element_name_ao_pair : material._mat_atomic_composition) {
    std::string element_name = element_name_ao_pair.first;

    /** Remove numbers from element name, isotope doesn't matter here as we're
      obtaining the atomic number, not the atomic mass **/
    element_name.erase(
        std::remove_if(element_name.begin(), element_name.end(),
                       [](unsigned char c) { return std::isdigit(c); }),
        element_name.end());

    atomic_numbers.push_back(
        fp::util::GetAtomicNumberFromElementName(monitor, element_name));
    percent.push_back(element_name_ao_pair.second);
  }

  input.setMass(atomic_numbers, percent);

  setFispactSchedule(input);
}

void FispactProblem::setFispactSchedule(fp::InputData &input) {
  FispactSchedule &schedule = getUserObject<FispactSchedule>(_schedule_uo_name);

  const std::vector<double> &flux_schedule = schedule.getFluxSchedule();
  const std::vector<double> &times = schedule.getTimes();
  input.setSchedule(times, flux_schedule);
}

FispactProblem::MaterialDefinition &
FispactProblem::getElementMaterial(int &elem_id) {
  libMesh::Elem *elem = _mesh.elemPtr(elem_id);

  const std::string &subdomain_name =
      _mesh.getSubdomainName(elem->subdomain_id());

  // Return material definition if it exists, otherwise throw error
  if (_mat_definitions.find(subdomain_name) != _mat_definitions.end()) {
    return _mat_definitions.at(subdomain_name);

  } else {
    mooseError("No FISPACT material named " + subdomain_name + " was found.");
  }
}

void FispactProblem::read_material_xml_data() {
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(_materials_xml_file.c_str());

  if (!result) {
    // mooseError();
    mooseError("No file called " + _materials_xml_file +
               " could be found, exiting.");
  }
  for (pugi::xml_node material : doc.child("materials").children()) {

    // Get material name
    std::string material_name = material.attribute("name").value();
    // Get material density
    double density =
        std::stod(material.child("density").attribute("value").value());

    // Get atomic composition of material
    std::vector<std::pair<std::string, double>> atomic_comp;

    for (pugi::xml_node nuclide : material.children("nuclide")) {
      std::pair<std::string, double> nuclide_symbol_and_percentage =
          std::make_pair(std::string(nuclide.attribute("name").value()),
                         std::stod(nuclide.attribute("ao").value()));
      atomic_comp.push_back(nuclide_symbol_and_percentage);
    }

    // Create material definition
    MaterialDefinition material_def{material_name, atomic_comp, density};
    // Insert material definition into material map
    _mat_definitions.insert(std::make_pair(material_name, material_def));
  }
}

void FispactProblem::setNeutronBins() {
  if (_neutron_bin_type == "G1102") {
    _neutron_bins = fp::groups::G1102();
    _num_neutron_bins = 1102;
  }
}

void FispactProblem::checkMaterialsExist() {
  std::vector<unsigned short> subdomain_ids;
  for (auto &subdomain_id : _mesh.meshSubdomains()) {
    subdomain_ids.push_back(subdomain_id);
  }

  std::vector<SubdomainName> subdomain_names =
      _mesh.getSubdomainNames(subdomain_ids);

  for (auto &subdomain_name : subdomain_names) {
    if (_mat_definitions.find(subdomain_name) == _mat_definitions.end()) {
      mooseError("Block " + subdomain_name +
                 " does not have a corresponding FISPACT material defined");
    }
  }
}
