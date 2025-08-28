#include "ExternalProblem.h"
#include "FispactProblem.h"

// Custom user object includes
#include "FispactSchedule.h"

/// HDF5 include
#include "H5Cpp.h"
#include <hdf5/openmpi/H5Dpublic.h>
#include <hdf5/openmpi/H5FDmpio.h>
#include <hdf5/openmpi/H5Fpublic.h>
#include <hdf5/openmpi/H5Ipublic.h>
#include <hdf5/openmpi/H5Ppublic.h>
#include <hdf5/openmpi/H5Spublic.h>
#include <hdf5/openmpi/H5Tpublic.h>
#include <hdf5/openmpi/H5public.h>
#include <hdf5/openmpi/H5version.h>

#include "MooseError.h"
#include "MooseTypes.h"

// Fispact includes
#include "fispactcompute.hpp"
#include "fispactconstantsapi.h"
#include "fispactgroupstructures.hpp"
#include "fispactinputdata.hpp"
#include "fispactnucleardata.hpp"
#include "fispactoutputdata.hpp"
#include "fispactoutputdataapi.h"
#include "fispactutil.hpp"

#include "libmesh/elem.h"
#include "mpi.h"

/// PugiXML include
#include "pugixml.hpp"

// Cpp includes
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
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

  params.addParam<bool>(
      "write_photon_flux", false,
      "Boolean value used to determine whether to wite photon spectra to hdf5 "
      "after Fizzy has finished running.");

  params.addParam<std::string>(
      "photon_flux_filename", "photon_flux.h5",
      "Filename for the h5 file containing the output photon spectra");

  params.addParam<bool>(
      "comm_photon_flux", false,
      "Boolean value used to indicate whether to use boost::interprocess to "
      "communicate photon spectra through IPC.");
  return params;
}

FispactProblem::FispactProblem(const InputParameters &params)
    : ExternalProblem(params), _fp_monitor(fispactLogName()),
      _fp_nuclear_data(_fp_monitor),
      _neutron_flux_filename(getParam<std::string>("neutron_flux_file")),
      _neutron_flux_hdf5_path(getParam<std::string>("neutron_flux_hdf5_path")),
      _photon_flux_filename(getParam<std::string>("photon_flux_filename")),
      _materials_from_xml(getParam<bool>("read_materials_from_xml")),
      _neutron_bin_type(getParam<std::string>("neutron_bin_type")),
      _schedule_uo_name(getParam<std::string>("fispact_schedule_uo")),
      _local_domain_strength(0), _total_domain_strength(0),
      _write_photon_flux(getParam<bool>("write_photon_flux")),
      _comm_photon_flux(getParam<bool>("comm_photon_flux")) {

  // Load materials from xml file if read_materials_from_xml is set to true, and
  // a materials xml filename has been passed
  if (_materials_from_xml) {
    if (!isParamSetByUser("materials_xml_file")) {
      mooseError("read_materials_from_xml is set to true, but "
                 "materials_xml_file is not set!");
    }
    // Set file to get materials from
    _materials_xml_file = getParam<std::string>("materials_xml_file");
    // Populate _mat_definitions with materials from openmc xml
    read_material_xml_data();
  }

  // If write_photon_flux was set to true then check that user input a filename,
  // if not use default
  if (_write_photon_flux && !isParamSetByUser("photon_flux_filename")) {
    _console << "write_photon_flux is set to true but photon_flux_filename is "
                "not set! Photon flux filename defaulting to " +
                    _photon_flux_filename;
  }

  // Check a corresponding material exists for all mesh blocks
  checkMaterialsExist();

  // Initialise FISPACT
  fp::GlobalInitialise(_fp_monitor);

  // Get the path to our nuclear data
  std::string fp_nuclear_data_path =
      getParam<std::string>("fispact_nuclear_data_path");

  // Set nuclear data paths
  setNuclearData(fp_nuclear_data_path);

  // Set neutron bin type
  setNeutronBins();
  // Read neutron flux from h5 file
  readNeutronFluxFromHDF5(_neutron_flux_filename, _neutron_flux_hdf5_path);
  if (_comm_photon_flux) {
#ifdef LIBMESH_HAVE_BOOST
    // Give the shared memory region a name based on current MPI rank to
    // prevent clashes
    std::string data_name = "SHARING_DATA_" + std::to_string(comm().rank());

    // Remove any shared memory region with a similar name just in case
    bi::shared_memory_object::remove(data_name.c_str());

    // Calculate space needed for shared mem region
    unsigned long shared_memory_size = calculateMemorySize();

    // Create shared memory region
    _segment = bi::managed_shared_memory(bi::create_only, data_name.c_str(),
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

      unsigned long shared_memory_size = calculateMemorySize();

      std::string data_name = "SHARING_DATA_" + std::to_string(comm().rank());

      // vector of all local domain strengths
      std::vector<double> local_domain_strengths;
      comm().allgather(_local_domain_strength, local_domain_strengths);

      // calculate TotalDomainStrength
      getTotalDomainStrength();

      // Create a shared instantiation of the photon sharing class
      PhotonSharingData *photon_sharing_instance =
          _segment.construct<PhotonSharingData>(
              "PhotonSharingData photon_sharing_instance")(
              _segment, _photon_fluxes, _element_strengths, 24,
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

  // Set up fispact input data
  fp::InputData fispact_input(_fp_monitor);
  fp::OutputData fispact_output(_fp_monitor);

  int counter = 1;
  for (auto element_iter : *_mesh.getActiveLocalElementRange()) {

    _console << "starting " + std::to_string(counter) << std::endl;

    // Get element id
    int elem_id = element_iter->id();
    _console << "Elem ID" << std::to_string(elem_id) << std::endl;
    std::vector<double> photon_spectra(24, 0);

    // Check if there is any neutron flux in current element
    bool is_flux = isFlux(elem_id);

    // If there is flux in the element, run FISPACT
    if (is_flux) {
      // Here we are assuming the input mesh is in centimeters
      double element_volume = (element_iter)->volume();

      // Get element material definition
      MaterialDefinition &el_mat = getElementMaterial(elem_id);

      setFispactInputData(_fp_monitor, fispact_input, el_mat,
                          _neutron_fluxes[elem_id], _neutron_bins,
                          element_volume);
      // Run FISPACT
      fp::Process(fispact_input, _fp_nuclear_data, fispact_output, _fp_monitor,
                  process_callback);

      convertGammaEvToCount(
          fispact_input, fispact_output.getGammaSpectrumBins(1),
          fispact_output.getGammaSpectrumBoundaries(0), photon_spectra);

      for (auto &ev : fispact_input.getGammaEnergyBounds()) {
        _console << std::to_string(ev) << " ";
      }
    }

    _photon_fluxes.insert(
        std::pair<int, std::vector<double>>(elem_id, photon_spectra));

    // Calculate photon source strength of this element and insert it into
    // _element_strengths map
    insertElementStrength(element_iter, _photon_fluxes.at(elem_id));

    // Add calculated element strength to _local_domain_strength
    updateLocalDomainStrength(element_iter);

    // Output how many elements have been checked
    _console << counter++ << "/" << _mesh.getMesh().n_active_local_elem()
             << std::endl;
  }

  _photon_bins = fispact_output.getGammaSpectrumBoundaries(0);

  for (auto &bin : _photon_bins) {
    if (bin < 0.01) {
      continue;
    }
    // Change bins to MeV;
    bin *= 1e6;
  }

  writePhotonFluxToHDF5(_photon_flux_filename, fispact_output);
  fp::GlobalFinalise(_fp_monitor);
}

std::string FispactProblem::fispactLogName() {
  std::string log_name = "FISPACT_app_" +
                         this->getMooseApp().getInputFileNames()[0] +
                         std::to_string(processor_id()) + ".log";
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

void FispactProblem::writePhotonFluxToHDF5(const std::string &filename,
                                           fp::OutputData &fispact_output) {
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Info info = MPI_INFO_NULL;
  _console << "Writing..." << std::endl;
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_all_coll_metadata_ops(plist_id, true);
  H5Pset_coll_metadata_write(plist_id, true);
  // Set HDF5 driver
  H5Pset_fapl_mpio(plist_id, comm, info);

  // Create HDF5 file
  auto file_id =
      H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

  // Loop over all MPI ranks
  std::string dataset_name = "photon_data";
  hsize_t dataspace_dims[2];

  // Set up dimensions for photon flux dataspace
  dataspace_dims[0] = _mesh.getMesh().n_active_elem();
  dataspace_dims[1] = 24;

  hid_t filespace = H5Screate_simple(2, dataspace_dims, NULL);
  hid_t h5_dataset =
      H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_DOUBLE, filespace,
                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Sclose(filespace);

  // Get dataspace corresponding to the hdf5 dataset
  filespace = H5Dget_space(h5_dataset);

  for (auto &element_flux_pair : _photon_fluxes) {

    // CONTEXT: A hyperslab is a section of a hdf5 dataspace, don't be put off
    // by the fancy name

    // Create hdf5 dataspace to be our memory space for writing
    hsize_t memory_dspace_dims[2] = {1, 24};
    hid_t memspace = H5Screate_simple(2, memory_dspace_dims, NULL);

    // Retrieve element id from _photon_fluxes map
    hsize_t element_id = (hsize_t)element_flux_pair.first;

    // offset and count are used to select our hyperslab.
    // Given the dimensions of the dataspace are num_elems * 24, our offset
    // selection should be the element_id who's flux we wish to write
    hsize_t offset[2] = {element_id, 0};

    // count specifies the number of entries we wish to write in each
    // dimension
    hsize_t count[2] = {1, 24};

    // Select the hyperslab of the dataspace that we wish to write to
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    // Set data transfer mode
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

    // Write the data
    H5Dwrite(h5_dataset, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id,
             element_flux_pair.second.data());
    H5Sclose(memspace);
  }
  writePhotonFluxBins(file_id, fispact_output);
  // Close all the HDF5 bits and pieces
  H5Sclose(filespace);
  H5Dclose(h5_dataset);
  H5Pclose(plist_id);
  H5Fclose(file_id);
}

void FispactProblem::writePhotonFluxBins(hid_t file_id,
                                         fp::OutputData &fispact_output) {

  // Get gamma bins from the FISPACT output data
  std::vector<double> photon_bins =
      fispact_output.getGammaSpectrumBoundaries(0);

  // Set up hsize_t object to hold dataset dimensions
  hsize_t bin_dataset_dims[1];
  bin_dataset_dims[0] = photon_bins.size();

  std::string dataset_name = "photon_bins";

  // Create a dataspace for the photon bins
  hid_t filespace_photon_bins = H5Screate_simple(1, bin_dataset_dims, NULL);

  // Create a dataset for the photon bins using the dataspace
  hid_t h5_dataset_photon_bins =
      H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_DOUBLE,
                filespace_photon_bins, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Set up transfer properties
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

  // Only write the photon bins out on rank 0, we don't need to write it from
  // all of them
  if (ExternalProblem::comm().rank() == 0) {
    H5Dwrite(h5_dataset_photon_bins, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
             plist_id, photon_bins.data());
  }
  H5Sclose(filespace_photon_bins);
  H5Dclose(h5_dataset_photon_bins);
  H5Pclose(plist_id);
}

void FispactProblem::readNeutronFluxFromHDF5(std::string filename,
                                             std::string tally_dir) {

  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Info info = MPI_INFO_NULL;
  // Set access properties
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);

  H5Pset_all_coll_metadata_ops(plist_id, true);
  H5Pset_coll_metadata_write(plist_id, true);
  // Set HDF5 driver
  H5Pset_fapl_mpio(plist_id, comm, info);
  // Open statepoint file as H5File
  hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, plist_id);
  H5Pclose(plist_id);
  // Create dataset and dataspace for the tally results
  hid_t dataset_tally_id = H5Dopen(file_id, tally_dir.c_str(), H5P_DEFAULT);
  hid_t space_tally_id = H5Dget_space(dataset_tally_id);

  // Read in dimensions of tally results array
  hsize_t tally_array_dims[3];
  H5Sget_simple_extent_dims(space_tally_id, tally_array_dims, NULL);

  // Check that the neutron flux we are reading is suitable for this mesh
  if (tally_array_dims[0] != _mesh.nElem() * _num_neutron_bins) {
    mooseError("Neutron flux file is incorrect dimension");
  }

  // Read in number of realizations to calculate mean
  // Set up dataset and dataspace for reading realization count
  int n_realizations;
  hid_t dataset_realizations_id =
      H5Dopen(file_id, "n_realizations", H5P_DEFAULT);
  hid_t space_realizations_id = H5Dget_space(dataset_realizations_id);

  // Set up memory space for reading realization (batch) count
  hsize_t memspace_dimensions_realizations[3] = {1, 1, 1};
  hid_t memspace_realizations_id =
      H5Screate_simple(1, memspace_dimensions_realizations, nullptr);

  H5Dread(dataset_realizations_id, H5T_STD_I32LE, memspace_realizations_id,
          space_realizations_id, H5P_DEFAULT, &n_realizations);

  // Close realizations hdf5 structures
  H5Sclose(space_realizations_id);
  H5Sclose(memspace_realizations_id);
  H5Dclose(dataset_realizations_id);

  // Set up vector of vectors to store neutron fluxes
  std::vector<std::vector<double>> neutron_fluxes(
      tally_array_dims[0] / _num_neutron_bins,
      std::vector<double>(_num_neutron_bins, 0));

  for (int i = 0; i < neutron_fluxes.size(); i++) {
    // Set up std::vector to store neutron flux data
    std::vector<double> neutron_flux_data(_num_neutron_bins, 0.0);

    // Set up counts and offsets for selecting hyperslab of tally array
    hsize_t data_offset[3] = {static_cast<hsize_t>((_num_neutron_bins * i)), 0,
                              0};
    hsize_t data_count[3] = {static_cast<hsize_t>(_num_neutron_bins), 1, 1};

    // Set up memory space for reading tally results
    hsize_t arr_len[3] = {static_cast<hsize_t>(_num_neutron_bins), 1, 1};
    hid_t memspace_tally_id = H5Screate_simple(1, arr_len, NULL);
    H5Sselect_hyperslab(space_tally_id, H5S_SELECT_SET, data_offset, NULL,
                        data_count, NULL);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    // Read in neutron flux tally results
    H5Dread(dataset_tally_id, H5T_IEEE_F64LE, memspace_tally_id, space_tally_id,
            H5P_DEFAULT, neutron_flux_data.data());

    for (auto &neutron_flux : neutron_flux_data) {
      neutron_flux = neutron_flux / n_realizations;
    }
    // Add neutron flux data to map
    neutron_fluxes[i] = neutron_flux_data;

    // Close memspace dataspace
    H5Sclose(memspace_tally_id);
    H5Pclose(plist_id);
  }

  // for now, only keep neutron fluxes for local elements
  for (MeshBase::element_iterator it = _mesh.activeLocalElementsBegin();
       it != _mesh.activeLocalElementsEnd(); it++) {
    int elem_id = (*it)->id();
    _neutron_fluxes.insert(std::make_pair(elem_id, neutron_fluxes[elem_id]));
  }
  H5Dclose(dataset_tally_id);
  H5Sclose(space_tally_id);
  H5Fclose(file_id);
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

  // Volume read in is in cm ^ 3, so we need to scale by 1e-6, as mass is in
  // kg
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

bool FispactProblem::isFlux(int elem_id) {
  // Check if there is any neutron flux in current element
  bool is_zero_flux = std::all_of(_neutron_fluxes.at(elem_id).begin(),
                                  _neutron_fluxes.at(elem_id).end(),
                                  [](double j) { return j == 0; });
  // return true if there is flux, false if there isn't, as the function name
  // implies
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

  // Fetch all subdomain ID's from libmesh
  std::vector<unsigned short> subdomain_ids;
  for (auto &subdomain_id : _mesh.meshSubdomains()) {
    subdomain_ids.push_back(subdomain_id);
  }

  // Fetch all subdomain names using libmesh
  std::vector<SubdomainName> subdomain_names =
      _mesh.getSubdomainNames(subdomain_ids);

  // Check that for all subdomain names there is an associated material
  for (auto &subdomain_name : subdomain_names) {
    if (_mat_definitions.find(subdomain_name) == _mat_definitions.end()) {
      mooseError("Block " + subdomain_name +
                 " does not have a corresponding FISPACT material defined");
    }
  }
}

void FispactProblem::convertGammaEvToCount(
    fp::InputData &input, const std::vector<double> &photon_spectra,
    const std::vector<double> &photon_flux_bins,
    std::vector<double> &photons_per_cc_per_s) {

  // Reserve memory for photons per cc per s vector
  photons_per_cc_per_s.resize(photon_spectra.size());

  // Get inventory density and mass
  double inv_density = input.getDensity();
  double inv_mass = input.getMassTotal();

  // Convert from MeV/s to per cc per s for each bin
  for (int i = 0; i < photon_spectra.size(); i++) {
    double bin_energy = photon_flux_bins[i] +
                        ((photon_flux_bins[i + 1] - photon_flux_bins[i]) / 2);

    // per_cc_per_s = MeV/s * (inventory_density/(inventory_mass *
    // energy_bin_midpoint))
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

  // Attempt to retrieve element strength from map, otherwise catch exception
  // and give a useful error
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
  // Get number of active local elements
  int n_local_elem = _mesh.getMesh().n_active_local_elem();

  unsigned long photon_flux_map_size =
      ((24 * sizeof(double)) + sizeof(int)) * n_local_elem;

  unsigned long element_strengths_map_size =
      (sizeof(int) + sizeof(double)) * n_local_elem;

  unsigned long memory_size = photon_flux_map_size +
                              element_strengths_map_size + (sizeof(int) * 2) +
                              (sizeof(double) * 2);
  // Really naive way of doing this, but currently giving a 20% buffer to
  // account for the memory space required by Boost allocators and such
  return memory_size * 2;
}
