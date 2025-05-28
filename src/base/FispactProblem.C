#include "FispactProblem.h"

InputParameters FispactProblem::validParams()
{
    InputParameters params = ExternalProblem::validParams();

    // Suppress parameters we are not using
    params.suppressParameter<bool>("allow_invalid_solution");
    params.suppressParameter<bool>("boundary_restricted_elem_integrity_check");
    params.suppressParameter<bool>("boundary_restricted_node_integrity_check");
    params.suppressParameter<bool>("check_uo_aux_state");
    params.suppressParameter<bool>("error_on_jacobian_nonzero_reallocation");
    params.suppressParameter<std::vector<std::vector<TagName>>>("extra_tag_matrices");
    params.suppressParameter<std::vector<TagName>>("extra_tag_solutions");
    params.suppressParameter<std::vector<std::vector<TagName>>>("extra_tag_vectors");
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
    params.addRequiredParam<std::string>("neutron_flux_file", "HDF5 file storing neutron flux");
    params.addRequiredParam<std::string>("neutron_flux_hdf5_path", "Path within HDF5 file for the vector storing neutron flux");


    params.addRequiredParam<std::string>("fispact_nuclear_data_path", "Path to FISPACT nuclear data");


    params.addParam<bool>("read_materials_from_xml", false, "Parameter determining whether user wishes to read materaial nuclide compositions from openMC XML file");
    params.addParam<std::string>("materials_xml_file", "materials.xml", "Path and name of material xml file user wishes to use.");

    return params;
}


FispactProblem::FispactProblem(const InputParameters& params) :
    ExternalProblem(params),
    _fp_monitor(fispactLogName()),
    _fp_nuclear_data(_fp_monitor),
    _neutron_flux_filename(getParam<std::string>('neutron_flux_file')),
    _neutron_flux_hdf5_path(getParam<std::string>('neutron_flux_hdf5_path'))
{
    // Initialise FISPACT
    fp::GlobalInitialise(monitor);

    // Get the path to our nuclear data
    std::string fp_nuclear_data_path = getParam<std::string>("fispact_nuclear_data_path");

    // Set nuclear data paths
    setNuclearData(fp_nuclear_data_path);

    // Read neutron flux from HDF5
    readNeutronFluxFromHDF5(neutron_flux_file, neutron_flux_hdf5_path);

    if(getParam<bool>("read_materials_from_xml"))
    {
        
    }
}


void FispactProblem::externalSolve()
{
    // Set up fispact input data
    fp::InputData fispact_input(_fp_monitor);
    fp::OutputData fispact_output(_fp_monitor);
    setFispactInputData(_fp_monitor, fispact_input, fispact_output);

    // Check if there is a
    bool is_zero_flux = std::all_of(neutron_fluxes[i].begin(), neutron_fluxes[i].end(), [](double j) { return j==0; });
    if(!is_flux_zero)
    {
        setFispactInputData(_fp_monitor, fispact_input, fispact_output);
    }
}    

std::vector<double> FispactProblem::readNeutronFluxFromH5(const std::string& filename)
{
    std::vector<double> neutron_flux;

    return neutron_flux;    
}

std::string FispactProblem::fispactLogName()
{
    std::string log_name = "FISPACT_app_" + std::to_string(processor_id()) + ".log";
    return log_name;
}

void FispactProblem::setNuclearData(std::string nd_base_path)
{
    fp::io::NuclearDataReader nd_reader(_fp_monitor);

    nd_reader.setPath(FISPACT_ND_IND_NUC_KEY, nd_base_path + "decay2020/decay_2020_index.txt");
    
    nd_reader.setPath(FISPACT_ND_XS_ENDF_KEY, nd_base_path + "TENDL2021data/gendf-1102");
    nd_reader.setPath(FISPACT_ND_PROB_TAB_KEY, nd_base_path + "TENDL2021data/tp-1102-294");

    nd_reader.setPath(FISPACT_ND_FY_ENDF_KEY, nd_base_path + "GEFY61data/gefy61_nfy");
    nd_reader.setPath(FISPACT_ND_SF_ENDF_KEY, nd_base_path + "GEFY61data/gefy61_sfy");

    nd_reader.setPath(FISPACT_ND_DK_ENDF_KEY, nd_base_path + "decay2020/decay_2020");
    nd_reader.setPath(FISPACT_ND_ABSORP_KEY, nd_base_path + "decay/abs_2012");

    nd_reader.setPath(FISPACT_ND_HAZARDS_KEY, nd_base_path + "/decay/hazards_2012");
    nd_reader.setPath(FISPACT_ND_CLEAR_KEY, nd_base_path + "/decay/clear_2012");
    nd_reader.setPath(FISPACT_ND_A2DATA_KEY, nd_base_path + "/decay/a2_2012");
    
    nd_reader.load(_fp_nuclear_data, &FispactProblem::load_callback);
}


std::vector<std::vector<double>> FispactProblem::readNeutronFluxFromHDF5(std::string filename, std::string tally_dir)
{
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
    hsize_t memspace_dimensions_realizations[3] = {1,1,1};
    H5::DataSpace memspace_realizations(1, memspace_dimensions_realizations);
    dataset_realizations.read(&n_realizations, H5::PredType::STD_I32LE, memspace_realizations, dspace_realizations);

    // Set up vector of vectors to store neutron fluxes
    std::vector<std::vector<double>> neutron_fluxes(tally_array_dims[0]/NUM_NEUTRON_BINS, std::vector<double>(NUM_NEUTRON_BINS, 0));
    for(int i = 0; i < neutron_fluxes.size(); i++)
    {
        // Set up std::vector to store neutron flux data
        std::vector<double> neutron_flux_data(NUM_NEUTRON_BINS, 0.0);

        // Set up counts and offsets for selecting hyperslab of tally array 
        hsize_t dataCount[3] = {NUM_NEUTRON_BINS, 1, 1};
        hsize_t dataOffset[3] = {(i * NUM_NEUTRON_BINS),0,0};

        // Set up memory space for reading tally results
        hsize_t arr_len[3] = {NUM_NEUTRON_BINS,1,1};
        H5::DataSpace memspace_tally (1, arr_len);
        dspace_tally.selectHyperslab(H5S_SELECT_SET, dataCount, dataOffset);

        // Read in neutron flux tally results
        dataset_tally.read(neutron_flux_data.data(), H5::PredType::IEEE_F64LE, memspace_tally, dspace_tally);

        for(auto& neutron_flux: neutron_flux_data)
        {
            neutron_flux = neutron_flux / n_realizations;
        }
        neutron_fluxes[i] = neutron_flux_data;  
    }

    int world_size = comm().size();
    // for now manual parallelisation, later switch to phdf5
    int count = neutron_fluxes.size() / world_size;
    int remainder = neutron_fluxes.size() % world_size;
    int start, stop;

    if (world_rank < remainder) {
        // The first 'remainder' ranks get 'count + 1' tasks each
        start = world_rank * (count + 1);
        stop = start + count;
    } else {
        // The remaining 'size - remainder' ranks get 'count' task each
        start = world_rank * count + remainder;
        stop = start + (count - 1);
    }

    neutron_fluxes.resize(stop - start);
    for(int i = start; i < stop; i++)
    {
        _neutron_fluxes.insert(i - start ,neutron_fluxes[i]);
    }
}


void FispactProblem::setFispactInputData(fp::Monitor& monitor, fp::InputData& input, std::string material, std::vector<double> neutron_flux, std::vector<double>& bins)
{
    input.setFlux(bins, neutron_flux);
    input.setFluxWallLoading(1.0);
    input.setFluxName("neutrons");

    input.setDensity(material_def.getDensity());
    input.setAtomsThreshold(1.0e3);

    // Volume read in is in cm ^ 3, so we need to scale by 1e-6
    double total_mass = material_def.getDensity() * volume * 1e-6;
    input.setMassTotal(total_mass);

    std::vector<int> atomic_numbers;
    std::vector<double> percent;
    for(auto& element_name_ao_pair: material_def.getAtomicComposition())
    {
        std::string element_name = element_name_ao_pair.first;
        element_name.erase(std::remove_if(element_name.begin(), element_name.end(), [](unsigned char c) { 
            return std::isdigit(c); 
        }), element_name.end()); 
        // std::cout << element_name_ao_pair.first << " " << fp::util::GetAtomicNumberFromElementName(monitor, "C13") << std::endl;
        atomic_numbers.push_back(fp::util::GetAtomicNumberFromElementName(monitor, element_name));
        percent.push_back(element_name_ao_pair.second);
    }
    
    input.setMass(atomic_numbers, percent);

    std::vector<double> irradiationtime = {5.0 * FISPACT_MIN_TO_SEC};
    std::vector<double> fluxamp = {1.116e10};
    input.setSchedule(irradiationtime, fluxamp);

    // std::vector<double> cooltimes = {36, 15, 16, 15, 15 ,26, 33, 36, 53, 66, 66, 97};
    std::vector<double> cooltimes = {36, 15, 16};

    for(double time : cooltimes)
    {
        input.appendSchedule(time, 0.0);
    }
}     

