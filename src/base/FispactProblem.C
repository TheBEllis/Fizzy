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

    return params;
}


FispactProblem::FispactProblem(const InputParameters& params) :
    ExternalProblem(params),
    _fp_monitor(fispactLogName()),
    _fp_nuclear_data(_fp_monitor)
{
    // Get the path to our nuclear data
    std::string fp_nuclear_data_path = getParam<std::string>("fispact_nuclear_data_path");


    setNuclearData(fp_nuclear_data_path);
}


void FispactProblem::externalSolve()
{

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