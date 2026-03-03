[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = "../../geometry/cube.e"
  []
[]

[AuxVariables]
  [Photon_flux]
    family = MONOMIAL
    order = CONSTANT
  []
[]


[Problem]
  type = FispactProblem 
  neutron_flux_file = './statepoint_neutrons.10.h5'
  neutron_flux_tally_id = 2
  neutron_bin_structure = 1102
  output_inventory_time = 450

  read_materials_from_xml = True
  materials_xml_file = "materials.xml"

  write_photon_flux = True 
  photon_flux_filename = 'photon_spectra_out'

  comm_photon_flux = True

  molar_mass_data = '../../molar_masses.h5'

  fispact_schedule_uo = 'Schedule'
  fispact_nuclear_data_uo = 'endf_nuclear_data'
[]

[UserObjects]
  [Schedule]
    type = FispactSchedule
    times = '300 30 30 30 30 30'
    flux_amplitude = '1e10 0 0 0 0 0'  
  []

  [endf_nuclear_data]
    type = FispactNuclearDataPaths
    base_path = "/Projects/FispactNuclearData/"
    ND_IND_NUC_KEY = "ENDFB80data/endfb80_index"
    #ND_XS_ENDF_KEY = "ENDFB80data/endfb80-n/gxs-709"
    ND_XS_ENDFB_KEY = "endfb80-n.bin"
    ND_FY_ENDF_KEY = "ENDFB80data/endfb80-n/endfb80nfy"
    ND_SF_ENDF_KEY = "ENDFB80data/endfb80-n/endfb80sfy"
    ND_DK_ENDF_KEY = "ENDFB80data/decay"
    ND_ABSORP_KEY = "decay/abs_2012"
  []
[]

[Executioner]
  type = Steady
[]







