[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = "../../geometry/cube.e"
  []
[]

[Problem]
  type = FispactProblem 

  neutron_flux_file = './statepoint_neutrons.10.h5'
  neutron_flux_tally_id = 2
  neutron_bin_structure = 1102

  write_photon_flux = True
  photon_flux_filename = 'photon_spectra_out'

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

  [steel]
    type = FispactMaterial
    material_type = MASS 
    nuclides = 'Fe Mo'
    nuclide_fraction = '0.7 0.3'
    block = 'steel' 
    density = 5
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
