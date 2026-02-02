[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = "../../example_geom/cube.e"
  []
[]

[Problem]
  type = FispactProblem 
  fispact_nuclear_data_path = '/home/bill/Projects/FISPACT/nuclear_data/'
  #fispact_nuclear_data_path = '/home/bill/Projects/FispactNuclearData/ENDFB80data/'

  neutron_flux_file = './../statepoint.2.h5'
  neutron_flux_tally_id = 1
  neutron_bin_type = 'G1102'

  write_photon_flux = True
  photon_flux_filename = 'photon_spectra_out'

  molar_mass_data = '../../molar_masses.h5'

  fispact_schedule_uo = 'Schedule'
[]

[UserObjects]
  [Schedule]
    type = FispactSchedule
    times = '300 30 30 30 30 30'
    flux_schedule = '1e10 0 0 0 0 0'  
  []

  [steel]
    type = FISPACTMaterial
    material_type = FUEL
    nuclides = 'Fe56 Fe57'
    nuclide_fraction = '0.7 0.3'
    block = 'steel' 
    density = 5
  []

  [nuclear_data]
    type = FispactNuclearDataPaths
    base_path = "/Projects/FispactNuclearData/"
    ND_IND_NUC_KEY = "/decay_2020_index.txt"
    ND_XS_ENDF_KEY = "TENDL2021data/gendf-1102"
    ND_PROB_TAB_KEY = "TENDL2021data/tp-1102-294"
    ND_FY_ENDF_KEY = "GEFY61data/gefy61_nfy"
    ND_SF_ENDF_KEY = "GEFY61data/gefy61_sfy"
    ND_DK_ENDF_KEY = "decay_2020"
    ND_ABSORP_KEY = "decay_2012/abs_2012"
    ND_HAZARDS_KEY = "decay_2012/hazards_2012"
    ND_CLEAR_KEY = "decay_2012/clear_2012"
    ND_A2DATA_KEY = "decay_2012/a2_2012"

    execute_on = INITIAL
  []
[]

[Executioner]
  type = Steady
[]
