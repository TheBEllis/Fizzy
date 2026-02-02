[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = "../../example_geom/cube.e"
  []
[]

[Problem]
  type = FispactProblem 

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
    ND_IND_NUC_KEY = "ENDFB80data/endfb80_index.txt"
    ND_XS_ENDF_KEY = "ENDFB80data/endfb80-n/gxs-709"
    ND_FY_ENDF_KEY = "ENDFB80data/endfb80-n/endfb80nfy"
    ND_SF_ENDF_KEY = "ENDFB80data/endfb80-n/endfb80sfy"
    ND_DK_ENDF_KEY = "ENDFB80data/decay"
    ND_ABSORP_KEY = "decay/abs_2012"
  []
[]

[Executioner]
  type = Steady
[]
