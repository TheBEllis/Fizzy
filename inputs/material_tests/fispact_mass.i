[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = "../../geometry/cube.e"
  []
[]

[Problem]
  type = FispactProblem 
  fispact_nuclear_data_path = '/home/bill/Projects/FISPACT/nuclear_data/'

  neutron_flux_file = './../statepoint_neutrons.10.h5'
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
    material_type = MASS 
    nuclides = 'Fe Mo'
    nuclide_fraction = '0.7 0.3'
    block = 'steel' 
    density = 5
  []

[]

[Executioner]
  type = Steady
[]
