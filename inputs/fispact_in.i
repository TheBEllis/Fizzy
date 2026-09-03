[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = "../example_geom/cube.e"
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
  fispact_nuclear_data_path = '/home/bill/Projects/FISPACT/nuclear_data/'

  neutron_flux_file = '/home/bill/Projects/PyFIS/cube_example/neutron_cube/statepoint.2.h5'
  neutron_flux_tally_id = 1
  neutron_bin_type = 'G1102'
  write_photon_flux = True
  photon_flux_filename = 'photon_spectra_test.h5'
  comm_photon_flux = True

  fispact_schedule_uo = 'Schedule'

  read_materials_from_xml = True
  materials_xml_file = '/home/bill/Projects/fizzy/materials.xml'
[]

[UserObjects]
  [Schedule]
    type = FispactSchedule
    times = '300 30 30 30 30 30'
    flux_schedule = '1e10 0 0 0 0 0'  
  []
[]


[Executioner]
  type = Steady
[]

[Outputs]
  exodus = True
[]
