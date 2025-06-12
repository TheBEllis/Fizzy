[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = "../librti-tets-scaled.e"
  []
[]

[Problem]
  type = FispactProblem 
  fispact_nuclear_data_path = '/home/bill/Projects/FISPACT/nuclear_data/'

  neutron_flux_file = '/home/bill/Projects/PyFIS/statepoint.3.h5'
  neutron_flux_hdf5_path = 'tallies/tally 1/results'
  neutron_bin_type = 'G1102'
  
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
