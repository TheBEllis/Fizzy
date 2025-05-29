[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = "../librti-tets-scaled.e"
  []
[]

[Problem]
  type = FispactProblem 
  fispact_nuclear_data_path = '/home/bill/Projects/FISPACT/nuclear_data/'

  neutron_flux_file = '../'
  neutron_flux_hdf5_path = '../tallies/tally 1/results'
  neutron_bin_type = 'G1102'

  read_material_xml_data = true
  materials_xml_file = 'materials.xml '
[]


[Executioner]
  type = Steady
[]
