[Mesh]
  type = FileMesh
  file = "../../geometry/cube.e"
[]

[Problem]
  type = OpenMCCellAverageProblem
  cell_level = 0
  source_strength = 1
  xml_directory = './neutrons_xml'
  [Tallies]
    [neutron_flux]
      type = MeshTally
      filters = 'neutron_filter energy_filter'
      name = 'neutron_flux'
      score = 'flux'
    []
  []

  [Filters]
    [neutron_filter]
      type = ParticleFilter
      particles = 'neutron' 
    []
    [energy_filter]
      type = EnergyFilter
      group_structure = UKAEA_1102
    []
  []
[]

[Executioner]
  type = Steady 
[]

[Outputs]
  exodus = True
[]
