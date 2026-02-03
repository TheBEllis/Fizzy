[Mesh]
  type = FileMesh
  file = "./geometry/cube_hex.e"
[]

[Problem]
  type = OpenMCCellAverageProblem
  cell_level = 0
  source_strength = 1
  xml_directory = './photons_xml/'
  [Tallies]
    [photon_flux]
      type = MeshTally
      filters = 'Photon'
      name = 'photon_flux'
      score = 'flux'
    []
  []

  [Filters]
    [Photon]
      type = ParticleFilter
      particles = 'photon' 
    []
  []
[]

[Executioner]
  type = Steady 
[]

