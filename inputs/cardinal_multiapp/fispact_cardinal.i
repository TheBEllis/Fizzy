[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = "../../example_geom/cube.e"
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
  neutron_flux_hdf5_path = 'tallies/tally 1/results'
  neutron_bin_type = 'G1102'
  write_photon_spectra = True
  photon_spectra_filename = 'photon_spectra_out'
  comm_photon_spectra = True

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

[MultiApps]
  [Cardinal]
    type = FullSolveMultiApp
    execute_on = timestep_end
    app_type = 'CardinalApp'
    input_files = '/home/bill/Projects/cardinal-fispact/inputs/compiled_source/card_input.i'
    library_path = '/home/bill/Projects/cardinal-fispact/lib/'
    library_name = 'libcardinal-opt.la'
  []
[]

[Transfers]
  [Photon_flux_from_cardinal]
    type = MultiAppCopyTransfer
    source_variable =  'photon_flux_photon'
    variable = Photon_flux
    from_multi_app='Cardinal'
  [../]
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = True
[]
