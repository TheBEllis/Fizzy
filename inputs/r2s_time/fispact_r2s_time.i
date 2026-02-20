[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = "../../geometry/cube.e"
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
  neutron_flux_file = './statepoint_neutrons.10.h5'
  neutron_flux_tally_id = 2
  neutron_bin_structure = 1102 

  write_photon_flux = True
  photon_flux_filename = 'photon_spectra_out'

  comm_photon_flux = True

  molar_mass_data = '../../molar_masses.h5'

  fispact_schedule_uo = 'Schedule'
  fispact_nuclear_data_uo = 'endf_nuclear_data'
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
    nuclides = "B10 B11 C12 C13 Si28 V50 V51 Cr50 Cr52 Cr53 Cr54 Mn55 Fe54 Fe56 Fe57 Fe58 Co59 Ni58 Ni60 Ni61 Ni62 Ni64 Mo92 Mo94 Mo95 Mo96 Mo97 Mo98 Mo100 Cu63 Cu65"
    nuclide_fraction = "0.000693 0.002807 0.039520273803427465 0.000479726196572539 0.47 0.0003921925373160739 0.15960780746268394 0.73 14.07 1.6 0.4 1.14 3.95 62.51 1.46 0.2 0.14 7.31 2.79 0.12 0.38 0.1 0.2974760972643644 0.1906159496693566 0.33284984504744625 0.3533051306488447 0.20516502475964524 0.5254922939545649 0.21509565865577795 0.06 0.03"
    block = 1
    density = 8
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

  [cendl_nuclear_data]
    type = FispactNuclearDataPaths
    base_path = "/Projects/FispactNuclearData/"
    ND_IND_NUC_KEY = "decay2020/decay_2020_index.txt"
    #ND_XS_ENDF_KEY = "CENDL32data/gendf-1102"
    ND_XS_ENDFB_KEY = "CENDL-n.bin"
    ND_FY_ENDF_KEY = "GEFY61data/gefy61_nfy"
    ND_SF_ENDF_KEY = "GEFY61data/gefy61_sfy"
    ND_PROB_TAB_KEY = "CENDL32data/tp-1102-294"
    ND_DK_ENDF_KEY = "decay2020/decay_2020"
    ND_ABSORP_KEY = "decay/abs_2012"
  []
  [tendl_nuclear_data]
    type = FispactNuclearDataPaths
    base_path = "/Projects/FispactNuclearData/"
    ND_IND_NUC_KEY = "TENDL2017data/tendl17_decay12_index"
    #ND_XS_ENDF_KEY = "TENDL2017data/tal2017-n/gxs-709"
    ND_XS_ENDFB_KEY = "TENDL-n.bin"
    ND_FY_ENDF_KEY = "GEFY61data/gefy61_nfy"
    ND_SF_ENDF_KEY = "GEFY61data/gefy61_sfy"
    ND_PROB_TAB_KEY = "TENDL2017data/tal2017-n/tp-709-294"
    ND_DK_ENDF_KEY = "decay2020/decay_2020"
    ND_ABSORP_KEY = "decay/abs_2012"
  []
[]

[MultiApps]
  [photons]
    type = TransientMultiApp
    execute_on = timestep_end
    app_type = 'CardinalApp'
    input_files = '/Projects/Fizzy/inputs/r2s_time/photons_input.i'
    library_path = '/Projects/cardinal/lib/'
    library_name = 'libcardinal-opt.la'
  []
[]

[Transfers]
  [Photon_flux_from_cardinal]
    type = MultiAppCopyTransfer
    source_variable =  'photon_flux_photon'
    variable = Photon_flux
    from_multi_app='photons'
  [../]
[]

[Times]
  [FizzyTimes]
    type = FispactScheduleTimes
    FispactScheduleName = Schedule
#    FispactScheduleTimeIndices = '0 1'

  []
[]

[Executioner]

  type = Transient

  [TimeStepper]
    type = TimeSequenceFromTimes
    times = FizzyTimes 
    use_last_t_for_end_time = True
  []
[]

[Outputs]
  exodus = True
[]






