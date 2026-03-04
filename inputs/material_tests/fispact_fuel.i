[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = "../../geometry/cube.e"
  []
[]

[Problem]
  type = FispactProblem 

  write_photon_flux = True
  photon_flux_filename = 'photon_spectra_out'

  molar_mass_data = '../../molar_masses.h5'

  fispact_schedule_uo = 'Schedule'
  fispact_nuclear_data_uo = 'endf_nuclear_data'
  fispact_input_flux_uo = 'InputFlux'
[]

[UserObjects]

  [InputFlux]
    type = OpenMCFluxInput
    statepoint_filename = './statepoint_neutrons.10.h5'
    energy_filter_id = 2
    flux_tally_id = 2
    wall_loading = 10
  []

  [Schedule]
    type = FispactSchedule
    times = '300 30 30 30 30 30'
    flux_amplitude = '1e10 0 0 0 0 0'
  []

  [steel]
    type = FispactMaterial
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
[]

[Executioner]
  type = Steady
[]




