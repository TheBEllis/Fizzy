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
  neutron_flux_file = './statepoint.10.h5'
  neutron_flux_tally_id = 2
  neutron_bin_type = 'G1102'

  write_photon_flux = False
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
    nuclides = "C13 C12 Mn55 P31 S34 S32 S33 S36 Si28 Si29 Si30 Cr53 Cr50 Cr52 Cr54 Ni62 Ni60 Ni64 Ni58 Ni61 Fe56 Fe58 Fe57 Fe54"
    nuclide_fraction = "1.5132548e-05 0.00135086745 0.01991 0.000795 2.14834688e-05 0.000486608589 3.8332928e-06 7.46496e-08 0.0179598856 0.000911951747 0.000601162667 0.0189874635 0.00868335215 0.167449803 0.00472638155 0.0033866271 0.0244346846 0.00086247408 0.0634340554 0.00106215882 0.610087944 0.00187506594 0.0140895912 0.0388643986"
    block = 1
    density = 8
  []

  [endf_nuclear_data]
    type = FispactNuclearDataPaths
    base_path = "/Projects/FispactNuclearData/"
    ND_IND_NUC_KEY = "ENDFB80data/endfb80_index"
    ND_XS_ENDF_KEY = "ENDFB80data/endfb80-n/gxs-709"
    ND_FY_ENDF_KEY = "ENDFB80data/endfb80-n/endfb80nfy"
    ND_SF_ENDF_KEY = "ENDFB80data/endfb80-n/endfb80sfy"
    ND_DK_ENDF_KEY = "ENDFB80data/decay"
    ND_ABSORP_KEY = "decay/abs_2012"
  []
[]

[MultiApps]
  [photons]
    type = FullSolveMultiApp
    execute_on = timestep_end
    app_type = 'CardinalApp'
    input_files = '/Projects/Fizzy/inputs/r2s_multiapp/photons_input.i'
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

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = True
[]






