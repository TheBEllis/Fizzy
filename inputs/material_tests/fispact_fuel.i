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

  output_inventory_time = 1e5

  fispact_schedule_uo = 'Schedule'
  fispact_nuclear_data_uo = 'endf_nuclear_data'
  fispact_input_flux_uo = 'InputFlux'

  atol = 1e-1
  rtol = 1e-1
[]

[AuxVariables]
  [B10_atoms]
    family = MONOMIAL
    order = CONSTANT
  []
  [Nb93_atoms]
    family = MONOMIAL
    order = CONSTANT
  []
  [atoms]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [B10_atoms_kernel]
    type = FispactNuclideKernel
    variable = B10_atoms
    nuclide = "B10"
    metric = "ATOMS"
 []
 [B12_atoms_kernel]
    type = FispactNuclideKernel
    variable = Nb93_atoms 
    nuclide = "Nb93"
    metric = "ATOMS"
 []
 [atoms_kernel]
    type = FispactElementKernel
    variable = atoms 
    metric = "INVENTORY_TOTAL_ATOMS"
  []
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
    times = '1e5 1e4 1e4 2e5 2e5 2e5 5e5 5e5 2e5'
    flux_amplitude = '1e10 0 0 0 0 0 0 0 0'
  []

  [steel]
    type = FispactMaterial
    material_type = FUEL
    nuclides = "B10 B11 C12 C13 Si28 V50 V51 Cr50 Cr52 Cr53 Cr54 Mn55 Fe54 Fe56 Fe57 Fe58 Co59 Ni58 Ni60 Ni61 Ni62 Ni64 Mo92 Mo94 Mo95 Mo96 Mo97 Mo98 Mo100 Cu63 Cu65"
    nuclide_fraction = "0.00000693 0.0000281 0.000395203 0.00000479726 0.0047 0.00000392193 0.00159608 0.0073 0.1407 0.016 0.004 0.0114 0.0395 0.6251 0.0146 0.002 0.0014 0.0731 0.0279 0.0012 0.0038 0.001 0.00297476 0.00190616 0.00332850 0.00353305 0.00205165 0.00525492 0.00215096 0.0006 0.0003"
    fraction_type = 'wo'
    density_units = 'g/cm3'
    block = 1
    density = 8
  []

  [endf_nuclear_data]
    type = FispactNuclearDataPaths
    base_path = "/home/bill/Projects/FispactNuclearData/"
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
    base_path = "/home/bill/Projects/FispactNuclearData/"
    ND_IND_NUC_KEY = "decay_2020_index.txt"
    ND_XS_ENDF_KEY = "CENDL32data/gendf-1102"
    #ND_XS_ENDFB_KEY = "CENDL-n.bin"
    ND_FY_ENDF_KEY = "GEFY61data/gefy61_nfy"
    ND_SF_ENDF_KEY = "GEFY61data/gefy61_sfy"
    ND_PROB_TAB_KEY = "CENDL32data/tp-1102-294"
    ND_DK_ENDF_KEY = "decay_2020"
    ND_ABSORP_KEY = "decay/abs_2012"
  []
[]

[Postprocessors]
  [photon_emission]
    type = ElementPhotonEmission
    element_ids = '0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47'
  []

  [nuclide_contribution]
    type = NuclideMetric
    nuclides = "Cr51"
    metric = "DOSE"
  []

[]

[VectorPostprocessors]
  [PhotonEmissionAllBlocks]
  []
[]

[Times]
  [FizzyTimes]
    type = FispactScheduleTimes
    FispactScheduleName = Schedule
#    FispactScheduleTimeIndices = '0 1'

  []
[]

[Executioner]
  #type = Steady 

  type = Transient

  [TimeStepper]
    type = TimeSequenceFromTimes
    times = FizzyTimes 
    use_last_t_for_end_time = True
  []
[]

[Outputs]
  [CSV]
    type = CSV
  []
  exodus = true
[]




