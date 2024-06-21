!include part_base.i
!include part_mono.i
!include part_htc_gap.i

[Variables]
  [temperature]
    initial_condition = ${T_initial}
  []
[]

[Kernels]
  [time_derivative]
    type = ADHeatConductionTimeDerivative
    variable = temperature
  []
  [heat_conduction]
    type = ADHeatConduction
    variable = temperature
  []
[]

[MultiApps]
  [sub]
    type = TransientMultiApp
    app_type = ThermalHydraulicsApp
    input_files = sub_2d.i
    positions = '0 0 0'
    max_procs_per_app = 1
    output_in_position = true
    output_sub_cycles = false
    sub_cycling = false
    catch_up = true
    execute_on = 'TIMESTEP_BEGIN'
  []
[]

[FunctorMaterials]
  [gap_heat_flux_fmat]
    type = ADConvectionHeatFluxFunctorMaterial
    heat_flux_name = gap_heat_flux
    htc = ${htc_gap}
    T_fluid = temperature
    T_solid = T_hp
  []
[]

[BCs]
  [gap_heat_flux_bc]
    type = FunctorNeumannBC
    variable = temperature
    boundary = mono_inner
    functor = gap_heat_flux
  []
[]

[AuxVariables]
  [T_hp]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = ${T_initial}
  []
[]

[UserObjects]
  [T_layered_avg_uo]
    type = NearestPointLayeredSideAverage
    variable = temperature
    boundary = mono_inner
    num_layers = ${n_axial_mono}
    direction = z
    direction_min = 0
    direction_max = ${length_mono}
    points = '0 0 0'
    execute_on = 'INITIAL TIMESTEP_END' # check that TIMESTEP_END executes during FP iteration
  []
[]

[Transfers]
  [gap_heat_rate_to_sub]
    type = MultiAppPostprocessorTransfer
    to_multi_app = sub
    from_postprocessor = gap_heat_rate
    to_postprocessor = gap_heat_rate_main
  []
  [T_to_sub]
    type = MultiAppGeneralFieldUserObjectTransfer
    variable = T_mono
    to_multi_app = sub
    source_user_object = T_layered_avg_uo
    to_boundaries = hp:evap:outer
    use_bounding_boxes = true
    error_on_miss = true
  []
  [T_from_sub]
    type = MultiAppGeneralFieldUserObjectTransfer
    variable = T_hp
    from_multi_app = sub
    source_user_object = T_hp_uo
    to_boundaries = mono_inner
    use_bounding_boxes = false
    greedy_search = true
    use_nearest_app = true
    error_on_miss = true
  []
[]

[Postprocessors]
  [gap_heat_rate]
    type = ADSideIntegralFunctorPostprocessor
    boundary = mono_inner
    functor = gap_heat_flux
    functor_argument = qp
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [dt]
    type = TimestepSize
  []
  [ss_err]
    type = RelativeSolutionDifferenceNorm
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [wall_time]
    type = PerfGraphData
    section_name = 'ThermalHydraulicsTestApp (main)'
    data_type = TOTAL
  []
  [n_fp_iterations]
    type = NumFixedPointIterations
  []
[]

[Executioner]
  fixed_point_max_its = 15
  fixed_point_abs_tol = 1e-8
  fixed_point_rel_tol = 1e-8
[]
