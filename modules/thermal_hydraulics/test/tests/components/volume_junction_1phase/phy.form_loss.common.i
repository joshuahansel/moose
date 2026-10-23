# This test measures the pressure drop across the volume junction with K=1.

A = 0.1

[GlobalParams]
  gravity_vector = '0 0 0'

  scaling_factor_1phase = '1 1 1e-5'
  scaling_factor_rhoV = 1
  scaling_factor_rhouV = 1
  scaling_factor_rhovV = 1
  scaling_factor_rhowV = 1
  scaling_factor_rhoEV = 1e-5

  initial_T = 300
  initial_p = 1e5
  initial_vel = 1

  n_elems = 20
  length = 1

  f = 0

  fp = fp
  closures = simple_closures
[]

[FluidProperties]
  [fp]
    type = StiffenedGasFluidProperties
    gamma = 1.4
    cv = 725
    q = 0
    q_prime = 0
    p_inf = 0
  []
[]

[Closures]
  [simple_closures]
    type = Closures1PhaseSimple
  []
[]

[Functions]
  [K_fn]
    type = TimeRampFunction
    initial_value = 0
    initial_time = 2
    ramp_duration = 5
    final_value = 1
  []
[]

[ChainControls]
  active = ''
  [K_crtl_fn]
    type = GetFunctionValueChainControl
    execute_on = 'INITIAL TIMESTEP_BEGIN'
    function = K_fn
  []
  [K_crtl]
    type = SetRealValueChainControl
    execute_on = 'INITIAL TIMESTEP_BEGIN'
    parameter = component::junction/K
    value = K_crtl_fn:value
  []
[]

[Postprocessors]
  [pJ_in]
    type = SideAverageValue
    variable = p
    boundary = pipe1:out
  []

  [pJ_out]
    type = SideAverageValue
    variable = p
    boundary = pipe2:in
  []

  [dpJ]
    type = DifferencePostprocessor
    value1 = pJ_in
    value2 = pJ_out
  []
[]

[Preconditioning]
  [pc]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  scheme = 'bdf2'

  start_time = 0
  end_time = 20
  dt = 0.5
  abort_on_solve_fail = true

  solve_type = 'NEWTON'
  nl_rel_tol = 0
  nl_abs_tol = 1e-8
  nl_max_its = 15

  l_tol = 1e-3
  l_max_its = 10

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
  execute_on = 'final'
  show = 'dpJ'
[]
