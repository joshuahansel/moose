[GlobalParams]
  initial_T = 393.15
  initial_vel = 0.0372
  f = 0
  fp = fp

  scaling_factor_1phase = '1e-2 1e-2 1e-5'

  closures = simple_closures
[]

[FluidProperties]
  [fp]
    type = StiffenedGasFluidProperties
    gamma = 2.35
    q = -1167e3
    q_prime = 0
    p_inf = 1.e9
    cv = 1816
  []
[]

[Closures]
  [simple_closures]
    type = Closures1PhaseSimple
  []
[]

[Functions]
  [pump_head_fn]
    type = PiecewiseLinear
    x = '0  0.5'
    y = '0  1  '
  []
[]

[ChainControls]
  [pump_head_ctrl_fn]
    type = GetFunctionValueChainControl
    function = pump_head_fn
  []
  [pump_head_ctrl]
    type = SetRealValueChainControl
    parameter = component::pump/head
    value = pump_head_ctrl_fn:value
  []
[]

[Postprocessors]
  [pump_head]
    type = RealComponentParameterValuePostprocessor
    component = pump
    parameter = head
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
  dt = 0.1
  num_steps = 10
  abort_on_solve_fail = true

  solve_type = 'PJFNK'
  line_search = 'basic'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-6
  nl_max_its = 15
  l_tol = 1e-4
  l_max_its = 10

  [Quadrature]
    type = GAUSS
    order = SECOND
  []
[]

[Outputs]
  [out]
    type = CSV
    show = 'pump_head'
  []
  print_linear_residuals = false
[]
