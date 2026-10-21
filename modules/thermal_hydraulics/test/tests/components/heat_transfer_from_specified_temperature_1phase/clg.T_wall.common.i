[GlobalParams]
  initial_p = 0.1e6
  initial_vel = 0
  initial_T = 300
  scaling_factor_1phase = '1e+0 1e-2 1e-4'
  closures = simple_closures
[]

[FluidProperties]
  [eos]
    type = StiffenedGasFluidProperties
    gamma = 2.35
    cv = 1816.0
    q = -1.167e6
    p_inf = 1.0e9
    q_prime = 0
  []
[]

[Closures]
  [simple_closures]
    type = Closures1PhaseSimple
  []
[]

[Preconditioning]
  [SMP_PJFNK]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  scheme = 'bdf2'

  start_time = 0
  dt = 0.05
  num_steps = 20
  abort_on_solve_fail = true

  solve_type = 'newton'
  line_search = 'basic'
  petsc_options_iname = '-pc_type'
  petsc_options_value = ' lu'
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-8
  nl_max_its = 10

  l_tol = 1e-3
  l_max_its = 30
[]

[Outputs]
  csv = true
[]

[Functions]
  [T_wall_fn]
    type = PiecewiseLinear
    x = '0 1'
    y = '310 320'
  []
[]

[Postprocessors]
  [T_wall]
    type = ElementAverageValue
    block = 'pipe1'
    variable = T_wall
  []
[]
