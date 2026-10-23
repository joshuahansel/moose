[GlobalParams]
  gravity_vector = '0 -9.81 0'

  initial_T = 310
  initial_p = 1e5
  initial_vel = 0

  scaling_factor_1phase = '1e0 1e-2 1e-4'

  closures = simple_closures
[]

[FluidProperties]
  [fp]
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

[Problem]
  solve = false
[]

[Executioner]
  type = Transient
  scheme = 'bdf2'

  start_time = 0
  dt = 1e-2
  num_steps = 1
  abort_on_solve_fail = true

  solve_type = 'PJFNK'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-8
  nl_max_its = 30

  l_tol = 1e-2
  l_max_its = 30
[]

[Outputs]
  exodus = true
  show = 'A'
[]
