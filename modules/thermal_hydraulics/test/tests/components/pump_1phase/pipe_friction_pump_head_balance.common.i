# This test balances the pipe friction pressure drop with the pump head pressure rise and runs to steady state.

[GlobalParams]
  initial_T = 393.15
  initial_vel = 0.0
  A = 0.567
  fp = fp
  scaling_factor_1phase = '0.04 0.04 0.04e-5'
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
  dt = 1.e-3
  num_steps = 38
  abort_on_solve_fail = true
  solve_type = 'PJFNK'
  line_search = 'basic'
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
  [out_x]
    type = Exodus
    show = 'p T vel'
    output_dimension = 3 # classic components write 3D meshes
  []
  velocity_as_vector = false
[]
