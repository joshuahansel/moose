[FluidProperties]
  [fp]
    type = IdealGasFluidProperties
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

[Convergence]
  [components_conv]
    type = ComponentsConvergence
    max_iterations = 15
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2

  start_time = 0
  dt = 1
  num_steps = 1

  solve_type = NEWTON
  nonlinear_convergence = components_conv
  l_tol = 1e-3
  l_max_its = 10

  verbose = true
[]
