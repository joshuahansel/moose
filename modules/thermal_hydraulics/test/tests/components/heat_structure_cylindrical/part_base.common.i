[Functions]
  [power_profile_fn]
    type = ParsedFunction
    expression = '1.570796326794897 * sin(x / 3.6576 * pi)'
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
  dt = 2
  num_steps = 10

  solve_type = 'NEWTON'
  nonlinear_convergence = nl_conv

  l_tol = 1e-4
  l_max_its = 300
[]

[Outputs]
  file_base = transient
  [exodus]
    type = Exodus
    output_dimension = 3 # classic components write 3D meshes
  []
  [console]
    type = Console
    execute_scalars_on = none
  []
[]
