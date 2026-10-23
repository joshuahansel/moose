[GlobalParams]
  initial_T = 300
  initial_p = 1e5
  initial_vel = 0

  initial_vel_x = 0
  initial_vel_y = 0
  initial_vel_z = 0

  scaling_factor_1phase = '1 1 1'
  scaling_factor_rhoV  = 1
  scaling_factor_rhouV = 1
  scaling_factor_rhovV = 1
  scaling_factor_rhowV = 1
  scaling_factor_rhoEV = 1

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
  [SMP_PJFNK]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  scheme = 'bdf2'

  start_time = 0
  num_steps = 10
  dt = 1
  abort_on_solve_fail = true

  solve_type = 'PJFNK'
  line_search = 'basic'

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-7
  nl_max_its = 10

  l_tol = 1e-3
  l_max_its = 100

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  [Quadrature]
    type = gauss
    order = second
  []
[]

[Outputs]
  [out]
    type = Exodus
    show = 'rhouA p'
    output_dimension = 3 # classic components write 3D meshes
    append_object_name = false
    execute_on = 'initial final'
  []
[]
