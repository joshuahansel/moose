T_initial = 300

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  end_time = 10000
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.0
    optimal_iterations = 5
    iteration_window = 0
    growth_factor = 1.2
    cutback_factor = 0.8
  []
  dtmin = 1e-3

  steady_state_detection = true
  steady_state_tolerance = 1e-8

  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-8
  nl_max_its = 20

  l_tol = 1e-3
  l_max_its = 10
[]

[Outputs]
  exodus = true
  csv = true
  perf_graph = true
[]
