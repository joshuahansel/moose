# Testing that T_solid gets properly projected onto a pipe
# That's why Hw in pipe1 is set to 0, so we do not have any heat exchange
# Note that the pipe and the heat structure have an opposite orientation, which
# is crucial for this test.

[GlobalParams]
  initial_p = 1.e5
  initial_vel = 0.
  initial_T = 300.

  closures = simple_closures
[]

[FluidProperties]
  [eos]
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

[SolidProperties]
  [wall-mat]
    type = ThermalFunctionSolidProperties
    k = 100.0
    rho = 100.0
    cp = 100.0
  []
[]

[Functions]
  [T_init]
    type = ParsedFunction
    expression = '290 + sin((1 - z) * pi * 1.4)'
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
  dt = 1
  abort_on_solve_fail = true

  solve_type = 'NEWTON'
  line_search = 'basic'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-6
  nl_max_its = 20

  l_tol = 1e-3
  l_max_its = 300

  start_time = 0.0
  num_steps = 1
[]

[Outputs]
  [out]
    type = Exodus
    show = 'T_wall T_solid'
  []
  print_linear_residuals = false
[]
