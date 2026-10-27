[GlobalParams]
  gravity_vector = '0 0 0'

  initial_p = 1e6
  initial_T = 353.1
  initial_vel = 0.0

  2nd_order_mesh = true

  closures = simple_closures
[]

[Closures]
  [simple_closures]
    type = Closures1PhaseSimple
  []
[]

[SolidProperties]
  [hs-mat]
    type = ThermalFunctionSolidProperties
    k = 1
    cp = 1
    rho = 1
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
  dt = 0.1
  dtmin = 1.e-7

  solve_type = 'PJFNK'
  line_search = 'basic'
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-6
  nl_max_its = 30

  l_tol = 1e-3
  l_max_its = 100

  start_time = 0.0
  end_time = 4.0

  [Quadrature]
    type = TRAP
    order = FIRST
  []
[]

[Outputs]
  [out]
    type = Exodus
  []
[]
