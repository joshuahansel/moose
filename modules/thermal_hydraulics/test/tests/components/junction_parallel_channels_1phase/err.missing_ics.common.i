[GlobalParams]
  gravity_vector = '0 0 0'

  A = 1e-4
  f = 0
  fp = fp

  closures = simple_closures
[]

[FluidProperties]
  [fp]
    type = StiffenedGasFluidProperties
    gamma = 1.4
    cv = 725
    p_inf = 0
    q = 0
    q_prime = 0
  []
[]

[Closures]
  [simple_closures]
    type = Closures1PhaseSimple
  []
[]

[Executioner]
  type = Transient
  abort_on_solve_fail = true
[]
