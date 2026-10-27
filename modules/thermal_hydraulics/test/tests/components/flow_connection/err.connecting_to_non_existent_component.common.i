# Tests that we report an error if users try to connect to a non-existent component

[GlobalParams]
  initial_p = 1e5
  initial_T = 300
  initial_vel = 0
  closures = simple_closures
[]

[FluidProperties]
  [water]
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

[Executioner]
  type = Transient
  dt = 0.01
[]
