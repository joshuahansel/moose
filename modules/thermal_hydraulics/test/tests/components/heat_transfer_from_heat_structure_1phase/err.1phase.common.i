[GlobalParams]
  initial_p = 1e5
  initial_vel = 0
  initial_T = 300

  closures = simple_closures
[]

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

[SolidProperties]
  [fuel-mat]
    type = ThermalFunctionSolidProperties
    k = 2.5
    cp = 300.
    rho = 1.032e4
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
  dt = 1.e-5
  solve_type = 'NEWTON'
  num_steps = 1
  abort_on_solve_fail = true
[]
