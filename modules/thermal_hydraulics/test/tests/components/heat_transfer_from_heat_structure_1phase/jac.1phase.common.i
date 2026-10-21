[GlobalParams]
  initial_p = 1.e5
  initial_vel = 2
  initial_T = 300

  scaling_factor_1phase = '1 1 1'
  scaling_factor_temperature = '1'

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

  start_time = 0
  dt = 1
  num_steps = 1
  abort_on_solve_fail = true

  solve_type = 'NEWTON'

  petsc_options_iname = '-snes_test_err'
  petsc_options_value = ' 1e-11'
[]
