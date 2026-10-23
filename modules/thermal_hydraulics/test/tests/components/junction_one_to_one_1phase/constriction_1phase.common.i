# This test is used to test the JunctionOneToOne1Phase1Phase component with unequal areas
# at the junction. The downstream flow channel has an area half that of the
# upstream pipe, so there should be a pressure increase just upstream of the
# junction due to the partial wall. The velocity should increase through the
# junction (approximately by a factor of 2, but there are compressibility effects).

[GlobalParams]
  gravity_vector = '0 0 0'

  fp = fp
  closures = simple_closures
  f = 0

  initial_T = 300
  initial_p = 1e5
  initial_vel = 1
[]

[FluidProperties]
  [fp]
    type = IdealGasFluidProperties
    gamma = 1.4
    molar_mass = 11.64024372
  []
[]

[Closures]
  [simple_closures]
    type = Closures1PhaseSimple
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient

  scheme = bdf2
  dt = 0.01
  num_steps = 10
  abort_on_solve_fail = true

  solve_type = NEWTON
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-8
  nl_max_its = 60

  l_tol = 1e-4
[]

[Outputs]
  exodus = true
  show = 'p T vel'
  execute_on = 'initial timestep_end'
  velocity_as_vector = false
[]
