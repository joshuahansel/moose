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
  [sp_ss316]
    type = ThermalSS316Properties
  []
[]

[Functions]
  [T0_fn]
    type = ParsedFunction
    expression = '290 + 20 * (y - 1)'
  []
[]

[Postprocessors]
  [T_avg]
    type = ElementAverageValue
    block = 'pipe'
    variable = T
    execute_on = 'INITIAL TIMESTEP_END'
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
  dt = 0.01
  num_steps = 1

  solve_type = NEWTON
  nl_rel_tol = 0
  nl_abs_tol = 1e-6
  nl_max_its = 15

  l_tol = 1e-3
  l_max_its = 10
[]

[Outputs]
  csv = true
[]
