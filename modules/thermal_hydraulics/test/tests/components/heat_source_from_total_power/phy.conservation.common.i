# Tests energy conservation for HeatGeneration component when a power component is used

n_units = 5
power = 1e5
power_fraction = 0.3
t = 1

energy_change = ${fparse power_fraction * power * t}

[GlobalParams]
  scaling_factor_temperature = 1e-3
[]

[Functions]
  [power_shape]
    type = ConstantFunction
    value = 0.4
  []
[]

[SolidProperties]
  [main-material]
    type = ThermalFunctionSolidProperties
    k = 1e4
    cp = 500.0
    rho = 100.0
  []
[]

[Postprocessors]
  [E_tot]
    type = ADHeatStructureEnergyRZ
    block = 'heat_structure:rgn1 heat_structure:rgn2 heat_structure:rgn3'
    n_units = ${n_units}
    execute_on = 'initial timestep_end'
  []
  [E_tot_change]
    type = ChangeOverTimePostprocessor
    change_with_respect_to_initial = true
    postprocessor = E_tot
    execute_on = 'initial timestep_end'
  []

  [E_tot_change_rel_err]
    type = RelativeDifferencePostprocessor
    value1 = E_tot_change
    value2 = ${energy_change}
    execute_on = 'initial timestep_end'
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

  solve_type = 'NEWTON'
  line_search = 'basic'
  petsc_options_iname = '-pc_type'
  petsc_options_value = ' lu'

  nl_rel_tol = 0
  nl_abs_tol = 1e-6
  nl_max_its = 15

  l_tol = 1e-3
  l_max_its = 10

  start_time = 0.0
  dt = ${t}
  num_steps = 1

  abort_on_solve_fail = true

  [Quadrature]
    type = GAUSS
    order = SECOND
  []
[]

[Outputs]
  csv = true
  show = 'E_tot_change_rel_err'
  execute_on = 'final'
[]
