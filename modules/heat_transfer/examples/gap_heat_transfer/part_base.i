rho = 8000
cp = 500
k = 15

[Variables]
  [temperature]
  []
[]

[Kernels]
  [time_derivative]
    type = ADHeatConductionTimeDerivative
    variable = temperature
  []
[]

[SolidProperties]
  [sp]
    type = ThermalFunctionSolidProperties
    rho = ${rho}
    cp = ${cp}
    k = ${k}
  []
[]

[Materials]
  [thermal_mat]
    type = ADThermalSolidPropertiesMaterial
    temperature = temperature
    sp = sp
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
  end_time = 1000
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.0
    optimal_iterations = 5
    iteration_window = 0
    growth_factor = 1.2
    cutback_factor = 0.8
  []

  steady_state_detection = true
  steady_state_tolerance = 1e-8

  solve_type = NEWTON

  nl_abs_tol = 1e-12
  nl_rel_tol = 1e-10
  nl_max_its = 20
[]

[Outputs]
  exodus = true
[]
