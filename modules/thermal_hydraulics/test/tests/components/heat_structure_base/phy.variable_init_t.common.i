# Tests that a function can be used to initialize temperature in a heat structure.

[GlobalParams]
[]

[Functions]
  [fn-initial_T]
    type = ParsedFunction
    expression = 'baseT + (dT * sin((pi * x) / length))'
    symbol_names = 'baseT   dT    length'
    symbol_values = '560.0  30.0  3.6576'
  []
[]

[SolidProperties]
  [fuel-mat]
    type = ThermalFunctionSolidProperties
    k = 3.65
    cp = 288.734
    rho = 1.0412e2
  []
  [gap-mat]
    type = ThermalFunctionSolidProperties
    k = 0.1
    cp = 1.0
    rho = 1.0
  []
  [clad-mat]
    type = ThermalFunctionSolidProperties
    k = 16.48672
    cp = 321.384
    rho = 6.6e1
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

  start_time = 0
  dt = 0.01
  num_steps = 10
  abort_on_solve_fail = true

  solve_type = 'PJFNK'
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-6
  nl_max_its = 8

  l_tol = 1e-4
  l_max_its = 10
[]


[Outputs]
  [out]
    type = Exodus
    output_dimension = 3 # classic components write 3D meshes
  []
  [console]
    type = Console
    execute_scalars_on = none
  []
[]
