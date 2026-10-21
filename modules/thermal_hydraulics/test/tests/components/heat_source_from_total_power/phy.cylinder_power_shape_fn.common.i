[GlobalParams]
  scaling_factor_temperature = 1e0
[]

[Functions]
  [psf]
    type = ParsedFunction
    expression = 1
  []
[]

[SolidProperties]
  [fuel-mat]
    type = ThermalFunctionSolidProperties
    k = 16
    cp = 191.67
    rho = 1.4583e4
  []
  [gap-mat]
    type = ThermalFunctionSolidProperties
    k = 64
    cp = 1272
    rho = 865
  []
  [clad-mat]
    type = ThermalFunctionSolidProperties
    k = 26
    cp = 638
    rho = 7.646e3
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
  dt = 1e-3
  num_steps = 1
  abort_on_solve_fail = true

  solve_type = 'PJFNK'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-7
  nl_max_its = 40

  l_tol = 1e-5
  l_max_its = 50
[]

[Outputs]
  [out]
    type = Exodus
  []
[]
