# Tests that error is generated when no initial temperature function is provided
# when not restarting.

[GlobalParams]
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
    k = 1.084498
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

  dt = 0.1
  dtmin = 1e-1

  solve_type = 'PJFNK'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  nl_max_its = 30

  l_tol = 1e-4
  l_max_its = 300

  start_time = 0.0
  end_time = 2.0
[]
