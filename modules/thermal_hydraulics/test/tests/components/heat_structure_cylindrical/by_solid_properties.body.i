[SolidProperties]
  [fuel_sp]
    type = ThermalFunctionSolidProperties
    rho = 1.0412e2
    cp = 288.734
    k = 3.65
  []
  [gap_sp]
    type = ThermalFunctionSolidProperties
    rho = 1.0
    cp = 1.0
    k = 1.084498
  []
  [clad_sp]
    type = ThermalFunctionSolidProperties
    rho = 6.6e1
    cp = 321.384
    k = 16.48672
  []
[]

[Convergence]
  [nl_conv]
    type = DefaultNonlinearConvergence
    nl_rel_tol = 1e-6
    nl_abs_tol = 1e-6
    nl_max_its = 30
  []
[]
