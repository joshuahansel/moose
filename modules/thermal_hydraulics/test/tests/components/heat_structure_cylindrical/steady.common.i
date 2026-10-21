# Tests that cylindrical heat structure geometry can be used with a steady executioner.

[Functions]
  [power_profile_fn]
    type = ParsedFunction
    expression = '1.570796326794897 * sin(x / 3.6576 * pi)'
  []
[]

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

[Preconditioning]
  [pc]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Steady

  solve_type = 'NEWTON'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  nl_max_its = 30

  l_tol = 1e-4
  l_max_its = 300
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
