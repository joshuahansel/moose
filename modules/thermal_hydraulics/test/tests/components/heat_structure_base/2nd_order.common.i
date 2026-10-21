# This tests ensures that 2nd-order meshes can be used; it checks for the
# "Solve Converged" string at the end of a time step.

[GlobalParams]
  2nd_order_mesh = true
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

  start_time = 0
  dt = 0.1
  num_steps = 1
  abort_on_solve_fail = true

  solve_type = 'PJFNK'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  nl_max_its = 30

  l_tol = 1e-4
  l_max_its = 300
[]
