# Sub input file.

L = 5.0
radius = 0.01

n_elems_axial = 10
n_elems_radial = 5

T_initial = 300.0
power = 1000.0
t = 10.0
E_change = ${fparse power * t}

rho = 8000.0
cp = 500.0
k = 15.0

[SolidProperties]
  [hs_mat]
    type = ThermalFunctionSolidProperties
    rho = ${rho}
    cp = ${cp}
    k = ${k}
  []
[]

[Postprocessors]
  [E_hs_change]
    type = ChangeOverTimePostprocessor
    postprocessor = E_hs
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [E_change_relerr]
    type = RelativeDifferencePostprocessor
    value1 = E_hs_change
    value2 = ${E_change}
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [integral_relerr]
    type = RelativeDifferencePostprocessor
    value1 = hs_boundary_integral
    value2 = ${power}
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

  scheme = bdf2
  dt = ${t}
  num_steps = 1
  abort_on_solve_fail = true

  solve_type = NEWTON
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-8
  nl_max_its = 10

  l_tol = 1e-3
  l_max_its = 10

  [Quadrature]
    type = GAUSS
    order = SECOND
  []
[]

[Outputs]
  csv = true
  show = 'E_change_relerr integral_relerr'
  execute_on = 'FINAL'
[]
