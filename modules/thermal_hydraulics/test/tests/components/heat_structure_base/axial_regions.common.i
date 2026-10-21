# This input file is used to test the ability to specify axial regions of a heat
# structure. A heat structure is split into 3 axial regions, and a boundary
# condition is applied to only the bottom 2 regions. A single time step is
# taken, and the output should show heat transfer only at the bottom 2
# boundaries.

[SolidProperties]
  [hs_mat]
    type = ThermalFunctionSolidProperties
    k = 5
    cp = 300
    rho = 100
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

  start_time = 0
  dt = 1
  num_steps = 1
  abort_on_solve_fail = true

  solve_type = PJFNK
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  nl_max_its = 15

  l_tol = 1e-3
  l_max_its = 10

  [Quadrature]
    type = GAUSS
    order = SECOND
  []
[]


[Outputs]
  exodus = true
[]
