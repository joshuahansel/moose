emissivity1 = 0.75
emissivity2 = 0.5

[SolidProperties]
  [hs_mat]
    type = ThermalFunctionSolidProperties
    k = 15
    cp = 500
    rho = 8000
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
  dt = 10
  num_steps = 10
  abort_on_solve_fail = true

  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10
  nl_max_its = 10

  l_tol = 1e-4
  l_max_its = 10
[]

[Outputs]
  file_base = 'heat_structure_2d_radiation_coupler_rz'
  [csv]
    type = CSV
    show = 'E_tot_change T1 T2'
  []
[]
