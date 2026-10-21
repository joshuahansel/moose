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

[Postprocessors]
[]

[Executioner]
  type = Transient
  scheme = 'bdf2'

  start_time = 0
  dt = 1000
  num_steps = 10
  abort_on_solve_fail = true

  solve_type = 'NEWTON'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  nl_max_its = 30

  l_tol = 1e-4
  l_max_its = 300
[]

[Outputs]
  file_base = 'cylindrical'
  [exodus]
    type = Exodus
    output_dimension = 3 # classic components write 3D meshes
    append_object_name = false
  []
[]
