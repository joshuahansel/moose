[SolidProperties]
  [ss316]
    type = ThermalFunctionSolidProperties
    rho = 8.0272e3
    cp = 502.1
    k = 16.26
  []
[]


[Executioner]
  type = Transient
  scheme = bdf2
  dt = 0.1
  abort_on_solve_fail = true
  solve_type = NEWTON
  line_search = basic

  nl_rel_tol = 1e-7
[]

[Outputs]
  [exodus]
    type = Exodus
    output_dimension = 3 # classic components write 3D meshes
    append_object_name = false # keep the sub-app file name that 'exodus = true' gives
    show = 'T_ext'
  []
[]
