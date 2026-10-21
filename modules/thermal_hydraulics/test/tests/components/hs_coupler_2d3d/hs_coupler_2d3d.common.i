# Tests physics and energy conservation for HSCoupler2D3D.

R_pipe = 0.005
length_matrix = 0.5
length_extend = 0.6

n_elems_radial = 3
n_elems_axial_matrix = 10
n_elems_axial_extend = 12

[Materials]
  [matrix_mat]
    type = ADGenericConstantMaterial
    block = 'hs3d:0 hs2d:pipe'
    prop_names = 'density specific_heat thermal_conductivity'
    prop_values = '8000 500 15'
  []
[]

[Functions]
  [initial_T_matrix_fn]
    type = ParsedFunction
    expression = '300 + 100*z - 1000*x'
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
  dt = 0.1
  num_steps = 10

  solve_type = NEWTON
  abort_on_solve_fail = true
  nl_abs_tol = 1e-8
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  [exodus]
    type = Exodus
    output_dimension = 3 # classic components write 3D meshes
    append_object_name = false
  []
[]
