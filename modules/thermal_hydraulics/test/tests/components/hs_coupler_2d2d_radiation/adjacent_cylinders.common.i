# This input file is used to test that HSCoupler2D2DRadiation can perform
# radiative heat transfer between multiple heat structures (surfaces 1 and 2)
# and the environment (surface 3).

emissivity1 = 0.8
emissivity2 = 0.5

orientation = '0 0 1'

length = 0.5
n_axial_elems = 10

radius = 0.1
n_radial_elems = 10

initial_T1 = 1200
initial_T2 = 1000
T3 = 300
T_ref = 300

y_shift = 0.5
position1 = '0 0 0'
position2 = '0 ${y_shift} 0'

view_factor_12 = ${fparse (pi - 2) / (2*pi)}
view_factor_13 = ${fparse 1.0 - view_factor_12}

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
  file_base = 'adjacent_cylinders'
  [exodus]
    type = Exodus
    output_dimension = 3 # classic components write 3D meshes
    append_object_name = false
  []
[]
