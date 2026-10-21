# This input file is used to test that HSCoupler2D2DRadiation conserves
# energy for a problem where 3 cylindrical heat structures (surfaces 1, 2, and 3)
# are enclosed by an annular heat structure (surface 4). Note that the mesh
# positions used in this input file do not reflect the real positions for this
# configuration, for convenience of viewing solutions.

emissivity1 = 0.8
emissivity2 = 0.5
emissivity3 = 0.2
emissivity4 = 0.9

orientation = '0 0 1'

length = 0.5
n_axial_elems = 10

radius_123 = 0.1
inner_radius_4 = 0.2
outer_radius_4 = 0.25
thickness_4 = ${fparse outer_radius_4 - inner_radius_4}
n_radial_elems_123 = 10
n_radial_elems_4 = 5

initial_T1 = 1200
initial_T2 = 1000
initial_T3 = 800
initial_T4 = 300
T_ref = 300

y_shift = 0.5
position1 = '0 0 0'
position2 = '0 ${y_shift} 0'
position3 = '0 ${fparse 2*y_shift} 0'
position4 = '0 ${fparse 3*y_shift} 0'

view_factor_12 = 0.15 # guessed some number < 1/6
view_factor_14 = ${fparse 1.0 - 2 * view_factor_12}
view_factor_41 = ${fparse radius_123 / inner_radius_4 * view_factor_14}
view_factor_44 = ${fparse 1.0 - 3 * view_factor_41}

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
  file_base = 'energy_conservation'
  [csv]
    type = CSV
    show = 'E_tot_err'
    execute_on = 'FINAL'
  []
[]
