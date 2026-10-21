# Tests conservation of energy for convection on the boundary of a cylindrical heat structure
# ActionComponent with an ambient temperature and heat transfer coefficient from an external
# application ('T_ext'/'htc_ext', set here by initial conditions instead of transfers).

T_hs = 300
T_ext = 500
htc_ext = 100
t = 0.001

L = 2
D_i = 0.2
thickness = 0.5

# SS 316
density = 8.0272e3
specific_heat_capacity = 502.1
conductivity = 16.26

R_i = ${fparse 0.5 * D_i}
D_o = ${fparse D_i + 2 * thickness}
A = ${fparse pi * D_o * L}
scale = 0.8
power = ${fparse scale * htc_ext * (T_ext - T_hs) * A}
E_change = ${fparse power * t}

[SolidProperties]
  [hs_mat]
    type = ThermalFunctionSolidProperties
    rho = ${density}
    cp = ${specific_heat_capacity}
    k = ${conductivity}
  []
[]

[ActionComponents]
  [hs]
    type = HeatStructureCylindrical
    orientation = '0 0 1'
    position = '0 0 0'
    length = ${L}
    n_elems = 10

    inner_radius = ${R_i}
    widths = '${thickness}'
    n_part_elems = '10'
    solid_properties = 'hs_mat'
    solid_properties_T_ref = '300'
    names = 'region'

    initial_T = ${T_hs}
  []

  [hs_boundary]
    type = HSBoundaryExternalAppConvection
    boundary = 'hs:outer'
    hs = hs
    scale = ${scale}
  []
[]

[ICs]
  [T_ext_ic]
    type = ConstantIC
    variable = T_ext
    value = ${T_ext}
  []
  [htc_ext_ic]
    type = ConstantIC
    variable = htc_ext
    value = ${htc_ext}
  []
[]

[Postprocessors]
  [E_hs]
    type = ADHeatStructureEnergy
    block = 'hs:region'
    plate_depth = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
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

  [heat_rate_pp_relerr]
    type = RelativeDifferencePostprocessor
    value1 = hs_boundary_integral
    value2 = ${power}
    execute_on = 'INITIAL'
  []
[]

[Executioner]
  type = Transient

  [TimeIntegrator]
    type = ActuallyExplicitEuler
    solve_type = lumped
  []
  dt = ${t}
  num_steps = 1
  abort_on_solve_fail = true
[]

[Outputs]
  [out]
    type = CSV
    show = 'E_change_relerr heat_rate_pp_relerr'
    execute_on = 'FINAL'
  []
[]
