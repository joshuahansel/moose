# A cylindrical heat structure ActionComponent with a specified temperature on its outer boundary and
# insulated elsewhere must reach that temperature everywhere at steady state.

T_specified = 350

[SolidProperties]
  [hs_mat]
    type = ThermalFunctionSolidProperties
    rho = 8000
    cp = 500
    k = 15
  []
[]

[ActionComponents]
  [hs]
    type = HeatStructureCylindrical
    orientation = '1 0 0'
    position = '0 0 0'
    length = 1
    n_elems = 5

    inner_radius = 0.1
    names = 'blk'
    widths = '0.1'
    n_part_elems = '3'
    solid_properties = 'hs_mat'
    solid_properties_T_ref = '300'

    initial_T = 300
  []

  [hs_boundary]
    type = HSBoundarySpecifiedTemperature
    boundary = 'hs:outer'
    hs = hs
    T = ${T_specified}
  []
[]

[Postprocessors]
  [T_avg]
    type = ElementAverageValue
    variable = T_solid
  []
  [T_relerr]
    type = RelativeDifferencePostprocessor
    value1 = T_avg
    value2 = ${T_specified}
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  nl_abs_tol = 1e-10
[]

[Outputs]
  [out]
    type = CSV
    show = 'T_relerr'
  []
[]
