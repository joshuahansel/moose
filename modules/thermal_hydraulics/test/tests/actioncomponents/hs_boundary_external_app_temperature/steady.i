# A plate heat structure ActionComponent whose outer boundary temperature is set from the 'T_ext'
# variable (normally transferred from an external application; set here by an initial condition)
# and that is insulated elsewhere must reach that temperature everywhere at steady state.

T_ext = 450

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
    type = HeatStructurePlate
    orientation = '1 0 0'
    position = '0 0 0'
    length = 1
    n_elems = 5

    depth = 0.5
    names = 'blk'
    widths = '0.1'
    n_part_elems = '3'
    solid_properties = 'hs_mat'
    solid_properties_T_ref = '300'

    initial_T = 300
  []

  [hs_boundary]
    type = HSBoundaryExternalAppTemperature
    boundary = 'hs:outer'
    hs = hs
  []
[]

[ICs]
  [T_ext_ic]
    type = ConstantIC
    variable = T_ext
    value = ${T_ext}
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
    value2 = ${T_ext}
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
