T_hs = 300
T_ambient1 = 500
htc1 = 100
T_ambient2 = 400
htc2 = 300
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
heat_flux_avg = ${fparse 0.5 * (htc1 * (T_ambient1 - T_hs) + htc2 * (T_ambient2 - T_hs))}
heat_flux_integral = ${fparse heat_flux_avg * A}
scale = 0.8
power = ${fparse scale * heat_flux_integral}
E_change = ${fparse power * t}

[FunctorMaterials]
  [test_fm]
    type = ADGenericFunctorMaterial
    prop_names = 'T_ambient_prop htc_ambient_prop'
    prop_values = 'T_ambient_fn htc_ambient_fn'
  []
[]

[Functions]
  [T_ambient_fn]
    type = PiecewiseConstant
    axis = z
    x = '0 1'
    y = '${T_ambient1} ${T_ambient2}'
  []
  [htc_ambient_fn]
    type = PiecewiseConstant
    axis = z
    x = '0 1'
    y = '${htc1} ${htc2}'
  []
[]

[SolidProperties]
  [hs_mat]
    type = ThermalFunctionSolidProperties
    rho = ${density}
    cp = ${specific_heat_capacity}
    k = ${conductivity}
  []
[]

[Postprocessors]
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
    value1 = ambient_convection_integral
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
