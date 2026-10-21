T_hs = 1200
T_ambient = 1500
emissivity = 0.3
view_factor = 0.6
t = 5.0

L = 2
thickness = 0.5
depth = 0.6

# SS 316
density = 8.0272e3
specific_heat_capacity = 502.1
conductivity = 16.26

stefan_boltzmann = 5.670367e-8
A = ${fparse L * depth}
heat_flux = ${fparse stefan_boltzmann * emissivity * view_factor * (T_ambient^4 - T_hs^4)}
scale = 0.8
E_change = ${fparse scale * heat_flux * A * t}

[SolidProperties]
  [hs_mat]
    type = ThermalFunctionSolidProperties
    rho = ${density}
    cp = ${specific_heat_capacity}
    k = ${conductivity}
  []
[]

[Postprocessors]
  [E_hs]
    type = ADHeatStructureEnergy
    block = 'hs:region'
    plate_depth = ${depth}
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
[]

[Executioner]
  type = Transient

  [TimeIntegrator]
    type = ActuallyExplicitEuler
  []
  dt = ${t}
  num_steps = 1
  abort_on_solve_fail = true

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  [out]
    type = CSV
    show = 'E_change_relerr'
    execute_on = 'FINAL'
  []
[]
