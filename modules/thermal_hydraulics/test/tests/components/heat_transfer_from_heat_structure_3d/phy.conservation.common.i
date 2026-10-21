# Testing energy conservation with fluid at rest

P_hf = ${fparse 0.6 * sin (pi/24)}

[GlobalParams]
  gravity_vector = '0 0 0'
[]

[Materials]
  [mat]
    type = ADGenericConstantMaterial
    block = 'blk:0'
    prop_names = 'density specific_heat thermal_conductivity'
    prop_values = '1000 100 30'
  []
[]

[FluidProperties]
  [fp]
    type = StiffenedGasFluidProperties
    gamma = 2.35
    cv = 1816.0
    q = -1.167e6
    p_inf = 1.0e9
    q_prime = 0
  []
[]

[Closures]
  [simple_closures]
    type = Closures1PhaseSimple
  []
[]

[Functions]
  [T_init]
    type = ParsedFunction
    expression = '1000*y+300+30*z'
  []
[]

[Postprocessors]
  [energy_hs]
    type = ADHeatStructureEnergy3D
    block = blk:0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [energy_fch1]
    type = ElementIntegralVariablePostprocessor
    block = fch1
    variable = rhoEA
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [energy_fch2]
    type = ElementIntegralVariablePostprocessor
    block = fch2
    variable = rhoEA
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [total_energy]
    type = SumPostprocessor
    values = 'energy_fch1 energy_fch2 energy_hs'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [energy_change]
    type = ChangeOverTimePostprocessor
    change_with_respect_to_initial = true
    postprocessor = total_energy
    compute_relative_change = true
    execute_on = 'INITIAL TIMESTEP_END'
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
  line_search = basic
  abort_on_solve_fail = true
  nl_abs_tol = 1e-8
[]

[Outputs]
  file_base = 'phy.conservation'
  [csv]
    type = CSV
    show = 'energy_change'
    execute_on = 'FINAL'
  []
[]
