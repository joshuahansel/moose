[GlobalParams]
  scaling_factor_1phase = '1 1 1e-3'
[]

[SolidProperties]
  [mat]
    type = ThermalFunctionSolidProperties
    rho = 1000
    cp = 100
    k = 30
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
    type = HeatStructureEnergy3D
    block = blk:0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [energy_fch]
    type = ElementIntegralVariablePostprocessor
    block = fch
    variable = rhoEA
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [total_energy]
    type = SumPostprocessor
    values = 'energy_fch energy_hs'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [energy_change]
    type = ChangeOverTimePostprocessor
    change_with_respect_to_initial = true
    postprocessor = total_energy
   compute_relative_change = false
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
  dt = 1
  solve_type = PJFNK
  line_search = basic
  num_steps = 1000
  steady_state_detection = true
  steady_state_tolerance = 1e-08
  nl_abs_tol = 1e-8
[]
