!include part_base.i
!include part_hp_3d.i

gap_htc = 3000

[Variables]
  [temperature]
    initial_condition = ${T_initial}
  []
[]

[AuxVariables]
  [T_mono]
    initial_condition = ${T_initial}
  []
[]

[FunctorMaterials]
  [gap_heat_flux_fmat]
    type = ADConvectionHeatFluxFunctorMaterial
    heat_flux_name = gap_heat_flux
    htc = ${gap_htc}
    T_fluid = temperature
    T_solid = T_mono
  []
[]

[Kernels]
  [time_derivative]
    type = ADHeatConductionTimeDerivative
    variable = temperature
  []
  [heat_conduction]
    type = ADHeatConduction
    variable = temperature
  []
[]

[BCs]
  [gap_heat_flux_bc]
    type = FunctorNeumannBC
    variable = temperature
    boundary = hp_coupled_boundary
    functor = gap_heat_flux
  []
[]

[Postprocessors]
  [gap_heat_rate_main]
    type = Receiver
  []
  [gap_heat_rate]
    type = ADSideIntegralFunctorPostprocessor
    boundary = hp_coupled_boundary
    functor = gap_heat_flux
    functor_argument = qp
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
