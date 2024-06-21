!include part_base.i
!include part_mono.i
!include part_hp_3d.i

emissivity_mono = 0.8
emissivity_clad = 0.8
gap_conductivity = 1e-5

[Mesh]
  [combiner]
    type = CombinerGenerator
    inputs = 'rename_mono_boundary_id rename_hp_block'
  []
[]

[Variables]
  [temperature]
    initial_condition = ${T_initial}
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

[ThermalContact]
  [thermal_contact]
    type = GapHeatTransfer
    variable = temperature
    primary = mono_inner
    secondary = hp_coupled_boundary
    emissivity_primary = ${emissivity_mono}
    emissivity_secondary = ${emissivity_clad}
    quadrature = true
    gap_conductivity = ${gap_conductivity}
    gap_geometry_type = cylinder
    cylinder_axis_point_1 = '0 0 0'
    cylinder_axis_point_2 = '0 0 1'
  []
[]
