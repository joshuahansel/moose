!include part_base.i
!include part_shell.i
!include part_cylinder.i

emissivity_shell = 0.8
emissivity_cylinder = 0.8
gap_conductivity = 1e-5

[Mesh]
  [combiner]
    type = CombinerGenerator
    inputs = 'rename_shell_boundary_id rename_cylinder_block'
  []
[]

[ThermalContact]
  [thermal_contact]
    type = GapHeatTransfer
    variable = temperature
    primary = shell_inner
    secondary = cylinder_outer
    emissivity_primary = ${emissivity_shell}
    emissivity_secondary = ${emissivity_cylinder}
    quadrature = true
    gap_conductivity = ${gap_conductivity}
    gap_geometry_type = cylinder
    cylinder_axis_point_1 = '0 0 0'
    cylinder_axis_point_2 = '0 0 1'
  []
[]
