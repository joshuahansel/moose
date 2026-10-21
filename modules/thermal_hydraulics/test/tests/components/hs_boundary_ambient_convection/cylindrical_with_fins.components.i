[pipe_with_fins]
  type = HeatStructureCylindrical
  orientation = '0 0 1'
  position = '0 0 0'
  length = '${L_uncooled} ${L_cooled}'
  n_elems = '${n_elems_uncooled} ${n_elems_cooled}'
  axial_region_names = 'uncooled cooled'

  names = 'body'
  widths = '${diameter}'
  n_part_elems = '${n_elems_radial}'
  solid_properties = 'sp_ss316'
  solid_properties_T_ref = '${T_ref}'

  initial_T = ${T_initial}
[]
[pipe_without_fins]
  type = HeatStructureCylindrical
  orientation = '0 0 1'
  position = '0 0.02 0'
  length = '${L_uncooled} ${L_cooled}'
  n_elems = '${n_elems_uncooled} ${n_elems_cooled}'
  axial_region_names = 'uncooled cooled'

  names = 'body'
  widths = '${diameter}'
  n_part_elems = '${n_elems_radial}'
  solid_properties = 'sp_ss316'
  solid_properties_T_ref = '${T_ref}'

  initial_T = ${T_initial}
[]

[pipe_with_fins_convection]
  type = HSBoundaryAmbientConvection
  boundary = 'pipe_with_fins:cooled:outer'
  hs = pipe_with_fins
  T_ambient = ${T_ambient}
  htc_ambient = ${htc_ambient}
  scale = fin_enhancement
[]
[pipe_without_fins_convection]
  type = HSBoundaryAmbientConvection
  boundary = 'pipe_without_fins:cooled:outer'
  hs = pipe_without_fins
  T_ambient = ${T_ambient}
  htc_ambient = ${htc_ambient}
[]
