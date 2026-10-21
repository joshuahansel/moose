  [hs]
    type = HeatStructureCylindrical
    position = '0 0 0'
    orientation = '0 0 1'
    length = ${L}
    n_elems = ${n_elems_axial}

    names = 'body'
    widths = '${radius}'
    n_part_elems = '${n_elems_radial}'
    solid_properties = 'hs_mat'
    solid_properties_T_ref = '300'

    initial_T = ${T_initial}
  []
  [hs_boundary]
    type = HSBoundaryExternalAppHeatFlux
    hs = hs
    boundary = 'hs:outer'
    heat_flux_name = q_ext
    heat_flux_is_inward = true
    perimeter_ext = P_ext
  []
