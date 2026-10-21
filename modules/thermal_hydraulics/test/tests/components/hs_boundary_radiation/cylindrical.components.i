  [hs]
    type = HeatStructureCylindrical
    orientation = '0 0 1'
    position = '0 0 0'
    length = ${L}
    n_elems = 10

    inner_radius = ${R_i}
    widths = '${thickness}'
    n_part_elems = '10'
    solid_properties = 'hs_mat'
    solid_properties_T_ref = '300'
    names = 'region'

    initial_T = ${T_hs}
  []

  [hs_boundary]
    type = HSBoundaryRadiation
    boundary = 'hs:outer'
    hs = hs
    T_ambient = T_ambient_prop
    emissivity = emissivity_prop
    view_factor = view_factor_prop
    scale = scale_prop
  []
