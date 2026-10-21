  [hs]
    type = HeatStructurePlate
    orientation = '0 0 1'
    position = '0 0 0'
    length = ${L}
    n_elems = 10

    depth = ${depth}
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
    T_ambient = ${T_ambient}
    emissivity = ${emissivity}
    view_factor = ${view_factor}
    scale = ${scale}
  []
