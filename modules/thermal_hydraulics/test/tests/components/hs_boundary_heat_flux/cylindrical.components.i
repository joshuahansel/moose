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

  [heat_flux_boundary]
    type = HSBoundaryHeatFlux
    boundary = 'hs:outer'
    hs = hs
    q = q_fn
    scale = ${scale}
  []
