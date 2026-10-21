  [hs]
    type = HeatStructureCylindrical
    orientation = '1 0 0'
    position = '0 0 0'
    length = 1
    n_elems = 10

    inner_radius = 0.1
    widths = '0.5'
    n_part_elems = '10'
    solid_properties = 'ss316'
    solid_properties_T_ref = '300'
    names = 'region'

    initial_T = 300
  []

  [ext_temperature]
    type = HSBoundaryExternalAppTemperature
    boundary = 'hs:outer'
    hs = hs
  []
