  [total_power]
    type = TotalPower
    power = 1234.
  []

  [ch1:solid]
    type = HeatStructureCylindrical
    position = '0 0 0'
    orientation = '1 0 0'
    length = 1
    n_elems = 1
    initial_T = 300
    names = '0'
    widths = '1'
    n_part_elems = '1'
    solid_properties = 'mat'
    solid_properties_T_ref = '300'
  []
