[motor]
  type = ShaftConnectedMotor
  inertia = inertia_fn
  torque = torque_fn
[]

[shaft]
  type = Shaft
  connected_components = 'motor'
  initial_speed = 0
[]

[hs]
  type = HeatStructureCylindrical
  position = '0 0 0'
  orientation = '1 0 0'
  length = 1
  n_elems = 1

  names = '0'
  n_part_elems = 1
  widths = '1'
  solid_properties = 'mat'
  solid_properties_T_ref = '300'

  initial_T = 300
[]
