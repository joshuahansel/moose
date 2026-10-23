[outlet]
  type = Outlet1Phase
  input = 'pipe:in'
  p = ${outlet_p}
[]
[pipe]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  length = 1.0
  n_elems = 5
  A = ${A_pipe}
  fp = fp
  closures = simple_closures
  f = 0
[]
[junction]
  type = VolumeJunction1Phase
  connections = 'pipe:out'
  position = '1.0 0 0'
  volume = ${V_junction}
[]
[junction_flux]
  type = VolumeJunctionCoupledFlux1Phase
  volume_junction = junction
  A_coupled = ${A_coupled}
  pp_suffix = test
  multi_app = child
  normal_from_junction = '1 0 0'
[]
