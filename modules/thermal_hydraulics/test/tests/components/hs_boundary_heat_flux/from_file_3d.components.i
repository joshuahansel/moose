[hs]
  type = HeatStructureFromFile3D
  file = box.e
  position = '0 0 0'
  initial_T = ${T_hs}
[]

[heat_flux_boundary]
  type = HSBoundaryHeatFlux
  boundary = 'hs:left'
  hs = hs
  q = q_fn
  scale = ${scale}
[]
