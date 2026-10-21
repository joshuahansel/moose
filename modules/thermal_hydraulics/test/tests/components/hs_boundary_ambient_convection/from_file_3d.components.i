[hs]
  type = HeatStructureFromFile3D
  file = box.e
  position = '0 0 0'
  initial_T = ${T_hs}
[]

[ambient_convection]
  type = HSBoundaryAmbientConvection
  boundary = 'hs:left'
  hs = hs
  T_ambient = T_ambient_fn
  htc_ambient = htc_ambient_fn
  scale = ${scale}
[]
