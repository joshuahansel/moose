!include ../../components/hs_boundary_radiation/cylindrical.common.i

[ActionComponents]
  !include ../../components/hs_boundary_radiation/cylindrical.components.i
[]

[Postprocessors]
  !include cylindrical.energy.i
[]
