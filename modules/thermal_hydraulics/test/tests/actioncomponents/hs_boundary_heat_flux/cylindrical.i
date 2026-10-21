!include ../../components/hs_boundary_heat_flux/cylindrical.common.i

[ActionComponents]
  !include ../../components/hs_boundary_heat_flux/cylindrical.components.i
[]

[Postprocessors]
  !include cylindrical.energy.i
[]
