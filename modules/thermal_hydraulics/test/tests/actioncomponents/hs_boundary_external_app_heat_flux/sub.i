!include ../../components/hs_boundary_external_app_heat_flux/sub.common.i

[ActionComponents]
  !include ../../components/hs_boundary_external_app_heat_flux/sub.components.i
[]

[Postprocessors]
  !include sub.energy.i
[]
