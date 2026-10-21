!include ../../components/heat_structure_2d_radiation_coupler_rz/heat_structure_2d_radiation_coupler_rz.common.i

[ActionComponents]
  !include ../../components/heat_structure_2d_radiation_coupler_rz/heat_structure_2d_radiation_coupler_rz.components.i
[]

[Postprocessors]
  !include heat_structure_2d_radiation_coupler_rz.energy.i
[]
