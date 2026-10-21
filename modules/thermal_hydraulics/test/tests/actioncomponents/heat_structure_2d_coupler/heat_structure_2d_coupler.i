!include ../../components/heat_structure_2d_coupler/heat_structure_2d_coupler.common.i

[ActionComponents]
  !include ../../components/heat_structure_2d_coupler/heat_structure_2d_coupler.components.i
[]

[Postprocessors]
  !include heat_structure_2d_coupler.energy.i
[]
