# This input file tests checks that T_ext and htc_ext are properly
# transferred from the master app.

T_hs = 300

L = 1
thickness = 0.5
depth = 0.6

density = 8.0272e3
specific_heat_capacity = 502.1
conductivity = 20

scale = 0.8


[SolidProperties]
  [hs_mat]
    type = ThermalFunctionSolidProperties
    rho = ${density}
    cp = ${specific_heat_capacity}
    k = ${conductivity}
  []
[]


[Executioner]
  type = Transient
  abort_on_solve_fail = true
[]

[Outputs]
  [exodus]
    type = Exodus
    output_dimension = 3 # classic components write 3D meshes
    append_object_name = false # keep the sub-app file name that 'exodus = true' gives
  []
[]
