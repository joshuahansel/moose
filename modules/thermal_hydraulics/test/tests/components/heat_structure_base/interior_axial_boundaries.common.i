# This input file is used to test that the interior axial boundaries of a
# heat structure are being created correctly.
#
# To test this, an arbitrary temperature distribution is imposed on the heat
# structure, and the average temperature on the interior axial boundaries is
# tested against expected values.
#
# The interior axial boundaries are located at x={20,40}, and radial boundaries
# are located at y={0,0.5,1,1.5}. The temperature is set to be T(x,y) = xy. The
# following table gives the resulting expected average temperature values on
# each face:
#   Boundary                     T_avg
#   -----------------------------------
#   hs:radial1:axial1:axial2     5
#   hs:radial1:axial2:axial3     10
#   hs:radial2:axial1:axial2     15
#   hs:radial2:axial2:axial3     30
#   hs:radial3:axial1:axial2     25
#   hs:radial3:axial2:axial3     50

[Functions]
  [initial_T_fn]
    type = ParsedFunction
    expression = 'x * y'
  []
[]

[SolidProperties]
  [hs_mat]
    type = ThermalFunctionSolidProperties
    k = 1
    cp = 1
    rho = 1
  []
[]

[Postprocessors]
  [T_avg_radial1_axial1_axial2]
    type = SideAverageValue
    variable = T_solid
    boundary = 'hs:radial1:axial1:axial2'
    execute_on = 'INITIAL'
  []
  [T_avg_radial1_axial2_axial3]
    type = SideAverageValue
    variable = T_solid
    boundary = 'hs:radial1:axial2:axial3'
    execute_on = 'INITIAL'
  []
  [T_avg_radial2_axial1_axial2]
    type = SideAverageValue
    variable = T_solid
    boundary = 'hs:radial2:axial1:axial2'
    execute_on = 'INITIAL'
  []
  [T_avg_radial2_axial2_axial3]
    type = SideAverageValue
    variable = T_solid
    boundary = 'hs:radial2:axial2:axial3'
    execute_on = 'INITIAL'
  []
  [T_avg_radial3_axial1_axial2]
    type = SideAverageValue
    variable = T_solid
    boundary = 'hs:radial3:axial1:axial2'
    execute_on = 'INITIAL'
  []
  [T_avg_radial3_axial2_axial3]
    type = SideAverageValue
    variable = T_solid
    boundary = 'hs:radial3:axial2:axial3'
    execute_on = 'INITIAL'
  []
[]

[Preconditioning]
  [pc]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  num_steps = 0
[]

[Outputs]
  csv = true
  execute_on = 'INITIAL'
[]
