# Tests HeatStructure2DCoupler when the heat structures are separated by some
# distance. The first heat structure has a larger coupling surface than the
# second heat structure. The component will be used to model a given energy
# transfer rate per unit temperature difference [W/K]. This test checks that:
#   a) heat transfer occurs in the correct direction
#   b) energy is conserved
#
# With a goal of transferring 5 W/K and a temperature difference of 200 K, and
# a transient time of 10 seconds, ~10 kJ should be transferred. Note that this
# estimate will not be exact since the temperature difference changes slightly
# over the transient.

initial_T1 = 500
initial_T2 = 300

R1 = 0.1
R2 = 0.05

P2 = ${fparse 2 * pi * R2}

power_per_K = 5.0
L_hs = 0.5
htc = ${fparse power_per_K / (L_hs * P2)}

[SolidProperties]
  [hs_mat]
    type = ThermalFunctionSolidProperties
    k = 15
    cp = 500
    rho = 8000
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
  scheme = 'bdf2'

  start_time = 0
  dt = 1.0
  num_steps = 10
  abort_on_solve_fail = true

  solve_type = 'NEWTON'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  nl_max_its = 30

  l_tol = 1e-4
  l_max_its = 300
[]

[Outputs]
  csv = true
  show = 'E_hs1_change E_hs2_change E_tot_change'
[]
