//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "StiffenedGasMatchTwoPhaseFluidPropertiesTest.h"
#include "SinglePhaseFluidPropertiesTestUtils.h"
#include "SodiumTwoPhaseFluidProperties.h"
#include "StiffenedGasMatchTwoPhaseFluidProperties.h"
#include "TwoPhaseFluidProperties.h"
#include "SinglePhaseFluidProperties.h"

StiffenedGasMatchTwoPhaseFluidPropertiesTest::StiffenedGasMatchTwoPhaseFluidPropertiesTest()
  : MooseObjectUnitTest("FluidPropertiesApp"), _T(1000)
{
  buildObjects();
}

void
StiffenedGasMatchTwoPhaseFluidPropertiesTest::buildObjects()
{
  const std::string fp_2phase_ref_name = "fp_2phase_ref";
  const std::string fp_2phase_name = "fp_2phase";

  // 2-phase reference
  {
    const std::string class_name = "SodiumTwoPhaseFluidProperties";
    InputParameters params = _factory.getValidParams(class_name);
    _fe_problem->addUserObject(class_name, fp_2phase_ref_name, params);
    _fp_2phase_ref = &_fe_problem->getUserObject<TwoPhaseFluidProperties>(fp_2phase_ref_name);
  }
  _fp_2phase_ref->initialSetup();

  // 2-phase fit
  {
    const std::string class_name = "StiffenedGasMatchTwoPhaseFluidProperties";
    InputParameters params = _factory.getValidParams(class_name);
    params.set<UserObjectName>("fluid_properties") = {fp_2phase_ref_name};
    params.set<Real>("T") = _T;
    _fe_problem->addUserObject(class_name, fp_2phase_name, params);
    _fp_2phase = &_fe_problem->getUserObject<TwoPhaseFluidProperties>(fp_2phase_name);
  }
  _fp_2phase->initialSetup();
}

TEST_F(StiffenedGasMatchTwoPhaseFluidPropertiesTest, test)
{
  const Real p = _fp_2phase_ref->p_sat(_T);

  // liquid

  SinglePhaseFluidProperties & fp_liquid = _fe_problem->getUserObject<SinglePhaseFluidProperties>(_fp_2phase->getLiquidName());
  fp_liquid.initialSetup();

  SinglePhaseFluidProperties & fp_liquid_ref = _fe_problem->getUserObject<SinglePhaseFluidProperties>(_fp_2phase_ref->getLiquidName());
  fp_liquid_ref.initialSetup();

  REL_TEST(fp_liquid.molarMass(), fp_liquid_ref.molarMass(), REL_TOL_SAVED_VALUE);
  REL_TEST(fp_liquid.cp_from_p_T(p, _T), fp_liquid_ref.cp_from_p_T(p, _T), REL_TOL_SAVED_VALUE);
  REL_TEST(fp_liquid.cv_from_p_T(p, _T), fp_liquid_ref.cv_from_p_T(p, _T), REL_TOL_SAVED_VALUE);
  REL_TEST(fp_liquid.h_from_p_T(p, _T), fp_liquid_ref.h_from_p_T(p, _T), REL_TOL_SAVED_VALUE);
  REL_TEST(fp_liquid.rho_from_p_T(p, _T), fp_liquid_ref.rho_from_p_T(p, _T), REL_TOL_SAVED_VALUE);
  REL_TEST(fp_liquid.s_from_p_T(p, _T), fp_liquid_ref.s_from_p_T(p, _T), REL_TOL_SAVED_VALUE);
  REL_TEST(fp_liquid.mu_from_p_T(p, _T), fp_liquid_ref.mu_from_p_T(p, _T), REL_TOL_SAVED_VALUE);
  REL_TEST(fp_liquid.k_from_p_T(p, _T), fp_liquid_ref.k_from_p_T(p, _T), REL_TOL_SAVED_VALUE);

  // vapor

  SinglePhaseFluidProperties & fp_vapor = _fe_problem->getUserObject<SinglePhaseFluidProperties>(_fp_2phase->getVaporName());
  fp_vapor.initialSetup();

  SinglePhaseFluidProperties & fp_vapor_ref = _fe_problem->getUserObject<SinglePhaseFluidProperties>(_fp_2phase_ref->getVaporName());
  fp_vapor_ref.initialSetup();

  REL_TEST(fp_vapor.molarMass(), fp_vapor_ref.molarMass(), REL_TOL_SAVED_VALUE);
  REL_TEST(fp_vapor.cv_from_p_T(p, _T), fp_vapor_ref.cv_from_p_T(p, _T), REL_TOL_SAVED_VALUE);
  REL_TEST(fp_vapor.h_from_p_T(p, _T), fp_vapor_ref.h_from_p_T(p, _T), REL_TOL_SAVED_VALUE);
  REL_TEST(fp_vapor.mu_from_p_T(p, _T), fp_vapor_ref.mu_from_p_T(p, _T), REL_TOL_SAVED_VALUE);
  REL_TEST(fp_vapor.k_from_p_T(p, _T), fp_vapor_ref.k_from_p_T(p, _T), REL_TOL_SAVED_VALUE);

  // 2-phase

  REL_TEST(_fp_2phase->p_critical(), _fp_2phase_ref->p_critical(), REL_TOL_SAVED_VALUE);
  REL_TEST(_fp_2phase->T_triple(), _fp_2phase_ref->T_triple(), REL_TOL_SAVED_VALUE);
  REL_TEST(_fp_2phase->L_fusion(), _fp_2phase_ref->L_fusion(), REL_TOL_SAVED_VALUE);
  REL_TEST(_fp_2phase->T_sat(p), _fp_2phase_ref->T_sat(p), REL_TOL_SAVED_VALUE);
  REL_TEST(_fp_2phase->p_sat(_T), _fp_2phase_ref->p_sat(_T), REL_TOL_SAVED_VALUE);
  REL_TEST(_fp_2phase->dT_sat_dp(p), _fp_2phase_ref->dT_sat_dp(p), REL_TOL_SAVED_VALUE);
  REL_TEST(_fp_2phase->sigma_from_T(_T), _fp_2phase_ref->sigma_from_T(_T), REL_TOL_SAVED_VALUE);
  REL_TEST(_fp_2phase->dsigma_dT_from_T(_T), _fp_2phase_ref->dsigma_dT_from_T(_T), REL_TOL_SAVED_VALUE);
}
