//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "IdealGasMatchFluidPropertiesTest.h"
#include "SinglePhaseFluidPropertiesTestUtils.h"
#include "SodiumVaporFluidProperties.h"
#include "IdealGasMatchFluidProperties.h"

IdealGasMatchFluidPropertiesTest::IdealGasMatchFluidPropertiesTest()
: MooseObjectUnitTest("FluidPropertiesApp"),
  _p(1.9629094537e+04),
  _T(1000)
  {
    buildObjects();
  }

void
IdealGasMatchFluidPropertiesTest::buildObjects()
{
  const std::string fp_ref_name = "fp_ref";
  const std::string fp_igfit_name = "fp_igfit";

  // Reference
  {
    const std::string class_name = "SodiumVaporFluidProperties";
    InputParameters params = _factory.getValidParams(class_name);
    _fe_problem->addUserObject(class_name, fp_ref_name, params);
    _fp_ref = &_fe_problem->getUserObject<SinglePhaseFluidProperties>(fp_ref_name);
  }

  // IG fit
  {
    const std::string class_name = "IdealGasMatchFluidProperties";
    InputParameters params = _factory.getValidParams(class_name);
    params.set<UserObjectName>("fluid_properties") = {fp_ref_name};
    params.set<Real>("p") = _p;
    params.set<Real>("T") = _T;
    _fe_problem->addUserObject(class_name, fp_igfit_name, params);
    _fp_igfit = &_fe_problem->getUserObject<SinglePhaseFluidProperties>(fp_igfit_name);
  }

  _fp_igfit->initialSetup();
}


TEST_F(IdealGasMatchFluidPropertiesTest, test)
{
  REL_TEST(_fp_igfit->molarMass(), _fp_ref->molarMass(), REL_TOL_SAVED_VALUE);
  REL_TEST(_fp_igfit->cv_from_p_T(_p, _T), _fp_ref->cv_from_p_T(_p, _T), REL_TOL_SAVED_VALUE);
  REL_TEST(_fp_igfit->h_from_p_T(_p, _T), _fp_ref->h_from_p_T(_p, _T), REL_TOL_SAVED_VALUE);
  REL_TEST(_fp_igfit->mu_from_p_T(_p, _T), _fp_ref->mu_from_p_T(_p, _T), REL_TOL_SAVED_VALUE);
  REL_TEST(_fp_igfit->k_from_p_T(_p, _T), _fp_ref->k_from_p_T(_p, _T), REL_TOL_SAVED_VALUE);
}
