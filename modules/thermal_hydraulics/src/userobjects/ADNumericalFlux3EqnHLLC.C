//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADNumericalFlux3EqnHLLC.h"
#include "THMIndicesVACE.h"
#include "Numerics.h"

registerMooseObject("ThermalHydraulicsApp", ADNumericalFlux3EqnHLLC);

InputParameters
ADNumericalFlux3EqnHLLC::validParams()
{
  InputParameters params = ADNumericalFlux3EqnBase::validParams();
  params += NaNInterface::validParams();

  MooseEnum wave_speed_formulation("einfeldt davis", "einfeldt");
  params.addParam<MooseEnum>(
      "wave_speed_formulation", wave_speed_formulation, "Method for computing wave speeds");

  params.addRequiredParam<UserObjectName>("fluid_properties",
                                          "Name for fluid properties user object");

  params.addClassDescription("Computes internal side flux for the 1-D, 1-phase, variable-area "
                             "Euler equations using the HLLC approximate Riemann solver.");

  return params;
}

ADNumericalFlux3EqnHLLC::ADNumericalFlux3EqnHLLC(const InputParameters & parameters)
  : ADNumericalFlux3EqnBase(parameters),
    NaNInterface(this),
    _fp(getUserObject<SinglePhaseFluidProperties>("fluid_properties")),
    _wave_speed_formulation(
        getParam<MooseEnum>("wave_speed_formulation").getEnum<WaveSpeedFormulation>())
{
}

void
ADNumericalFlux3EqnHLLC::calcFlux(const std::vector<ADReal> & UL_3d,
                                  const std::vector<ADReal> & UR_3d,
                                  const RealVectorValue & nLR,
                                  const RealVectorValue & t1,
                                  const RealVectorValue & t2,
                                  std::vector<ADReal> & FL,
                                  std::vector<ADReal> & FR) const
{
  // compute the primitive variables

  ADReal rhoL, eL, pL, cL, unL, ut1L, ut2L, AL;
  computeREPCUA3D(UL_3d, nLR, t1, t2, rhoL, eL, pL, cL, unL, ut1L, ut2L, AL);

  ADReal rhoR, eR, pR, cR, unR, ut1R, ut2R, AR;
  computeREPCUA3D(UR_3d, nLR, t1, t2, rhoR, eR, pR, cR, unR, ut1R, ut2R, AR);

  const auto n_passives = UL_3d.size() - THMVACE3D::N_FLUX_INPUTS;
  std::vector<ADReal> passivesL(n_passives, 0.0), passivesR(n_passives, 0.0);
  for (const auto i : make_range(n_passives))
  {
    passivesL[i] = UL_3d[THMVACE3D::N_FLUX_INPUTS + i] / AL;
    passivesR[i] = UR_3d[THMVACE3D::N_FLUX_INPUTS + i] / AR;
  }

  // compute wave speeds
  ADReal sL, sR, sm;
  computeWaveSpeeds(
      rhoL, unL, ut1L, ut2L, eL, pL, cL, rhoR, unR, ut1R, ut2R, eR, pR, cR, sL, sR, sm);

  // compute flow area
  const auto UL_1d = convert3Dto1D(UL_3d);
  const auto UR_1d = convert3Dto1D(UR_3d);
  const ADReal A_flow = computeFlowArea(UL_1d, UR_1d);

  // compute the fluxes
  FL.resize(THMVACE3D::N_FLUX_OUTPUTS + n_passives);
  if (sL > 0.0)
  {
    const ADReal rhoEL = UL_3d[THMVACE3D::RHOEA] / AL;

    FL[THMVACE3D::MASS] = unL * rhoL * A_flow;
    FL[THMVACE3D::MOM_NORM] = (unL * rhoL * unL + pL) * A_flow;
    FL[THMVACE3D::MOM_TAN1] = rhoL * unL * ut1L * A_flow;
    FL[THMVACE3D::MOM_TAN2] = rhoL * unL * ut2L * A_flow;
    FL[THMVACE3D::ENERGY] = unL * (rhoEL + pL) * A_flow;
    for (const auto i : make_range(n_passives))
      FL[THMVACE3D::N_FLUX_OUTPUTS + i] = passivesL[i] * unL * A_flow;

    _last_region_index = 0;
  }
  else if (sL <= 0.0 && sm > 0.0)
  {
    const ADReal ps = starPressure(rhoL, unL, pL, sL, sm);
    const ADReal omegL = 1.0 / (sL - sm);
    const ADReal rhoLs = omegL * (sL - unL) * rhoL;
    const ADReal rhounLs = omegL * ((sL - unL) * rhoL * unL + ps - pL);
    const ADReal rhoEL = UL_3d[THMVACE3D::RHOEA] / AL;
    const ADReal rhoELs = omegL * ((sL - unL) * rhoEL - pL * unL + ps * sm);

    FL[THMVACE3D::MASS] = sm * rhoLs * A_flow;
    FL[THMVACE3D::MOM_NORM] = (sm * rhounLs + ps) * A_flow;
    FL[THMVACE3D::MOM_TAN1] = rhounLs * ut1L * A_flow;
    FL[THMVACE3D::MOM_TAN2] = rhounLs * ut2L * A_flow;
    FL[THMVACE3D::ENERGY] = sm * (rhoELs + ps) * A_flow;
    for (const auto i : make_range(n_passives))
    {
      const auto passiveLs = omegL * (sL - unL) * passivesL[i];
      FL[THMVACE3D::N_FLUX_OUTPUTS + i] = passiveLs * sm * A_flow;
    }

    _last_region_index = 1;
  }
  else if (sm <= 0.0 && sR >= 0.0)
  {
    const ADReal ps = starPressure(rhoL, unL, pL, sL, sm);
    const ADReal omegR = 1.0 / (sR - sm);
    const ADReal rhoRs = omegR * (sR - unR) * rhoR;
    const ADReal rhounRs = omegR * ((sR - unR) * rhoR * unR + ps - pR);
    const ADReal rhoER = UR_3d[THMVACE3D::RHOEA] / AR;
    const ADReal rhoERs = omegR * ((sR - unR) * rhoER - pR * unR + ps * sm);

    FL[THMVACE3D::MASS] = sm * rhoRs * A_flow;
    FL[THMVACE3D::MOM_NORM] = (sm * rhounRs + ps) * A_flow;
    FL[THMVACE3D::MOM_TAN1] = rhounRs * ut1R * A_flow;
    FL[THMVACE3D::MOM_TAN2] = rhounRs * ut2R * A_flow;
    FL[THMVACE3D::ENERGY] = sm * (rhoERs + ps) * A_flow;
    for (const auto i : make_range(n_passives))
    {
      const auto passiveRs = omegR * (sR - unR) * passivesR[i];
      FL[THMVACE3D::N_FLUX_OUTPUTS + i] = passiveRs * sm * A_flow;
    }

    _last_region_index = 2;
  }
  else if (sR < 0.0)
  {
    const ADReal rhoER = UR_3d[THMVACE3D::RHOEA] / AR;

    FL[THMVACE3D::MASS] = unR * rhoR * A_flow;
    FL[THMVACE3D::MOM_NORM] = (unR * rhoR * unR + pR) * A_flow;
    FL[THMVACE3D::MOM_TAN1] = rhoR * unR * ut1R * A_flow;
    FL[THMVACE3D::MOM_TAN2] = rhoR * unR * ut2R * A_flow;
    FL[THMVACE3D::ENERGY] = unR * (rhoER + pR) * A_flow;
    for (const auto i : make_range(n_passives))
      FL[THMVACE3D::N_FLUX_OUTPUTS + i] = passivesR[i] * unR * A_flow;

    _last_region_index = 3;
  }
  else
    std::fill(FL.begin(), FL.end(), getNaN());

  FR = FL;

  const ADReal A_wall_L = AL - A_flow;
  FL[THMVACE3D::MOM_NORM] += pL * A_wall_L;

  const ADReal A_wall_R = AR - A_flow;
  FR[THMVACE3D::MOM_NORM] += pR * A_wall_R;

  // std::cout<<"AL="<<AL.value()<<", AR="<<AR.value()<<", FL="<<FL[THMVACE3D::MOM_NORM].value()<<",
  // FR="<<FR[THMVACE3D::MOM_NORM].value()<<std::endl;
}

ADReal
ADNumericalFlux3EqnHLLC::computeFlowArea(const std::vector<ADReal> & UL,
                                         const std::vector<ADReal> & UR) const
{
  return std::min(UL[THMVACE1D::AREA], UR[THMVACE1D::AREA]);
}

ADReal
ADNumericalFlux3EqnHLLC::computeRiemannPressure1D(const std::vector<ADReal> & UL_1d,
                                                  const std::vector<ADReal> & UR_1d,
                                                  const Real nLR_dot_d) const
{
  ADReal rhoL, eL, pL, cL, unL, AL;
  computeREPCUA1D(UL_1d, nLR_dot_d, rhoL, eL, pL, cL, unL, AL);
  const ADReal ut1L = 0;
  const ADReal ut2L = 0;

  ADReal rhoR, eR, pR, cR, unR, AR;
  computeREPCUA1D(UR_1d, nLR_dot_d, rhoR, eR, pR, cR, unR, AR);
  const ADReal ut1R = 0;
  const ADReal ut2R = 0;

  ADReal sL, sR, sm;
  computeWaveSpeeds(
      rhoL, unL, ut1L, ut2L, eL, pL, cL, rhoR, unR, ut1R, ut2R, eR, pR, cR, sL, sR, sm);

  if (sL > 0.0)
    return pL;
  else if (sR < 0.0)
    return pR;
  else
    return starPressure(rhoL, unL, pL, sL, sm);
}

void
ADNumericalFlux3EqnHLLC::computeWaveSpeeds(const ADReal & rhoL,
                                           const ADReal & unL,
                                           const ADReal & ut1L,
                                           const ADReal & ut2L,
                                           const ADReal & eL,
                                           const ADReal & pL,
                                           const ADReal & cL,
                                           const ADReal & rhoR,
                                           const ADReal & unR,
                                           const ADReal & ut1R,
                                           const ADReal & ut2R,
                                           const ADReal & eR,
                                           const ADReal & pR,
                                           const ADReal & cR,
                                           ADReal & sL,
                                           ADReal & sR,
                                           ADReal & sm) const
{
  // compute left and right wave speeds
  if (_wave_speed_formulation == WaveSpeedFormulation::EINFELDT)
    computeWaveSpeedsEinfeldt(
        rhoL, unL, ut1L, ut2L, eL, pL, cL, rhoR, unR, ut1R, ut2R, eR, pR, cR, sL, sR);
  else if (_wave_speed_formulation == WaveSpeedFormulation::DAVIS)
    computeWaveSpeedsDavis(unL, cL, unR, cR, sL, sR);
  else
  {
    mooseAssert(false, "Invalid 'wave_speed_formulation'.");
  }

  // compute middle wave speed
  sm = (rhoR * unR * (sR - unR) - rhoL * unL * (sL - unL) + pL - pR) /
       (rhoR * (sR - unR) - rhoL * (sL - unL));
}

void
ADNumericalFlux3EqnHLLC::computeWaveSpeedsEinfeldt(const ADReal & rhoL,
                                                   const ADReal & unL,
                                                   const ADReal & ut1L,
                                                   const ADReal & ut2L,
                                                   const ADReal & eL,
                                                   const ADReal & pL,
                                                   const ADReal & cL,
                                                   const ADReal & rhoR,
                                                   const ADReal & unR,
                                                   const ADReal & ut1R,
                                                   const ADReal & ut2R,
                                                   const ADReal & eR,
                                                   const ADReal & pR,
                                                   const ADReal & cR,
                                                   ADReal & sL,
                                                   ADReal & sR) const
{
  ADReal un_roe, c_roe;
  computeRoeSpeeds(rhoL, unL, ut1L, ut2L, eL, pL, rhoR, unR, ut1R, ut2R, eR, pR, un_roe, c_roe);

  sL = min(unL - cL, un_roe - c_roe);
  sR = max(unR + cR, un_roe + c_roe);
}

void
ADNumericalFlux3EqnHLLC::computeWaveSpeedsDavis(const ADReal & unL,
                                                const ADReal & cL,
                                                const ADReal & unR,
                                                const ADReal & cR,
                                                ADReal & sL,
                                                ADReal & sR) const
{
  sL = min(unL - cL, unR - cR);
  sR = max(unL + cL, unR + cR);
}

std::vector<ADReal>
ADNumericalFlux3EqnHLLC::convert3Dto1D(const std::vector<ADReal> & U_3d) const
{
  std::vector<ADReal> U_1d(THMVACE1D::N_FLUX_INPUTS);
  U_1d[THMVACE1D::RHOA] = U_3d[THMVACE3D::RHOA];
  U_1d[THMVACE1D::RHOUA] = U_3d[THMVACE3D::RHOUA];
  U_1d[THMVACE1D::RHOEA] = U_3d[THMVACE3D::RHOEA];
  U_1d[THMVACE1D::AREA] = U_3d[THMVACE3D::AREA];
  return U_1d;
}

void
ADNumericalFlux3EqnHLLC::computeREPCUA3D(const std::vector<ADReal> & U_3d,
                                         const RealVectorValue & nLR,
                                         const RealVectorValue & t1,
                                         const RealVectorValue & t2,
                                         ADReal & rho,
                                         ADReal & e,
                                         ADReal & p,
                                         ADReal & c,
                                         ADReal & un,
                                         ADReal & ut1,
                                         ADReal & ut2,
                                         ADReal & A) const
{
  const ADReal & rhoA = U_3d[THMVACE3D::RHOA];
  const ADReal & rhouA = U_3d[THMVACE3D::RHOUA];
  const ADReal & rhovA = U_3d[THMVACE3D::RHOVA];
  const ADReal & rhowA = U_3d[THMVACE3D::RHOWA];
  const ADReal & rhoEA = U_3d[THMVACE3D::RHOEA];
  A = U_3d[THMVACE3D::AREA];

  rho = rhoA / A;
  const ADRealVectorValue uvec(rhouA / rhoA, rhovA / rhoA, rhowA / rhoA);
  un = uvec * nLR;
  ut1 = uvec * t1;
  ut2 = uvec * t2;
  const ADReal v = 1.0 / rho;
  e = rhoEA / rhoA - 0.5 * uvec * uvec;
  p = _fp.p_from_v_e(v, e);
  c = _fp.c_from_v_e(v, e);
}

void
ADNumericalFlux3EqnHLLC::computeREPCUA1D(const std::vector<ADReal> & U_1d,
                                         const Real nLR_dot_d,
                                         ADReal & rho,
                                         ADReal & e,
                                         ADReal & p,
                                         ADReal & c,
                                         ADReal & un,
                                         ADReal & A) const
{
  const ADReal rhoA = U_1d[THMVACE1D::RHOA];
  const ADReal rhouA = U_1d[THMVACE1D::RHOUA];
  const ADReal rhoEA = U_1d[THMVACE1D::RHOEA];
  A = U_1d[THMVACE1D::AREA];

  rho = rhoA / A;
  un = rhouA / rhoA * nLR_dot_d;
  const ADReal E = rhoEA / rhoA;
  e = E - 0.5 * un * un;
  const ADReal v = 1.0 / rho;
  p = _fp.p_from_v_e(v, e);
  c = _fp.c_from_v_e(v, e);
}

void
ADNumericalFlux3EqnHLLC::computeRoeSpeeds(const ADReal & rhoL,
                                          const ADReal & unL,
                                          const ADReal & ut1L,
                                          const ADReal & ut2L,
                                          const ADReal & eL,
                                          const ADReal & pL,
                                          const ADReal & rhoR,
                                          const ADReal & unR,
                                          const ADReal & ut1R,
                                          const ADReal & ut2R,
                                          const ADReal & eR,
                                          const ADReal & pR,
                                          ADReal & un_roe,
                                          ADReal & c_roe) const
{
  const ADReal sqrt_rhoL = sqrt(rhoL);
  const ADReal sqrt_rhoR = sqrt(rhoR);
  un_roe = (sqrt_rhoL * unL + sqrt_rhoR * unR) / (sqrt_rhoL + sqrt_rhoR);
  const ADReal ut1_roe = (sqrt_rhoL * ut1L + sqrt_rhoR * ut1R) / (sqrt_rhoL + sqrt_rhoR);
  const ADReal ut2_roe = (sqrt_rhoL * ut2L + sqrt_rhoR * ut2R) / (sqrt_rhoL + sqrt_rhoR);
  const ADRealVectorValue uvecL(unL, ut1L, ut2L);
  const ADRealVectorValue uvecR(unR, ut1R, ut2R);
  const ADReal EL = eL + 0.5 * uvecL * uvecL;
  const ADReal ER = eR + 0.5 * uvecR * uvecR;
  const ADReal HL = EL + pL / rhoL;
  const ADReal HR = ER + pR / rhoR;
  const ADReal H_roe = (sqrt_rhoL * HL + sqrt_rhoR * HR) / (sqrt_rhoL + sqrt_rhoR);
  const ADRealVectorValue uvec_roe(un_roe, ut1_roe, ut2_roe);
  const ADReal h_roe = H_roe - 0.5 * uvec_roe * uvec_roe;
  const ADReal rho_roe = sqrt(rhoL * rhoR);
  const ADReal v_roe = 1.0 / rho_roe;
  const ADReal e_roe = _fp.e_from_v_h(v_roe, h_roe);
  c_roe = _fp.c_from_v_e(v_roe, e_roe);
}

ADReal
ADNumericalFlux3EqnHLLC::starPressure(const ADReal & rhoL,
                                      const ADReal & unL,
                                      const ADReal & pL,
                                      const ADReal & sL,
                                      const ADReal & sm) const
{
  return rhoL * (sL - unL) * (sm - unL) + pL;
}
