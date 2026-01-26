//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADNumericalFlux3EqnBase.h"
#include "SinglePhaseFluidProperties.h"
#include "NaNInterface.h"

/**
 * Computes internal side flux for the 1-D, 1-phase, variable-area Euler equations using
 * the HLLC approximate Riemann solver.
 *
 * The approach in the following reference for the 3-D Euler equations was
 * extended to the 1-D, variable-area Euler equations:
 *
 * Batten, P., Leschziner, M. A., & Goldberg, U. C. (1997).
 * Average-state Jacobians and implicit methods for compressible viscous and turbulent flows.
 * Journal of computational physics, 137(1), 38-78.
 */
class ADNumericalFlux3EqnHLLC : public ADNumericalFlux3EqnBase, public NaNInterface
{
public:
  ADNumericalFlux3EqnHLLC(const InputParameters & parameters);

  virtual void calcFlux(const std::vector<ADReal> & UL_3d,
                        const std::vector<ADReal> & UR_3d,
                        const RealVectorValue & nLR,
                        const RealVectorValue & t1,
                        const RealVectorValue & t2,
                        std::vector<ADReal> & FL,
                        std::vector<ADReal> & FR) const override;

  virtual unsigned int getNumberOfRegions() const override { return 4; }

  /**
   * Computes the Riemann pressure for 1D inputs
   *
   * @param[in] UL_1d  Left 1D solution vector
   * @param[in] UR_1d  Right 1D solution vector
   * @param[in] nLR_dot_d  Dot product of left-to-right direction with flow channel direction
   */
  ADReal computeRiemannPressure1D(const std::vector<ADReal> & UL_1d,
                                  const std::vector<ADReal> & UR_1d,
                                  const Real nLR_dot_d) const;

  /**
   * Computes left, right, and middle wave speeds
   */
  void computeWaveSpeeds(const ADReal & rhoL,
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
                         ADReal & sm) const;

  /**
   * Computes left and right wave speeds using the Einfeldt formulation
   */
  void computeWaveSpeedsEinfeldt(const ADReal & rhoL,
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
                                 ADReal & sR) const;

  /**
   * Computes left and right wave speeds using the Davis formulation
   */
  void computeWaveSpeedsDavis(const ADReal & unL,
                              const ADReal & cL,
                              const ADReal & unR,
                              const ADReal & cR,
                              ADReal & sL,
                              ADReal & sR) const;

protected:
  /// Type for how to compute left and right wave speeds
  enum class WaveSpeedFormulation
  {
    EINFELDT,
    DAVIS
  };

  /**
   * Computes the flow area that is used in the numerical flux
   */
  virtual ADReal computeFlowArea(const std::vector<ADReal> & UL,
                                 const std::vector<ADReal> & UR) const;

  /**
   * Converts 3D solution vector to 1D solution vector
   *
   * @param[in] 1D solution vector
   */
  std::vector<ADReal> convert3Dto1D(const std::vector<ADReal> & U_3d) const;

  /**
   * Computes primitive variables from 3D solution vector
   */
  void computeREPCUA3D(const std::vector<ADReal> & U_3d,
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
                       ADReal & A) const;

  /**
   * Computes primitive variables from 1D solution vector
   */
  void computeREPCUA1D(const std::vector<ADReal> & U_1d,
                       const Real nLR_dot_d,
                       ADReal & rho,
                       ADReal & e,
                       ADReal & p,
                       ADReal & c,
                       ADReal & un,
                       ADReal & A) const;

  /**
   * Computes the Roe-averaged normal velocity and sound speed
   */
  void computeRoeSpeeds(const ADReal & rhoL,
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
                        ADReal & c_roe) const;

  /**
   * Computes the pressure in left or right subsonic region
   *
   * @param[in] rhoL  Left density
   * @param[in] unL  Left normal velocity
   * @param[in] pL  Left pressure
   * @param[in] sL  Left wave speed
   * @param[in] sm  Middle wave speed
   */
  ADReal starPressure(const ADReal & rhoL,
                      const ADReal & unL,
                      const ADReal & pL,
                      const ADReal & sL,
                      const ADReal & sm) const;

  /// fluid properties user object
  const SinglePhaseFluidProperties & _fp;

  /// How to compute left and right wave speeds
  const WaveSpeedFormulation _wave_speed_formulation;

public:
  static InputParameters validParams();
};
