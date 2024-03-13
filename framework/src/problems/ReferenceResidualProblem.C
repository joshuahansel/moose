//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "ReferenceResidualProblem.h"
#include "ReferenceResidualConvergence.h"

#include "AuxiliarySystem.h"
#include "MooseApp.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "MooseVariableScalar.h"
#include "NonlinearSystem.h"

// libMesh includes
#include "libmesh/enum_norm_type.h"
#include "libmesh/utility.h"

registerMooseObject("MooseApp", ReferenceResidualProblem);

InputParameters
ReferenceResidualProblem::validParams()
{
  InputParameters params = FEProblem::validParams();
  params += ReferenceResidualConvergence::validCommonReferenceResidualProblemParams();
  params.addClassDescription("Problem that checks for convergence relative to "
                             "a user-supplied reference quantity rather than "
                             "the initial residual");

  return params;
}

ReferenceResidualProblem::ReferenceResidualProblem(const InputParameters & params)
  : FEProblem(params),
    _use_group_variables(false),
    _reference_vector(nullptr),
    _converge_on(getParam<std::vector<NonlinearVariableName>>("converge_on")),
    _zero_ref_type(
        params.get<MooseEnum>("zero_reference_residual_treatment").getEnum<ZeroReferenceType>()),
    _reference_vector_tag_id(Moose::INVALID_TAG_ID)
{
  // mooseError("The `ReferenceResidualProblem` class is deprecated and will be removed "
  //" please the convergence class instead.");
  if (params.isParamValid("solution_variables"))
  {
    if (params.isParamValid("reference_vector"))
      mooseDeprecated("The `solution_variables` parameter is deprecated, has no effect when "
                      "the tagging system is used, and will be removed on January 1, 2020. "
                      "Please simply delete this parameter from your input file.");
    _soln_var_names = params.get<std::vector<NonlinearVariableName>>("solution_variables");
  }

  if (params.isParamValid("reference_residual_variables") &&
      params.isParamValid("reference_vector"))
    mooseError(
        "For `ReferenceResidualProblem` you can specify either the `reference_residual_variables` "
        "or `reference_vector` parameter, not both. `reference_residual_variables` is deprecated "
        "so we recommend using `reference_vector`");

  if (params.isParamValid("reference_residual_variables"))
  {
    mooseDeprecated(
        "The save-in method for composing reference residual quantities is deprecated "
        "and will be removed on January 1, 2020. Please use the tagging system instead; "
        "specifically, please assign a TagName to the `reference_vector` parameter");

    _ref_resid_var_names = params.get<std::vector<AuxVariableName>>("reference_residual_variables");

    if (_soln_var_names.size() != _ref_resid_var_names.size())
      mooseError("In ReferenceResidualProblem, size of solution_variables (",
                 _soln_var_names.size(),
                 ") != size of reference_residual_variables (",
                 _ref_resid_var_names.size(),
                 ")");
  }
  else if (params.isParamValid("reference_vector"))
  {
    if (numNonlinearSystems() > 1)
      paramError(
          "nl_sys_names",
          "reference residual problem does not currently support multiple nonlinear systems");
    _reference_vector_tag_id = getVectorTagID(getParam<TagName>("reference_vector"));
  }
  else
    mooseInfo("Neither the `reference_residual_variables` nor `reference_vector` parameter is "
              "specified for `ReferenceResidualProblem`, which means that no reference "
              "quantites are set. Because of this, the standard technique of comparing the "
              "norm of the full residual vector with its initial value will be used.");

  if (params.isParamValid("group_variables"))
  {
    _group_variables =
        params.get<std::vector<std::vector<NonlinearVariableName>>>("group_variables");
    _use_group_variables = true;
  }

  _accept_mult = params.get<Real>("acceptable_multiplier");
  _accept_iters = params.get<int>("acceptable_iterations");

  const auto norm_type_enum =
      params.get<MooseEnum>("normalization_type").getEnum<NormalizationType>();
  setNormType(norm_type_enum);
  if (norm_type_enum == NormalizationType::LOCAL_L2)
  {
    _norm_type = DISCRETE_L2;
    _local_norm = true;
  }
  else if (norm_type_enum == NormalizationType::GLOBAL_L2)
  {
    _norm_type = DISCRETE_L2;
    _local_norm = false;
  }
  else if (norm_type_enum == NormalizationType::LOCAL_LINF)
  {
    _norm_type = DISCRETE_L_INF;
    _local_norm = true;
  }
  else if (norm_type_enum == NormalizationType::GLOBAL_LINF)
  {
    _norm_type = DISCRETE_L_INF;
    _local_norm = false;
  }
  else
    mooseError("Internal error");

  if (_local_norm && !params.isParamValid("reference_vector"))
    paramError("reference_vector", "If local norm is used, a reference_vector must be provided.");
}

void
ReferenceResidualProblem::addDefaultConvergence()
{
  const std::string class_name = "ReferenceResidualConvergence";
  InputParameters params = _factory.getValidParams(class_name);
  params.applyParameters(parameters());
  addConvergence(class_name, _nonlinear_convergence_name, params);
}
