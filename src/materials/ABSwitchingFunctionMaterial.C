//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ABSwitchingFunctionMaterial.h"
registerMooseObject("yourApp", ABSwitchingFunctionMaterial);

//template<>
InputParameters
ABSwitchingFunctionMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<MaterialPropertyName>(
      "h_name", "Name of the switching function material property for the given phase");
  params.addRequiredCoupledVarWithAutoBuild(
    "phase_etas", "var_name_base", "op_num", "Supply the order parameter representing the grains");
  params.addRequiredCoupledVarWithAutoBuild("all_etas","var_name_base","eta_num","Vector of all order parameters for all phases");
  return params;
}

ABSwitchingFunctionMaterial::ABSwitchingFunctionMaterial(const InputParameters & parameters)
 :SwitchingFunctionMultiPhaseMaterial(parameters)
 {
 }

void
ABSwitchingFunctionMaterial::computeQpProperties()
{
  SwitchingFunctionMultiPhaseMaterial::computeQpProperties();
}
