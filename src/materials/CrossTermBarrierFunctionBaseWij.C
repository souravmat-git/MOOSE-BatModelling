//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CrossTermBarrierFunctionBaseWij.h"
registerMooseObject("sequoiaApp", CrossTermBarrierFunctionBaseWij);


InputParameters
CrossTermBarrierFunctionBaseWij::validParams()
{
  InputParameters params = Material::validParams();
  params.addParam<std::string>("function_name", "g", "actual name for g(eta_i)");
  MooseEnum g_order("SIMPLE=0 LOW", "SIMPLE");
  params.addParam<MooseEnum>("g_order", g_order, "Polynomial order of the barrier function g(eta)");
  params.addRequiredCoupledVarWithAutoBuild("etas","var_name_base", "op_num" , "eta_i order parameters, one for each h");
  return params;
}

CrossTermBarrierFunctionBaseWij::CrossTermBarrierFunctionBaseWij(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _function_name(getParam<std::string>("function_name")),
    _g_order(getParam<MooseEnum>("g_order")),
    _num_eta(coupledComponents("etas")),
    _eta_names(coupledNames("etas")),
    _eta(coupledValues("etas")),
    _prop_g(declareProperty<Real>(_function_name)),
    _prop_dg(_num_eta),
    _prop_d2g(_num_eta)
{

  // declare g derivative properties, fetch eta values
  for (unsigned int i = 0; i < _num_eta; ++i)
    _prop_d2g[i].resize(_num_eta);

  for (unsigned int i = 0; i < _num_eta; ++i)
  {
    _prop_dg[i] = &declarePropertyDerivative<Real>(_function_name, _eta_names[i]);
    for (unsigned int j = i; j < _num_eta; ++j)
    {
      _prop_d2g[i][j] = _prop_d2g[j][i] =
          &declarePropertyDerivative<Real>(_function_name, _eta_names[i], _eta_names[j]);
    }
  }
}

void
CrossTermBarrierFunctionBaseWij::computeQpProperties()
{
  // Initialize properties to zero before accumulating
  _prop_g[_qp] = 0.0;
  for (unsigned int i = 0; i < _num_eta; ++i)
  {
    (*_prop_dg[i])[_qp] = 0.0;
    (*_prop_d2g[i][i])[_qp] = 0.0;
  }
}
