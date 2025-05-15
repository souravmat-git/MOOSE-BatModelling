//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "DerivativeMaterialInterface.h"

// Forward Declarations

/**
 * CrossTermBarrierFunctionBaseWij is the base to a set of free energy penalties that
 * set the phase interface barriers for arbitrary pairs of phases.
 */
class CrossTermBarrierFunctionBaseWij : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  CrossTermBarrierFunctionBaseWij(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  /// name of the function of eta (used to generate the material property names)
  std::string _function_name;

  /// polynomial order of the switching function \f$ g(\eta) \f$
  unsigned int _g_order;

  /// order parameters
  unsigned int _num_eta;
  const std::vector<VariableName> _eta_names;
  const std::vector<const VariableValue *> _eta;

  ///@{ Barrier function and its derivatives
  MaterialProperty<Real> & _prop_g;
  std::vector<MaterialProperty<Real> *> _prop_dg;
  std::vector<std::vector<MaterialProperty<Real> *>> _prop_d2g;
  ///@}

};
