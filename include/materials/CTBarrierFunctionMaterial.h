//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
#include "CrossTermBarrierFunctionBaseWij.h"

// Forward Declarations
class CTBarrierFunctionMaterial;

class CTBarrierFunctionMaterial : public CrossTermBarrierFunctionBaseWij
{
public:
  CTBarrierFunctionMaterial(const InputParameters & parameters);

  static InputParameters validParams();

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

};
