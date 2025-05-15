//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CTBarrierFunctionMaterial.h"
registerMooseObject("yourApp", CTBarrierFunctionMaterial);

//template<>
InputParameters
CTBarrierFunctionMaterial::validParams()
{
  InputParameters params = CrossTermBarrierFunctionBaseWij::validParams();
  return params;
}

CTBarrierFunctionMaterial::CTBarrierFunctionMaterial(const InputParameters & parameters)
 :CrossTermBarrierFunctionBaseWij(parameters)
 {
   // error out if W_ij is not symmetric
   //for (unsigned int i = 0; i < _num_eta; ++i)
   //for (unsigned int j = 0; j < _num_eta; ++j)
   // {
  //    if (i==j)
  //        _W_ij[_num_eta * i + j] = 0;
  //      else
  //        _W_ij[_num_eta * i + j] = 1;
  //  }


 }

void
CTBarrierFunctionMaterial::computeQpProperties()
{
  // Initialize properties to zero before accumulating
  CrossTermBarrierFunctionBaseWij::computeQpProperties();

  // Sum the components of our W_ij matrix to get constant used in our g function
  for (unsigned int i = 0; i < _num_eta; ++i)
    for (unsigned int j = i + 1; j < _num_eta; ++j)
    {
      const Real ni = (*_eta[i])[_qp];
      const Real nj = (*_eta[j])[_qp];

      switch (_g_order)
      {
        case 0: // SIMPLE
          _prop_g[_qp] += 16.0 * (ni * ni * nj * nj);
          // first derivatives
          (*_prop_dg[i])[_qp] += 16.0 * (2 * ni * nj * nj);
          (*_prop_dg[j])[_qp] += 16.0 * (2 * ni * ni * nj);
          // second derivatives (diagonal)
          (*_prop_d2g[i][i])[_qp] += 16.0 * (2 * nj * nj);
          (*_prop_d2g[j][j])[_qp] += 16.0 * (2 * ni * ni);
          // second derivatives (off-diagonal)
          (*_prop_d2g[i][j])[_qp] = 16.0  * (4 * ni * nj);
          break;

        case 1: // LOW
          _prop_g[_qp] += 4.0  * (ni * nj);
          // first derivatives
          (*_prop_dg[i])[_qp] += 4.0 * nj;
          (*_prop_dg[j])[_qp] += 4.0  * ni;
          // second derivatives (diagonal) vanish
          // second derivatives (off-diagonal)
          (*_prop_d2g[i][j])[_qp] = 4.0 ;
          break;

        default:
          mooseError("Internal error");
      }
    }
}
