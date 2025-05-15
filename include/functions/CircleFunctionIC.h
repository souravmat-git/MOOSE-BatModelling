//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

//MOOSE includes
#include "Function.h"
#include "DelimitedFileReader.h"

class CircleFunctionIC;

//template<>
//InputParameters validParams<CircleFunctionIC>();

class CircleFunctionIC : public Function
{
  public:

    static InputParameters validParams();

    CircleFunctionIC(const InputParameters & parameters);

    //This member function is re-dedfined  from the base class
    virtual Real value(Real t, const Point & p) const override;

   protected:
    const Real _inside_val;
    const Real _outside_val;

    //variable to hold the string type table name from input file
    FileName _table_name;

    //Moose utility csv_reader defined for reading the table
    //Note that DelimitedFileReader is a class defined within

    MooseUtils::DelimitedFileReader _table_reader;
    const std::vector<std::string>& _col_names;

    std::vector<Real> _radii;  //Radius of the circle;
    std::vector<Real> _x;     // x-Coordinates
    std::vector<Real> _y;     // y-coordinates
};
