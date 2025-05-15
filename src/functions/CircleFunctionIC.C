//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CircleFunctionIC.h"
registerMooseObject("yourApp", CircleFunctionIC);

InputParameters
CircleFunctionIC::validParams()
{
  InputParameters params = Function::validParams();
  params.addParam<Real>("inside_value", 1.0, "Inside value");
  params.addParam<Real>("outside_value", 0.0,  "Outside value");
  params.addRequiredParam<FileName>("table_name", "Table with x y r");
  return params;
}

CircleFunctionIC::CircleFunctionIC(const InputParameters & parameters)
  : Function(parameters),
   _inside_val(getParam<Real>("inside_value")),
   _outside_val(getParam<Real>("outside_value")),
   _table_name(getParam<FileName>("table_name")),
   _table_reader(_table_name, &_communicator),
   _col_names(_table_reader.getNames())
{

  //Lines begining with # are comments
  _table_reader.setComment("#");
  //tab is a delimiter
  _table_reader.setDelimiter(" ");

  std::ifstream file(_table_name.c_str());
  if (file.good()){
    _console << "Reading the file.... \n" << _table_name << "\n";
    //read the data in the table and convert into double
    _table_reader.read();

    //Check the first column is always mole_fraction of component B
      if (_col_names[0] != "x")
        mooseError("First column must be x-coordinates !!");
   }

    _x = _table_reader.getData(_col_names[0]);
    _y = _table_reader.getData(_col_names[1]);
    _radii = _table_reader.getData(_col_names[2]);


  //Check if the input is equal to x
  const Real num = _radii.size();
  if (num != _x.size())
    mooseError("Try again");

}

Real
CircleFunctionIC::value(Real /*t*/, const Point & pt) const
{
    Real x = pt(0);
    Real y = pt(1);

    Real _value = _outside_val;

    for (unsigned i = 0; i< _radii.size(); i++){
     if ( std::pow((x -_x[i]),2.0) +  std::pow((y - _y[i]),2.0) <= _radii[i]*_radii[i])
         _value =  _inside_val;
    }
    return _value;
}
