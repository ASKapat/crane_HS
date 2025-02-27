/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#pragma once

#include "ADKernel.h"

class ADEEDFEnergyLog_HS : public ADKernel
{
public:
  static InputParameters validParams();
  ADEEDFEnergyLog_HS(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual();

  
  std::string _reaction_name;
  std::string _reaction_coeff_name;
  
 

  const ADMaterialProperty<Real> & _reaction_coefficient;
  const ADVariableValue & _em;
  const ADVariableValue & _target;
  const ADVariableValue & _threshold_energy;
  const Real _rxn_sgn;
  // Threshold energy is just a parameter generally, though elastic collisions require a material property.

};
