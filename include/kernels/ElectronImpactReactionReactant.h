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

#ifndef ELECTRONIMPACTREACTIONREACTANT_H
#define ELECTRONIMPACTREACTIONREACTANT_H

#include "Kernel.h"

class ElectronImpactReactionReactant;

template <>
InputParameters validParams<ElectronImpactReactionReactant>();

class ElectronImpactReactionReactant : public Kernel
{
public:
  ElectronImpactReactionReactant(const InputParameters & parameters);
  virtual ~ElectronImpactReactionReactant();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  Real _r_units;
  std::string _reaction_coeff_name;
  std::string _reaction_name;

  const MaterialProperty<Real> & _diffem;
  const MaterialProperty<Real> & _muem;
  const MaterialProperty<Real> & _alpha;
  const MaterialProperty<Real> & _d_iz_d_actual_mean_en;
  const MaterialProperty<Real> & _d_muem_d_actual_mean_en;
  const MaterialProperty<Real> & _d_diffem_d_actual_mean_en;

  const VariableValue & _mean_en;
  const VariableGradient & _grad_potential;
  const VariableValue & _em;
  const VariableGradient & _grad_em;
  unsigned int _mean_en_id;
  unsigned int _potential_id;
  unsigned int _em_id;
  const VariableValue & _target;
  unsigned int _target_id;
};

#endif /* ELECTRONIMPACTREACTIONREACTANT_H */
