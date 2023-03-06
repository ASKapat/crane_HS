#include "ADEEDFEnergyLog_HS.h"

// MOOSE includes
#include "MooseUtils.h"
#include "MooseVariable.h"

registerADMooseObject("CraneApp", ADEEDFEnergyLog_HS);

InputParameters
ADEEDFEnergyLog_HS::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addRequiredCoupledVar("electrons", "The electron density.");
  params.addRequiredCoupledVar("target", "The target species.");
  params.addRequiredCoupledVar("threshold_energy", "Energy required for reaction to take place.");
  params.addRequiredParam<std::string>("reaction", "Stores the full reaction equation.");
  params.addRequiredParam<Real>("gain_loss", "Whether reaction is gain or loss for species");
  
  params.addParam<std::string>(
      "number",
      "",
      "The reaction number. Optional, just for material property naming purposes. If a single "
      "reaction has multiple different rate coefficients (frequently the case when multiple "
      "species are lumped together to simplify a reaction network), this will prevent the same "
      "material property from being declared multiple times.");
  return params;
}

ADEEDFEnergyLog_HS::ADEEDFEnergyLog_HS(const InputParameters & parameters)
  : ADKernel(parameters),
    _reaction_name(getParam<std::string>("reaction")),
    _reaction_coefficient(getADMaterialProperty<Real>("k" + getParam<std::string>("number") + "_" +
                                                      getParam<std::string>("reaction"))),
    _em(adCoupledValue("electrons")),
    _target(adCoupledValue("target")),
    _threshold_energy(adCoupledValue("threshold_energy")),
    _rxn_sgn(getParam<Real>("gain_loss"))
{
}

ADReal
ADEEDFEnergyLog_HS::computeQpResidual()
{       
  return -_test[_i][_qp] *_rxn_sgn *_reaction_coefficient[_qp] * (2./3.) * std::exp(_em[_qp] + _threshold_energy[_qp]);
}
