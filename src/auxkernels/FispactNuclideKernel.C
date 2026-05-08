#include "FispactNuclideKernel.h"
#include "InputParameters.h"
#include "Registry.h"

registerMooseObject("FizzyApp", FispactNuclideKernel);

InputParameters FispactNuclideKernel::validParams() {
  InputParameters params = FispactAuxKernel::validParams();

  params.addRequiredParam<std::string>(
      "nuclide", "Which nuclide should this kernel track.");

  params.addRequiredParam<MooseEnum>("metric", getNuclideMetricsEnum(),
                                     "Which metric should be tracked.");

  return params;
}

FispactNuclideKernel::FispactNuclideKernel(const InputParameters &params)
    : FispactAuxKernel(params), _nuclide(getParam<std::string>("nuclide")),
      _metric(getParam<MooseEnum>("metric")
                  .getEnum<nuclide_quantities::NuclideQuantitiesEnum>()) {}

Real FispactNuclideKernel::computeValue() { return 1; }
