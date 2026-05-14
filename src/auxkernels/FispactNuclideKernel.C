#include "FispactNuclideKernel.h"
#include "FispactProblem.h"
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

Real FispactNuclideKernel::computeValue() {
  const FispactInventoryManager &inv_manager =
      getUserObjectByName<FispactInventoryManager>("inv_manager");

  FispactProblem &fis_problem = static_cast<FispactProblem &>(_c_fe_problem);

  size_t inventory_idx = fis_problem.getFispactInventoryIndexFromTime();

  return inv_manager.getNuclideMetric(_current_elem->id(), _nuclide,
                                      inventory_idx, _metric);
}
