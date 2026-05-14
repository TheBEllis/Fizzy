#include "FispactElementKernel.h"
#include "FispactProblem.h"
#include "FizzyEnums.h"
#include "InputParameters.h"
#include "Registry.h"

registerMooseObject("FizzyApp", FispactElementKernel);

InputParameters FispactElementKernel::validParams() {
  InputParameters params = FispactAuxKernel::validParams();

  params.addRequiredParam<MooseEnum>("metric", getInventoryMetricsEnum(),
                                     "Which metric should be tracked.");

  return params;
}

FispactElementKernel::FispactElementKernel(const InputParameters &params)
    : FispactAuxKernel(params),
      _metric(getParam<MooseEnum>("metric")
                  .getEnum<inventory_outputs::InventoryOutputsEnum>()) {}

Real FispactElementKernel::computeValue() {
  const FispactInventoryManager &inv_manager =
      getUserObjectByName<FispactInventoryManager>("inv_manager");

  FispactProblem &fis_problem = static_cast<FispactProblem &>(_c_fe_problem);

  size_t inventory_idx = fis_problem.getFispactInventoryIndexFromTime();

  return inv_manager.getElementMetric(_current_elem->id(), inventory_idx,
                                      _metric);
}
