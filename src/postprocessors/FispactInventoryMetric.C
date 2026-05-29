#include "FispactElementPostprocessor.h"
#include "FispactInventoryManager.h"
#include "FispactInventoryMetric.h"
#include "FispactProblem.h"
#include "FizzyEnums.h"
#include "UserObject.h"
#include <numeric>

registerMooseObject("FizzyApp", FispactInventoryMetric);

InputParameters FispactInventoryMetric::validParams() {

  InputParameters params = FispactElementPostprocessor::validParams();

  params.addClassDescription(
      "Postprocessor to calculate a user defined inventory "
      "metric, for all nuclides in the inventory, over the entire "
      "domain (or specified subdomains).");

  params.addRequiredParam<MooseEnum>("metric", getInventoryMetricsEnum(),
                                     "Which metric should be tracked");
  return params;
}

FispactInventoryMetric::FispactInventoryMetric(const InputParameters &params)
    : FispactElementPostprocessor(params),
      _metric(getParam<MooseEnum>("metric")
                  .getEnum<inventory_outputs::InventoryOutputsEnum>()) {}

void FispactInventoryMetric::initialize() { _sum = 0; }

void FispactInventoryMetric::threadJoin(const UserObject &y) {};

void FispactInventoryMetric::execute() {

  size_t inv_index = getFispactProblem().getFispactInventoryIndexFromTime();

  const FispactInventoryManager &inv_manager =
      getUserObjectByName<FispactInventoryManager>("inv_manager");

  _sum += inv_manager.getElementMetric(_current_elem->id(), inv_index, _metric);
}

void FispactInventoryMetric::finalize() { comm().sum(_sum); }

PostprocessorValue FispactInventoryMetric::getValue() const { return _sum; }
