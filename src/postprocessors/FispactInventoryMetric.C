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

void FispactInventoryMetric::initialize() {
  _sum = 0;
  _total_mass = 0;
}

void FispactInventoryMetric::threadJoin(const UserObject &y) {};

void FispactInventoryMetric::execute() {

  size_t inv_index = getFispactProblem().getFispactInventoryIndexFromTime();

  const FispactInventoryManager &inv_manager =
      getUserObjectByName<FispactInventoryManager>("inv_manager");

  double metric_value =
      inv_manager.getElementMetric(_current_elem->id(), inv_index, _metric);

  if (_metric == inventory_outputs::INVENTORY_DOSE_RATE) {

    double element_mat_density = getFispactProblem()
                                     .getElementMaterial(_current_elem->id())
                                     .getDensity();

    double element_mass = _current_elem_volume * element_mat_density;
    metric_value *= element_mass;
    _total_mass += element_mass;
  }

  _sum += metric_value;
}

void FispactInventoryMetric::finalize() {

  if (_metric == inventory_outputs::INVENTORY_DOSE_RATE) {
    comm().sum(_total_mass);

    // This changes our dose values back to sieverts
    _sum /= _total_mass;
  }
  comm().sum(_sum);
}

PostprocessorValue FispactInventoryMetric::getValue() const { return _sum; }
