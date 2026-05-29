#include "FispactElementPostprocessor.h"
#include "FispactInventoryManager.h"
#include "FispactNuclideMetric.h"
#include "FispactProblem.h"
#include "FizzyEnums.h"
#include "UserObject.h"
#include <numeric>

registerMooseObject("FizzyApp", FispactNuclideMetric);

InputParameters FispactNuclideMetric::validParams() {

  InputParameters params = FispactElementPostprocessor::validParams();

  params.addClassDescription(
      "Postprocessor to calculate a user defined inventory "
      "metric, for a given nuclide(s) in the inventory, over the entire "
      "domain.");

  params.addRequiredParam<std::vector<std::string>>(
      "nuclides", "Which nuclide should this kernel track");

  params.addRequiredParam<MooseEnum>("metric", getNuclideMetricsEnum(),
                                     "Which metric should be tracked");
  return params;
}

FispactNuclideMetric::FispactNuclideMetric(const InputParameters &params)
    : FispactElementPostprocessor(params),
      _nuclides(getParam<std::vector<std::string>>("nuclides")),
      _metric(getParam<MooseEnum>("metric")
                  .getEnum<nuclide_quantities::NuclideQuantitiesEnum>()) {}

void FispactNuclideMetric::initialize() { _sum = 0; }

void FispactNuclideMetric::threadJoin(const UserObject &y) {};

void FispactNuclideMetric::execute() {

  size_t inv_index = getFispactProblem().getFispactInventoryIndexFromTime();

  const FispactInventoryManager &inv_manager =
      getUserObjectByName<FispactInventoryManager>("inv_manager");

  for (std::string &nuclide : _nuclides) {

    _sum += inv_manager.getNuclideMetric(_current_elem->id(), nuclide,
                                         inv_index, _metric);
  }
}

void FispactNuclideMetric::finalize() { comm().sum(_sum); }

PostprocessorValue FispactNuclideMetric::getValue() const { return _sum; }
