#include "ActionWarehouse.h"
#include "AddAuxKernelAction.h"
#include "AddInventoryManager.h"
#include "FispactInventoryManager.h"
#include "InputParameters.h"
#include "MooseMeshUtils.h"
#include "MooseTypes.h"
#include "Registry.h"

registerMooseAction("FizzyApp", AddInventoryManager,
                    "add_vector_postprocessor");

typedef std::pair<nuclide_quantities::NuclideQuantitiesEnum,
                  std::set<SubdomainID>>
    NuclideMetricsBlocksPair;

InputParameters AddInventoryManager::validParams() {
  InputParameters params = Action::validParams();
  return params;
}

AddInventoryManager::AddInventoryManager(const InputParameters &params)
    : Action(params) {}

void AddInventoryManager::act() {

  std::unordered_map<std::string, std::vector<NuclideMetricsBlocksPair>>
      nuclide_metrics;

  std::vector<const AddAuxKernelAction *> aux_input_blocks =
      _awh.getActions<AddAuxKernelAction>();

  for (auto &aux_kernel_action : aux_input_blocks) {
    const InputParameters &params = aux_kernel_action->parameters();
    if (params.get<std::string>("type") == "FispactNuclideKernel") {

      const nuclide_quantities::NuclideQuantitiesEnum &metric =
          params.get<nuclide_quantities::NuclideQuantitiesEnum>("metric");

      const std::string &nuclide = params.get<std::string>("nuclide");

      std::vector<SubdomainName> block_names =
          params.get<std::vector<SubdomainName>>("blocks");

      std::vector<SubdomainID> block_ids =
          MooseMeshUtils::getSubdomainIDs(*_mesh, block_names);

      std::set<SubdomainID> block_ids_set(block_ids.begin(), block_ids.end());

      nuclide_metrics[nuclide].push_back(
          NuclideMetricsBlocksPair(metric, block_ids_set));
    }
  }

  if (!nuclide_metrics.empty()) {
    InputParameters inventory_params =
        _app.getFactory().getValidParams("FispactInventoryManager");
    _problem->addUserObject("FispactInventoryManager", "inv_manager",
                            inventory_params);

    FispactInventoryManager &inv_manager =
        _problem->getUserObject<FispactInventoryManager>("inv_manager");

    for (auto &[nuclide, metric_block_pairs] : nuclide_metrics) {
      for (NuclideMetricsBlocksPair &metric_block_pair : metric_block_pairs) {

        inv_manager.registerNuclideMetricRequest(
            nuclide, metric_block_pair.first, metric_block_pair.second);
      }
    }
  }
}
