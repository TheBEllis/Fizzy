#include "ActionWarehouse.h"
#include "AddInventoryManager.h"
#include "AddKernelAction.h"
#include "Attributes.h"
#include "AuxKernel.h"
#include "AuxiliarySystem.h"
#include "FispactInventoryManager.h"
#include "FizzyEnums.h"
#include "InputParameters.h"
#include "MooseMeshUtils.h"
#include "MooseTypes.h"
#include "NuclideMetric.h"
#include "Registry.h"
#include <memory>

registerMooseAction("FizzyApp", AddInventoryManager, "add_inventory_manager");

typedef std::pair<nuclide_quantities::NuclideQuantitiesEnum,
                  std::set<SubdomainID>>
    NuclideMetricsBlocksPair;

InputParameters AddInventoryManager::validParams() {
  InputParameters params = Action::validParams();
  params.addClassDescription(
      "A MOOSE action automating the inclusion of a FispactInventoryManager "
      "when aux kernels and postprocessors requiring inventory quantities are "
      "included in the input file.");
  return params;
}

AddInventoryManager::AddInventoryManager(const InputParameters &params)
    : Action(params) {}

void AddInventoryManager::act() {
  std::unordered_map<std::string, std::vector<NuclideMetricsBlocksPair>>
      nuclide_metrics;

  std::unordered_map<inventory_outputs::InventoryOutputsEnum,
                     std::set<SubdomainID>>
      element_metrics;

  // Get all aux kernels in simulation
  const std::vector<std::shared_ptr<AuxKernel>> aux_kernels =
      _problem->getAuxiliarySystem().elemAuxWarehouse().getObjects();

  // Get all postprocessors in simulation
  std::vector<Postprocessor *> postprocessors;
  _problem->theWarehouse()
      .query()
      .condition<AttribInterfaces>(Interfaces::Postprocessor)
      .queryInto(postprocessors);

  /// Add all auxiliary kernel requests for nuclide/inventory metrics
  for (auto &aux_kernel : aux_kernels) {
    const InputParameters &params = aux_kernel->parameters();

    std::set<SubdomainID> block_ids_set;

    if (aux_kernel->blocks().empty()) {
      block_ids_set = _mesh->meshSubdomains();
    } else {

      std::vector<SubdomainName> block_names = aux_kernel->blocks();

      std::vector<SubdomainID> block_ids =
          MooseMeshUtils::getSubdomainIDs(*_mesh, block_names);

      block_ids_set = std::set<SubdomainID>(block_ids.begin(), block_ids.end());
    }

    if (params.get<std::string>("_type") == "FispactNuclideKernel") {

      const nuclide_quantities::NuclideQuantitiesEnum &metric =
          params.get<MooseEnum>("metric")
              .getEnum<nuclide_quantities::NuclideQuantitiesEnum>();

      const std::string &nuclide = params.get<std::string>("nuclide");

      nuclide_metrics[nuclide].push_back(
          NuclideMetricsBlocksPair(metric, block_ids_set));
    }

    if (params.get<std::string>("_type") == "FispactElementKernel") {

      const inventory_outputs::InventoryOutputsEnum &metric =
          params.get<MooseEnum>("metric")
              .getEnum<inventory_outputs::InventoryOutputsEnum>();

      element_metrics[metric] = block_ids_set;
    }
  }

  for (Postprocessor *postprocessor : postprocessors) {

    if (NuclideMetric *fispact_postprocessor =
            dynamic_cast<NuclideMetric *>(postprocessor)) {

      const InputParameters &params = fispact_postprocessor->parameters();

      std::set<SubdomainID> block_ids_set;

      if (fispact_postprocessor->blocks().empty()) {
        block_ids_set = _mesh->meshSubdomains();
      } else {
        std::vector<SubdomainName> block_names =
            fispact_postprocessor->blocks();

        std::vector<SubdomainID> block_ids =
            MooseMeshUtils::getSubdomainIDs(*_mesh, block_names);

        block_ids_set =
            std::set<SubdomainID>(block_ids.begin(), block_ids.end());
      }

      const nuclide_quantities::NuclideQuantitiesEnum &metric =
          params.get<MooseEnum>("metric")
              .getEnum<nuclide_quantities::NuclideQuantitiesEnum>();

      const std::vector<std::string> &nuclides =
          params.get<std::vector<std::string>>("nuclides");

      for (const std::string &nuclide : nuclides) {

        nuclide_metrics[nuclide].push_back(
            NuclideMetricsBlocksPair(metric, block_ids_set));
      }
    }
  }

  if (!(nuclide_metrics.empty() && element_metrics.empty())) {
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

    for (auto &[metric, blocks] : element_metrics) {

      inv_manager.registerElementMetricRequest(metric, blocks);
    }
  }
}
