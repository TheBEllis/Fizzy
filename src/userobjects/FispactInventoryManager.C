#include "FispactContextBase.h"
#include "FispactInventoryManager.h"
#include "FispactProblem.h"
#include "FizzyEnums.h"
#include <memory>
#include <stdexcept>
#include <string>
registerMooseObject("FizzyApp", FispactInventoryManager);

FispactInventoryManager::FispactInventoryManager(
    const InputParameters &parameters)
    : FispactUserObject(parameters),
      _n_fispact_inventories(
          getFispactProblem().getSchedule().getNumInventories()) {

  initialiseNuclearInventory();
}

void FispactInventoryManager::registerNuclideMetricRequest(
    const std::string &nuclide,
    const nuclide_quantities::NuclideQuantitiesEnum metric,
    const std::set<SubdomainID> &blocks) {

  if (_nuclide_ids.find(nuclide) == _nuclide_ids.end()) {

    _nuclide_ids[nuclide] = ++_max_nuclide_id;

    _nuclide_names[_nuclide_ids[nuclide]] = nuclide;
  }

  for (const libMesh::Elem *elem :
       *getSubProblem().mesh().getActiveLocalElementRange()) {
    if (blocks.count(elem->subdomain_id())) {
      getElementInventory(elem->id())
          .registerNuclideMetricRequest(nuclide, _nuclide_ids[nuclide],
                                        _n_fispact_inventories, metric);
    }
  }
}

void FispactInventoryManager::registerElementMetricRequest(
    const inventory_outputs::InventoryOutputsEnum metric,
    const std::set<SubdomainID> &blocks) {

  for (const libMesh::Elem *elem :
       *getSubProblem().mesh().getActiveLocalElementRange()) {
    if (blocks.count(elem->subdomain_id())) {
      getElementInventory(elem->id())
          .registerElementMetricRequest(_n_fispact_inventories, metric);
    }
  }
}

void FispactInventoryManager::extractInventoryData(
    FispactContextBase &fp_context, size_t elem_id) {

  for (size_t inv_index = 0; inv_index < _n_fispact_inventories; inv_index++) {

    std::vector<std::unique_ptr<FispactOutputNuclideDataBase>>
        fispact_inventory_nuclides =
            fp_context.getOutput().getInventoryNuclides(inv_index);

    ElementInventory &elem_inv = getElementInventory(elem_id);

    for (auto &[metric, values] : elem_inv.getInventoryMetricMap()) {
      values[inv_index] =
          fp_context.getOutput().getInventoryValue(inv_index, metric);
    }

    for (size_t &nuclide : elem_inv.getNuclides()) {

      NuclideInventory &nuclide_inv = elem_inv.getNuclide(nuclide);

      int zai = fp_context.getUtils().GetZai(_nuclide_names[nuclide]);
      // Check nuclide data exists for this inventory step
      int fispact_nuclide_index =
          fp_context.getOutput().findInventoryExists(inv_index, zai)
              ?

              fp_context.getOutput().findInventoryIndex(inv_index, zai)
              : -1;

      // If the FISPACT inventory does not contain this nuclide, set the
      // quantity value to 0
      if (fispact_nuclide_index != -1) {

        for (nuclide_quantities::NuclideQuantitiesEnum quantity :
             elem_inv.getNuclide(nuclide).getQuantities()) {

          nuclide_inv.setQuantity(
              quantity, inv_index,
              fispact_inventory_nuclides[fispact_nuclide_index]->getQuantity(
                  quantity));
        }
      }
    }
  }
}

double FispactInventoryManager::getNuclideMetric(
    libMesh::dof_id_type elem_id, std::string &nuclide, int inv_index,
    nuclide_quantities::NuclideQuantitiesEnum metric) const {
  return getElementInventory(elem_id)
      .getNuclide(_nuclide_ids.at(nuclide))
      .getQuantity(metric, inv_index);
}

double FispactInventoryManager::getElementMetric(
    libMesh::dof_id_type elem_id, int inv_index,
    inventory_outputs::InventoryOutputsEnum metric) const {
  return getElementInventory(elem_id).getElementMetric(metric, inv_index);
}

void FispactInventoryManager::initialiseNuclearInventory() {
  for (const libMesh::Elem *elem :
       *getSubProblem().mesh().getActiveLocalElementRange()) {

    _element_inventories.emplace_back(elem->id());
  }
}

FispactInventoryManager::ElementInventory &
FispactInventoryManager::getElementInventory(size_t elem_id) {

  _console << elem_id << std::endl;
  size_t index = _fispact_problem.getLocalElemIndexMap()[elem_id];
  _console << index << std::endl;
  return _element_inventories[index];
}

const FispactInventoryManager::ElementInventory &
FispactInventoryManager::getElementInventory(size_t elem_id) const {

  _console << elem_id << std::endl;
  size_t index = _fispact_problem.getLocalElemIndexMap()[elem_id];
  return _element_inventories[index];
}
