#include "FispactContextBase.h"
#include "FispactInventoryManager.h"
#include "FizzyEnums.h"
#include <memory>
registerMooseObject("FizzyApp", FispactInventoryManager);

FispactInventoryManager::FispactInventoryManager(
    const InputParameters &parameters)
    : FispactUserObject(parameters) {

  _n_fispact_inventories = getFispactProblem().getSchedule()->getTimes().size();
}

void FispactInventoryManager::registerNuclideMetrics(
    std::string &nuclide,
    std::vector<nuclide_quantities::NuclideQuantitiesEnum> metric,
    std::set<SubdomainID> &blocks) {

  if (_nuclide_ids.find(nuclide) == _nuclide_ids.end()) {

    _nuclide_ids[nuclide] = ++_max_nuclide_id;

    _nuclide_names[_nuclide_ids[nuclide]] = nuclide;
  }

  for (const libMesh::Elem *elem :
       *getSubProblem().mesh().getActiveLocalElementRange()) {
    if (blocks.count(elem->subdomain_id())) {
      getElementInventory(elem->id())
          .insertNuclideMetricRequest(nuclide, _nuclide_ids[nuclide],
                                      _n_fispact_inventories, metric);
    }
  }
}

void FispactInventoryManager::extractInventoryData(
    FispactContextBase &fp_context, size_t elem_id, int inv_index) {

  std::vector<std::unique_ptr<FispactOutputNuclideDataBase>>
      fispact_inventory_nuclides =
          fp_context.getOutput().getInventoryNuclides(inv_index);

  ElementInventory &elem_inv = getElementInventory(elem_id);

  for (size_t &nuclide : elem_inv.getNuclides()) {

    NuclideInventory &nuclide_inv = elem_inv.getNuclide(nuclide);

    int zai = fp_context.getUtils().GetZai(_nuclide_names[nuclide]);
    // Check nuclide data exists for this inventory step
    if (fp_context.getOutput().findInventoryExists(inv_index, zai)) {

      int fispact_nuclide_index =
          fp_context.getOutput().findInventoryIndex(inv_index, zai);

      for (nuclide_quantities::NuclideQuantitiesEnum quantity :
           elem_inv.getNuclide(nuclide).getQuantities()) {

        nuclide_inv.getQuantity(quantity) =
            fispact_inventory_nuclides[fispact_nuclide_index]->getQuantity(
                quantity);
      }
    }
  }
}

double FispactInventoryManager::getNuclideMetric(
    libMesh::dof_id_type elem_id, std::string &nuclide, int inv_index,
    nuclide_quantities::NuclideQuantitiesEnum metric) const {
  return _element_inventories[elem_id].getNuclide(_nuclide_ids[nuclide]);
}

void FispactInventoryManager::initialiseNuclearInventory() {
  for (const libMesh::Elem *elem :
       *getSubProblem().mesh().getActiveLocalElementRange()) {

    _element_inventories.emplace_back(elem->id());
  }
}

FispactInventoryManager::ElementInventory &
FispactInventoryManager::getElementInventory(size_t elem_id) {
  for (int i = 0; i < _element_inventories.size(); i++) {
    if (_element_inventories[i].getElemID() == elem_id) {
      return _element_inventories[i];
    }
  }
}
