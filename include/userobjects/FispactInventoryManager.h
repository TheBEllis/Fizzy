#pragma once

#include "FispactContextBase.h"
#include "FispactUserObject.h"
#include "FizzyEnums.h"
#include "InputParameters.h"
#include "UserObject.h"
#include "libmesh/id_types.h"
#include <stdexcept>
#include <unordered_map>

class FispactInventoryManager : public FispactUserObject {
public:
  /// Forward declaration
  class ElementInventory;

  static InputParameters validParams() {

    InputParameters params = FispactUserObject::validParams();

    return params;
  }

  FispactInventoryManager(const InputParameters &parameters);

  virtual void initialize() {}
  virtual void finalize() {}
  virtual void execute() {}

  void registerNuclideMetricRequest(
      const std::string &nuclide,
      const nuclide_quantities::NuclideQuantitiesEnum metric,
      const std::set<SubdomainID> &blocks);

  void registerElementMetricRequest(
      const inventory_outputs::InventoryOutputsEnum metric,
      const std::set<SubdomainID> &blocks);

  void extractInventoryData(FispactContextBase &fp_context, size_t elem_id);

  double
  getNuclideMetric(libMesh::dof_id_type elem_id, std::string &nuclide,
                   int inv_index,
                   nuclide_quantities::NuclideQuantitiesEnum metric) const;

  double getElementMetric(libMesh::dof_id_type elem_id, int inv_index,
                          inventory_outputs::InventoryOutputsEnum metric) const;

  void initialiseNuclearInventory();

  ElementInventory &getElementInventory(size_t elem_id);

  const ElementInventory &getElementInventory(size_t elem_id) const;

protected:
  /// Mapping from nuclide names to a vector containing the quantities to be
  /// stored for that nuclide

  std::unordered_map<std::string, size_t> _nuclide_ids;

  std::unordered_map<size_t, std::string> _nuclide_names;

  std::vector<ElementInventory> _element_inventories;

  size_t _n_fispact_inventories;

  size_t _max_nuclide_id = 0;

public:
  class NuclideInventory {
  public:
    NuclideInventory(
        const std::string &nuclide_name, size_t n_inventories,
        const std::vector<nuclide_quantities::NuclideQuantitiesEnum>
            &quantities)
        : _nuclide_name(nuclide_name), _n_inventories(n_inventories) {
      for (auto &quantity : quantities) {
        _metric_data[quantity] = std::vector<double>(n_inventories, 0.0);
      }

      _n_metrics = _metric_data.size();
    };

    void addQuantity(nuclide_quantities::NuclideQuantitiesEnum quantity) {
      // _quantities.push_back(quantity);
      // _data.resize(_quantities.size() * _n_inventories);
      _metric_data[quantity] = std::vector<double>(_n_inventories, 0.0);
    }

    std::vector<nuclide_quantities::NuclideQuantitiesEnum> getQuantities() {
      std::vector<nuclide_quantities::NuclideQuantitiesEnum> quantities;

      for (auto &[quantity, value] : _metric_data) {
        quantities.push_back(quantity);
      }
      return quantities;
    }

    double &getQuantity(nuclide_quantities::NuclideQuantitiesEnum quantity,
                        size_t inv_index) {

      return _metric_data.at(quantity)[inv_index];

      // return _data[(_quantities.size() * inv_index) +
      // findQuantity(quantity)];
    }

    const double &
    getQuantity(nuclide_quantities::NuclideQuantitiesEnum quantity,
                size_t inv_index) const {

      return _metric_data.at(quantity)[inv_index];
      // return _data[(_quantities.size() * inv_index) +
      // findQuantity(quantity)];
    }

    void setQuantity(nuclide_quantities::NuclideQuantitiesEnum quantity,
                     size_t inv_index, const double &value) {

      _metric_data.at(quantity)[inv_index] = value;
      // _data[(_quantities.size() * inv_index) + findQuantity(quantity)] =
      // value;
    }

    size_t
    findQuantity(nuclide_quantities::NuclideQuantitiesEnum quantity) const {
      for (int i = 0; i < _quantities.size(); i++) {
        if (quantity == _quantities[i]) {
          return i;
        }
      }
      throw std::runtime_error(
          "Searching for nuclide quantity that is not stored.");
    }

    // double &operator()(size_t inv_index,
    //                    nuclide_quantities::NuclideQuantitiesEnum quantity) {
    //   size_t quantity_id = findQuantity(quantity);
    //   return _data[(inv_index * _n_metrics) + quantity_id];
    // }

  protected:
    std::string _nuclide_name;
    size_t _n_inventories;
    size_t _n_metrics;

    std::vector<nuclide_quantities::NuclideQuantitiesEnum> _quantities;
    std::vector<double> _data;

    std::unordered_map<nuclide_quantities::NuclideQuantitiesEnum,
                       std::vector<double>>
        _metric_data;
  };

  ///
  class ElementInventory {
  public:
    ElementInventory(size_t elem_id) : _elem_id(elem_id) {}

    void registerNuclideMetricRequest(
        const std::string &name, const size_t id, const size_t n_inventories,
        const nuclide_quantities::NuclideQuantitiesEnum metric) {

      if (std::find(_nuclides.begin(), _nuclides.end(), id) ==
          _nuclides.end()) {

        _nuclides.push_back(id);
        _nuclide_data.emplace_back(
            name, n_inventories,
            std::vector<nuclide_quantities::NuclideQuantitiesEnum>(metric));
        _nuclide_data.back().addQuantity(metric);
      } else {
        getNuclide(id).addQuantity(metric);
      }
    }

    void registerElementMetricRequest(
        const size_t n_inventories,
        const inventory_outputs::InventoryOutputsEnum metric) {

      _inventory_metrics[metric] = std::vector<double>(n_inventories, 0.0);
    }

    bool has_nuclide(size_t nuclide_id) const {
      if (find_nuclide(nuclide_id == -1)) {
        return false;
      }
      return true;
    }

    int find_nuclide(size_t nuclide_id) const {
      for (size_t i = 0; i < _nuclides.size(); ++i)
        if (_nuclides[i] == nuclide_id)
          return i;
      return -1;
    }

    NuclideInventory &getNuclide(size_t n) {
      return _nuclide_data[find_nuclide(n)];
    }

    const NuclideInventory &getNuclide(size_t n) const {
      return _nuclide_data[find_nuclide(n)];
    }

    size_t getElemID() { return _elem_id; }

    std::vector<size_t> &getNuclides() { return _nuclides; }

    std::unordered_map<inventory_outputs::InventoryOutputsEnum,
                       std::vector<double>> &
    getInventoryMetricMap() {
      return _inventory_metrics;
    }

    const double
    getElementMetric(inventory_outputs::InventoryOutputsEnum metric,
                     int inv_index) const {
      return _inventory_metrics.at(metric)[inv_index];
    }

  protected:
    std::vector<NuclideInventory> _nuclide_data;
    std::vector<size_t> _nuclides;
    libMesh::dof_id_type _elem_id;

    size_t _n_element_metrics;

    std::unordered_map<inventory_outputs::InventoryOutputsEnum,
                       std::vector<double>>
        _inventory_metrics;
  };
};
