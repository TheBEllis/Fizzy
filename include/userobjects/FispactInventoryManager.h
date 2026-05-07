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

  void registerNuclideMetrics(
      std::string &nuclide,
      std::vector<nuclide_quantities::NuclideQuantitiesEnum> metric,
      std::set<SubdomainID> &blocks);

  void extractInventoryData(FispactContextBase &fp_context, size_t elem_id,
                            int inv_index);

  double
  getNuclideMetric(libMesh::dof_id_type elem_id, std::string &nuclide,
                   int inv_index,
                   nuclide_quantities::NuclideQuantitiesEnum metric) const;

  void initialiseNuclearInventory();

  ElementInventory &getElementInventory(size_t elem_id);

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
        std::string &nuclide_name, size_t n_inventories,
        std::vector<nuclide_quantities::NuclideQuantitiesEnum> &quantities)
        : _nuclide_name(nuclide_name), _n_inventories(n_inventories),
          _quantities(quantities) {
      _n_metrics = _quantities.size();
    };

    void addQuantity(nuclide_quantities::NuclideQuantitiesEnum quantity) {
      _quantities.push_back(quantity);
      _data.resize(_quantities.size() * _n_inventories);
    }

    std::vector<nuclide_quantities::NuclideQuantitiesEnum> &getQuantities() {
      return _quantities;
    }

    double &getQuantity(nuclide_quantities::NuclideQuantitiesEnum quantity,
                        size_t inv_index) {

      return _data[(_quantities.size() * inv_index) + findQuantity(quantity)];
    }

    const double &
    getQuantity(nuclide_quantities::NuclideQuantitiesEnum quantity,
                size_t inv_index) const {

      return _data[(_quantities.size() * inv_index) + findQuantity(quantity)];
    }

    void setQuantity(nuclide_quantities::NuclideQuantitiesEnum quantity,
                     size_t inv_index, const double &value) {
      _data[(_quantities.size() * inv_index) + findQuantity(quantity)] = value;
    }

    size_t
    findQuantity(nuclide_quantities::NuclideQuantitiesEnum quantity) const {
      for (int i = 0; i < _quantities.size(); i++) {
        if (quantity == _quantities[i]) {
          return i;
        }
        throw std::runtime_error(
            "Searching for nuclide quantity that is not stored.");
      }
    }

    double &operator()(size_t inv_index,
                       nuclide_quantities::NuclideQuantitiesEnum quantity) {
      size_t quantity_id = findQuantity(quantity);
      return _data[(inv_index * _n_metrics) + quantity_id];
    }

  protected:
    std::string _nuclide_name;
    size_t _n_inventories;
    size_t _n_metrics;

    std::vector<nuclide_quantities::NuclideQuantitiesEnum> &_quantities;
    std::vector<double> _data;
  };

  ///
  class ElementInventory {
  public:
    ElementInventory(size_t elem_id) : _elem_id(elem_id) {}

    void insertNuclideMetricRequest(
        std::string &name, size_t id, size_t n_inventories,
        std::vector<nuclide_quantities::NuclideQuantitiesEnum> quantities) {

      if (std::find(_nuclides.begin(), _nuclides.end(), id) ==
          _nuclides.end()) {

        _nuclides.push_back(id);
        _nuclide_data.emplace_back(name, n_inventories, quantities);
      } else {
        for (auto &quantity : quantities) {
          getNuclide(id).addQuantity(quantity);
        }
      }
    }

    size_t find_nuclide(size_t nuclide_id) const {
      for (size_t i = 0; i < _nuclides.size(); ++i)
        if (_nuclides[i] == nuclide_id)
          return i;

      throw std::runtime_error("Nuclide not present in element");
    }

    NuclideInventory &getNuclide(size_t n) {
      return _nuclide_data[find_nuclide(n)];
    }

    size_t getElemID() { return _elem_id; }

    std::vector<size_t> &getNuclides() { return _nuclides; }

  protected:
    std::vector<NuclideInventory> _nuclide_data;
    std::vector<size_t> _nuclides;
    libMesh::dof_id_type _elem_id;
  };
};
