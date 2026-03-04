#pragma once

#include "BlockRestrictable.h"
#include "GeneralUserObject.h"
#include <unordered_map>

class FispactMaterial : public GeneralUserObject, public BlockRestrictable {
public:
  static InputParameters validParams();

  FispactMaterial(const InputParameters &params);

  virtual void initialize() {}
  virtual void finalize() {}
  virtual void execute() {}

  /**
   * Return a vector of Z or Zai numbers, depending on whether we are using
   * setMass or setFuel, corresponding to the nuclides in _nuclides.
   */
  // std::vector<int> getZsorZais(const fispact::FispactMonitor &monitor) const;

  /*
   * Getter for _nuclides
   */
  const std::vector<std::string> getNuclides() const;

  /*
   * Getter for _nuclide_fraction
   */
  const std::vector<double> getNuclideFractions() const;

  const std::unordered_map<std::string, double> &getNuclideFractionMap() const {
    return _nuclide_fraction_map;
  }
  /*
   * Getter for _nuclide_fraction
   */
  const double &getDensity() const { return _density; }

  /*
   * Getter for _material_type
   */
  const std::string &getMaterialType() const { return _material_type; }

protected:
  //
  std::unordered_map<std::string, double> _nuclide_fraction_map;

  // Density of the fispact material
  double _density;

  //
  const std::string _material_type;
};
