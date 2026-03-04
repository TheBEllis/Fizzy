#include "FispactMaterial.h"

#include <algorithm>
#include <numeric>
#include <stdexcept>

registerMooseObject("FizzyApp", FispactMaterial);

InputParameters FispactMaterial::validParams() {
  InputParameters params = GeneralUserObject::validParams();
  params += BlockRestrictable::validParams();

  params.addRequiredParam<std::vector<std::string>>(
      "nuclides", "a list of nuclides present in the material");

  params.addRequiredParam<std::vector<double>>(
      "nuclide_fraction", "List of percentages dictating the proportion of "
                          "the nuclides in the material.");

  params.addRequiredParam<double>("density",
                                  "Total material density in g/cm^3");

  // Enum to determine whether we are using FISPACT::setMass or
  // FISPACT::setFUEL
  MooseEnum material_type_enum{"MASS FUEL"};
  params.addRequiredParam<MooseEnum>(
      "material_type", material_type_enum,
      "Setting to determine if this is a FISPACT Mass or Fuel material");
  return params;
}

FispactMaterial::FispactMaterial(const InputParameters &parameters)
    : GeneralUserObject(parameters), BlockRestrictable(this),
      _density(getParam<double>("density")),
      _material_type(getParam<MooseEnum>("material_type")) {

  std::vector<std::string> nuclides =
      getParam<std::vector<std::string>>("nuclides");
  std::vector<double> nuclide_fraction =
      getParam<std::vector<double>>("nuclide_fraction");

  // Check _nuclides and _nuclide_fraction are of equal length
  if (nuclides.size() != nuclide_fraction.size()) {
    mooseError("Nuclides list length does not match nuclide fraction list "
               "length in material " +
               name());
  }

  // Check that nuclide fractions sum to equal to or less than 100
  if (!(std::accumulate(nuclide_fraction.begin(), nuclide_fraction.end(), 0) <=
        100)) {
    mooseError("FispactMaterial " + name() +
               " has nuclide fractions summing to more than 100.");
  }

  // Check that the user has either: NOT set isotope numbers if the material
  // type is set to MASS OR: HAS provided isotope number if material type is set
  // to FUEL
  for (auto &nuclide : nuclides) {
    if (std::any_of(nuclide.begin(), nuclide.end(),
                    [](unsigned char c) { return std::isdigit(c); })) {

      if (_material_type == "MASS") {
        mooseError("Material type is set to MASS, but the provided nuclides "
                   "have isotope numbers specified. Please provide only "
                   "nuclide symbols, not isotope numbers.");
      }
    } else if (_material_type == "FUEL") {
      mooseError("Material to FUEL, but provided nuclides do not have "
                 "isotope numbers defined. Please define nuclide isotopes");
    }
  }

  // Actually initialise map now we've checked all the data makes sense
  for (int i = 0; i < nuclides.size(); i++) {
    auto [it, inserted] =
        _nuclide_fraction_map.try_emplace(nuclides[i], nuclide_fraction[i]);

    // If we attempt to give the same isotope two mass fraction values, freak
    // out
    if (!inserted) {
      throw std::runtime_error(
          "Duplicate element in material definition for material " + name());
    }
  }
}

const std::vector<std::string> FispactMaterial::getNuclides() const {

  std::vector<std::string> nuclides;

  for (auto &pair : _nuclide_fraction_map) {
    nuclides.push_back(pair.first);
  }
  return nuclides;
}

const std::vector<double> FispactMaterial::getNuclideFractions() const {

  std::vector<double> nuclide_fractions;

  for (auto &pair : _nuclide_fraction_map) {
    nuclide_fractions.push_back(pair.second);
  }
  return nuclide_fractions;
}
