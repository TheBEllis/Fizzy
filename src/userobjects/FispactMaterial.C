#include "FispactFactory.h"
#include "FispactMaterial.h"
#include "FispactProblem.h"

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
      "nuclide_fraction",
      "List of decimal fractions dictating the proportions of "
      "the nuclides in the material.");

  params.addRequiredParam<MooseEnum>(
      "fraction_type", MooseEnum("wo ao"),
      "Determines whether supplied material fractions are mass fractions "
      "or atomic fractions.");

  params.addRequiredParam<double>("density", "Total material density");

  params.addRequiredParam<MooseEnum>(
      "density_units", MooseEnum("g/cm3 kg/m3 atom/b-cm atom/cm3"),
      "Material density units");

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
      _fispact_problem(static_cast<FispactProblem &>(_fe_problem)),
      _molar_mass_map(getFispactProblem().getMolarMassMap()),
      _density(getParam<double>("density")),
      _material_type(getParam<MooseEnum>("material_type")),
      _avogadro(getFispactProblem().avogadroNumber()) {

  std::vector<std::string> nuclide_names =
      getParam<std::vector<std::string>>("nuclides");
  std::vector<double> nuclide_fractions =
      getParam<std::vector<double>>("nuclide_fraction");

  // Check _nuclides and _nuclide_fractions are of equal length
  if (nuclide_names.size() != nuclide_fractions.size()) {
    mooseError("Nuclides list length does not match nuclide fraction list "
               "length in material " +
               name());
  }

  // Check that nuclide fractions sum to equal to or less than 1
  if (!(std::accumulate(nuclide_fractions.begin(), nuclide_fractions.end(),
                        0) <= 1)) {
    mooseError("FispactMaterial " + name() +
               " has nuclide fractions summing to more than 1.");
  }

  // Check that the user has either: NOT set isotope numbers if the material
  // type is set to MASS OR: HAS provided isotope number if material type is set
  // to FUEL
  for (auto &nuclide : nuclide_names) {
    if (std::any_of(nuclide.begin(), nuclide.end(),
                    [](unsigned char c) { return std::isdigit(c); })) {

      if (_material_type == "MASS") {
        mooseError("Material type is set to MASS, but the provided nuclides "
                   "have isotope numbers specified. Please provide only "
                   "nuclide symbols, not isotope numbers.");
      }
    } else if (_material_type == "FUEL") {
      mooseError("Material set to FUEL, but provided nuclides do not have "
                 "isotope numbers defined. Please define nuclide isotopes");
    }
  }

  std::string fraction_type = getParam<MooseEnum>("fraction_type");
  std::string density_units = getParam<MooseEnum>("density_units");

  if ((_material_type == "MASS") && (fraction_type != "wo")) {

    mooseError("Fispact MASS material only supports mass fractions!");
  }

  if ((_material_type == "MASS") &&
      !((density_units == "g/cm3") || (density_units == "kg/cm3"))) {
    mooseError(
        "Fispact MASS material only supports density units of g/cm3 or kg/cm3");
  }

  /// If using a FUEL material, calculate mean molar mass for the material as
  /// we'll probably need it
  double average_molar_mass = 0;
  if (_material_type == "FUEL") {

    average_molar_mass =
        getAverageMolarMass(nuclide_names, nuclide_fractions, fraction_type);

    if (fraction_type == "ao") {
      convertFromAtomToMassFraction(nuclide_names, nuclide_fractions,
                                    average_molar_mass);
    }
  }

  /// If density is not given in g/cm3, we have some conversions to do
  if (density_units == "g/cm3") {

  } else if (density_units == "kg/m3") {
    _density /= 1000;

  } else if (density_units == "atom/cm3") {
    _density *= average_molar_mass / _avogadro;
  }

  else if (density_units == "atom/b-cm") {
    _density *= average_molar_mass * 1e24 / _avogadro;
  }

  // Actually initialise map now we've checked all the data makes sense
  for (int i = 0; i < nuclide_names.size(); i++) {
    auto [it, inserted] = _nuclide_fraction_map.try_emplace(
        nuclide_names[i], nuclide_fractions[i]);

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

double FispactMaterial::getAverageMolarMass(
    const std::vector<std::string> &nuclide_names,
    const std::vector<double> &nuclide_fractions,
    const std::string &fraction_type) const {
  double average_molar_mass = 0;
  if (fraction_type == "wo") {
    double sum_mass_frac_over_molar_mass = 0;
    for (int i = 0; i < nuclide_names.size(); i++) {

      sum_mass_frac_over_molar_mass +=
          nuclide_fractions[i] / _molar_mass_map.at(nuclide_names[i]);
    }
    average_molar_mass = 1 / sum_mass_frac_over_molar_mass;

  } else if (fraction_type == "ao") {
    for (int i = 0; i < nuclide_names.size(); i++) {

      average_molar_mass +=
          nuclide_fractions[i] * _molar_mass_map.at(nuclide_names[i]);
    }
  }
  return average_molar_mass;
}

void FispactMaterial::convertFromAtomToMassFraction(
    const std::vector<std::string> &nuclide_names,
    std::vector<double> &nuclide_fractions, const double average_molar_mass) {

  /// While we're here, now that we have the average molar mass, convert to
  /// mass fraction
  for (int i = 0; i < nuclide_names.size(); i++) {

    nuclide_fractions[i] =
        (nuclide_fractions[i] * _molar_mass_map.at(nuclide_names[i])) /
        average_molar_mass;
  }
}

FispactProblem &FispactMaterial::getFispactProblem() const {
  return _fispact_problem;
}
