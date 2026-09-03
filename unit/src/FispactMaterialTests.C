#include "FispactMaterialTest.h"
#include "gtest/gtest.h"

TEST_F(FispactMaterialTests, get_average_molar_mass) {

  std::vector<std::string> nuclide_names = {"H1", "He3"};
  std::vector<double> nuclide_fractions = {0.5, 0.5};
  std::string fraction_type = "wo";

  EXPECT_DOUBLE_EQ(_mat1->getAverageMolarMass(nuclide_names, nuclide_fractions,
                                              fraction_type),
                   1.5);

  fraction_type = "ao";

  EXPECT_DOUBLE_EQ(_mat1->getAverageMolarMass(nuclide_names, nuclide_fractions,
                                              fraction_type),
                   2);

  nuclide_names = {"H1", "H2", "He3"};
  nuclide_fractions = {0.1, 0.4, 0.6};
  fraction_type = "wo";

  EXPECT_DOUBLE_EQ(_mat1->getAverageMolarMass(nuclide_names, nuclide_fractions,
                                              fraction_type),
                   2);

  fraction_type = "ao";

  EXPECT_DOUBLE_EQ(_mat1->getAverageMolarMass(nuclide_names, nuclide_fractions,
                                              fraction_type),
                   2.7);
}

TEST_F(FispactMaterialTests, convertAtomToMassFraction) {

  std::vector<std::string> nuclide_names = {"H1", "He3"};
  std::vector<double> nuclide_fractions = {0.5, 0.5};
  std::string fraction_type = "ao";

  double average_molar_mass = _mat1->getAverageMolarMass(
      nuclide_names, nuclide_fractions, fraction_type);

  _mat1->convertFromAtomToMassFraction(nuclide_names, nuclide_fractions,
                                       average_molar_mass);

  std::vector<double> expected_mass_frac = {0.25, 0.75};

  EXPECT_EQ(nuclide_fractions, expected_mass_frac);
}
