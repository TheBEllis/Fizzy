#include "AverageElementSize.h"
#include "FizzyObjectUnitTest.h"
#include "gtest/gtest.h"
#include <stdexcept>

class FispactProblemTest : public FizzyObjectUnitTest {
public:
  FispactProblemTest() : FizzyObjectUnitTest("FizzyApp") {}
};

TEST_F(FispactProblemTest, isFlux) {

  std::vector<double> no_flux({0, 0, 0, 0});
  std::vector<double> flux({1, 0, 0, 0});
  EXPECT_TRUE(_fe_problem->isFlux(flux));
  EXPECT_FALSE(_fe_problem->isFlux(no_flux));
}

TEST_F(FispactProblemTest, convertGammaEvToCount) {
  _fe_problem->setPhotonBins({0, 2, 6, 10});
  std::vector<double> input = {1, 2, 3};
  std::vector<double> output;
  _fe_problem->convertGammaEvToCount(input, output);
  EXPECT_EQ(output, std::vector<double>({1e6, 5e5, 3.75e5}));
}

TEST_F(FispactProblemTest, loadMolarMasses) {

  const std::unordered_map<std::string, double> &molar_masses =
      _fe_problem->getMolarMassMap();

  std::vector<int> zais;

  std::vector<double> masses;

  EXPECT_EQ(molar_masses.at("H1"), 1);
  EXPECT_EQ(molar_masses.at("H2"), 2);
  EXPECT_EQ(molar_masses.at("He3"), 3);

  EXPECT_NE(molar_masses.at("H1"), 4);

  EXPECT_THROW(molar_masses.at("H33"), std::out_of_range);
}

// class FispactProblemTestUserObjects : public FizzyObjectUnitTest {
// public:
//   FispactProblemTestUserObjects() : FizzyObjectUnitTest("FizzyApp") {
//     buildObjects();
//   }
//
// protected:
//   void buildObjects() {
//
//     InputParameters pars_mat1 = _factory.getValidParams("FispactMaterial");
//
//     pars_mat1.set<MooseEnum>("material_type") = "FUEL";
//     pars_mat1.set<std::vector<std::string>>("nuclides") = {"H1"};
//     pars_mat1.set<std::vector<double>>("nuclide_fraction") = {100};
//     pars_mat1.set<double>("density") = 1;
//     pars_mat1.set<MooseEnum>("density_units") = "g/cm3";
//     pars_mat1.set<MooseEnum>("fraction_type") = "wo";
//     pars_mat1.set<std::vector<SubdomainName>>("block") = {"1"};
//     _fe_problem->addObject<FispactMaterial>("FispactMaterial", "Hydrogen",
//                                             pars_mat1);
//     _mat1 = &_fe_problem->getUserObject<FispactMaterial>("Hydrogen");
//   }
//
//   const FispactMaterial *_mat1;
//   const FispactSchedule *_schedule;
//   // const FispactFluxInput *;
//   // const FispactNuclearDataPaths *;
// };

/// Test fixture to make sure getElementMaterial throws if two materials are
/// defined on the same block
class FispactProblemTestThrowMaterials : public FizzyObjectUnitTest {
public:
  FispactProblemTestThrowMaterials() : FizzyObjectUnitTest("FizzyApp") {
    buildObjects();
  }

protected:
  void buildObjects() {

    InputParameters pars_mat1 = _factory.getValidParams("FispactMaterial");
    InputParameters pars_mat2 = _factory.getValidParams("FispactMaterial");

    pars_mat1.set<MooseEnum>("material_type") = "FUEL";
    pars_mat1.set<std::vector<std::string>>("nuclides") = {"H1"};
    pars_mat1.set<std::vector<double>>("nuclide_fraction") = {1};
    pars_mat1.set<double>("density") = 1;
    pars_mat1.set<MooseEnum>("density_units") = "g/cm3";
    pars_mat1.set<MooseEnum>("fraction_type") = "wo";
    pars_mat1.set<std::vector<SubdomainName>>("block") = {"1"};
    _fe_problem->addObject<FispactMaterial>("FispactMaterial", "Hydrogen",
                                            pars_mat1);
    _mat1 = &_fe_problem->getUserObject<FispactMaterial>("Hydrogen");

    pars_mat2.set<MooseEnum>("material_type") = "FUEL";
    pars_mat2.set<std::vector<std::string>>("nuclides") = {"He3"};
    pars_mat2.set<std::vector<double>>("nuclide_fraction") = {1};
    pars_mat2.set<double>("density") = 3;
    pars_mat2.set<MooseEnum>("density_units") = "g/cm3";
    pars_mat2.set<MooseEnum>("fraction_type") = "wo";
    pars_mat2.set<std::vector<SubdomainName>>("block") = {"1"};
    _fe_problem->addUserObject("FispactMaterial", "Helium", pars_mat2);
    _mat2 = &_fe_problem->getUserObject<FispactMaterial>("Helium");
  }

  const FispactMaterial *_mat1;
  const FispactMaterial *_mat2;
};

TEST_F(FispactProblemTestThrowMaterials, getElementMaterial) {
  dof_id_type elem_id = _fe_problem->mesh().elemPtr(0)->id();
  EXPECT_THROW(_fe_problem->getElementMaterial(elem_id), std::runtime_error);
}
