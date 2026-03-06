#include "FizzyObjectUnitTest.h"
#include "gtest/gtest.h"

class FispactProblemTest : public FizzyObjectUnitTest {
public:
  FispactProblemTest() : FizzyObjectUnitTest("FizzyApp") {}
};

TEST_F(FispactProblemTest, test) {

  std::vector<double> no_flux({0, 0, 0, 0});
  std::vector<double> flux({1, 0, 0, 0});
  EXPECT_TRUE(_fe_problem->isFlux(flux));
  EXPECT_FALSE(_fe_problem->isFlux(no_flux));
}

class FispactProblemTestWithUserObjects : public FizzyObjectUnitTest {
public:
  FispactProblemTestWithUserObjects() : FizzyObjectUnitTest("FizzyApp") {
    buildObjects();
  }

protected:
  void buildObjects() {

    InputParameters pars_mat1 = _factory.getValidParams("FispactMaterial");
    InputParameters pars_mat2 = _factory.getValidParams("FispactMaterial");

    std::cout << _fe_problem->mesh().nSubdomains() << std::endl;

    pars_mat1.set<MooseEnum>("material_type") = "FUEL";
    pars_mat1.set<std::vector<std::string>>("nuclides") = {"H1"};
    pars_mat1.set<std::vector<double>>("nuclide_fraction") = {100};
    pars_mat1.set<double>("density") = 1;
    _fe_problem->addObject<FispactMaterial>("FispactMaterial", "Hydrogen",
                                            pars_mat1);
    _mat1 = &_fe_problem->getUserObject<FispactMaterial>("Hydrogen");

    // pars_mat2.set<MooseEnum>("material_type") = "FUEL";
    // pars_mat2.set<std::vector<std::string>>("nuclides") = {"He3"};
    // pars_mat2.set<std::vector<double>>("nuclide_fraction") = {100};
    // pars_mat2.set<std::vector<SubdomainName>>("block") = {"2"};
    // pars_mat2.set<double>("density") = 3;
    // _fe_problem->addUserObject("FispactMaterial", "Helium", pars_mat2);
    // _mat2 = &_fe_problem->getUserObject<FispactMaterial>("Helium");
  }

  const FispactMaterial *_mat1;
  const FispactMaterial *_mat2;
};

TEST_F(FispactProblemTestWithUserObjects, getElementMaterial) {
  dof_id_type elem_id = _fe_problem->mesh().elemPtr(0)->id();
  const FispactMaterial &mat = _fe_problem->getElementMaterial(elem_id);
  EXPECT_EQ(_mat1, &mat);
}
