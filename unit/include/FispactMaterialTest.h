#include "FizzyObjectUnitTest.h"

class FispactMaterialTests : public FizzyObjectUnitTest {
public:
  FispactMaterialTests() : FizzyObjectUnitTest("FizzyApp") { buildObjects(); }

protected:
  void buildObjects() {

    InputParameters pars_mat1 = _factory.getValidParams("FispactMaterial");

    pars_mat1.set<MooseEnum>("material_type") = "FUEL";
    pars_mat1.set<std::vector<std::string>>("nuclides") = {"H1"};
    pars_mat1.set<std::vector<double>>("nuclide_fraction") = {1};
    pars_mat1.set<double>("density") = 1;
    pars_mat1.set<MooseEnum>("density_units") = "g/cm3";
    pars_mat1.set<MooseEnum>("fraction_type") = "wo";
    pars_mat1.set<std::vector<SubdomainName>>("block") = {"1"};

    _fe_problem->addObject<FispactMaterial>("FispactMaterial", "test_mat_1",
                                            pars_mat1);

    _mat1 = &_fe_problem->getUserObject<FispactMaterial>("test_mat_1");
  }

  FispactMaterial *_mat1;
  // const FispactFluxInput *;
  // const FispactNuclearDataPaths *;
};
