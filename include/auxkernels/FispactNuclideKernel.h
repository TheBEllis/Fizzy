#include "FispactAuxKernel.h"
#include "InputParameters.h"

class FispactNuclideKernel : public FispactAuxKernel {

public:
  static InputParameters validParams();

  FispactNuclideKernel(const InputParameters &params);

protected:
  virtual Real computeValue() override;

  std::string _nuclide;

  nuclide_quantities::NuclideQuantitiesEnum _metric;
};
