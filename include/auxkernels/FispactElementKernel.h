#include "FispactAuxKernel.h"
#include "FizzyEnums.h"
#include "InputParameters.h"

class FispactElementKernel : public FispactAuxKernel {

public:
  static InputParameters validParams();

  FispactElementKernel(const InputParameters &params);

protected:
  virtual Real computeValue() override;

  std::string _nuclide;

  inventory_outputs::InventoryOutputsEnum _metric;
};
