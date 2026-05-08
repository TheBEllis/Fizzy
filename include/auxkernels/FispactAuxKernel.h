#include "AuxKernel.h"
#include "FispactInventoryManager.h"
#include "FizzyEnums.h"

class FispactAuxKernel : public AuxKernel {

public:
  static InputParameters validParams() {
    InputParameters params = AuxKernel::validParams();
    return params;
  }

  FispactAuxKernel(const InputParameters &params) : AuxKernel(params) {
    if (mooseVariableBase()->feType() != libMesh::FEType(CONSTANT, MONOMIAL))
      paramError("variable", "Must be of type CONSTANT MONOMIAL");
  }

protected:
  virtual Real computeValue() override = 0;

  const FispactInventoryManager &getInventoryManager() {
    return getUserObject<FispactInventoryManager>("InvManager");
  }
};
