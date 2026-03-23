#pragma once
#include "Action.h"
#include "InputParameters.h"

class AddPhotonSpectraVectorPPAllBlocks : public Action {
public:
  static InputParameters validParams();

  AddPhotonSpectraVectorPPAllBlocks(const InputParameters &params);

  virtual void act() override;
};
