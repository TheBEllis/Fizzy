#pragma once
#include "Action.h"
#include "InputParameters.h"

class AddInventoryManager : public Action {
public:
  static InputParameters validParams();

  AddInventoryManager(const InputParameters &params);

  virtual void act() override;
};
