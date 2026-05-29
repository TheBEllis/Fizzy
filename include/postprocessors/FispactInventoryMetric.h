#pragma once
#include "FispactElementPostprocessor.h"
#include "FizzyEnums.h"
#include "UserObject.h"

class FispactInventoryMetric : public FispactElementPostprocessor {

public:
  static InputParameters validParams();

  FispactInventoryMetric(const InputParameters &params);

  virtual void initialize();
  virtual void execute();
  virtual void finalize();

  virtual void threadJoin(const UserObject &y);

  virtual PostprocessorValue getValue() const;

  Real _sum;

  inventory_outputs::InventoryOutputsEnum _metric;
};
