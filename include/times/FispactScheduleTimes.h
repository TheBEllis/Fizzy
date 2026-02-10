#pragma once
// Moose includes
#include "Times.h"

/**
 * Simple times from an input parameter
 */
class FispactScheduleTimes : public Times {
public:
  static InputParameters validParams();
  FispactScheduleTimes(const InputParameters &parameters);
  virtual ~FispactScheduleTimes() = default;

protected:
  virtual void initialize() override {}
};
