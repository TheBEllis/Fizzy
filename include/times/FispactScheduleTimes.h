#pragma once
// Moose includes
#include "FispactSchedule.h"
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

  const FispactSchedule &_schedule;

  std::vector<Real> _fispact_times;
};
