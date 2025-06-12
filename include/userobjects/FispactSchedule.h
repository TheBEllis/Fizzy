#pragma once

#include "GeneralUserObject.h"

class FispactSchedule : public GeneralUserObject {
public:
  static InputParameters validParams();

  FispactSchedule(const InputParameters &params);

  virtual void initialize() {}
  virtual void finalize() {}
  virtual void execute() {}

protected:
  /// Vector of flux amplitudes
  const std::vector<int32_t> &_schedule;

  /// List of times corresponding to the list of flux amplitudes on or off
  const std::vector<double> &_times;
};
