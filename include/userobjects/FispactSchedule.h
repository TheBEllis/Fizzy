#pragma once

#include "GeneralUserObject.h"

class FispactSchedule : public GeneralUserObject {
public:
  static InputParameters validParams();

  FispactSchedule(const InputParameters &params);

  virtual void initialize() {}
  virtual void finalize() {}
  virtual void execute() {}

  const std::vector<double> &getFluxSchedule() const;
  const std::vector<double> &getTimes() const;
  const std::vector<double> &getCumulativeTimes() const;

protected:
  /// Vector of flux amplitudes
  const std::vector<double> _flux_schedule;

  /// List of times corresponding to the list of flux amplitudes on or off
  const std::vector<double> _times;

  /// Cumulative version of _times
  std::vector<double> _cumulative_times;
};
