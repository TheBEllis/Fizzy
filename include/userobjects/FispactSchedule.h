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

protected:
  /// Vector of flux amplitudes
  const std::vector<double> _flux_schedule;

  /// List of times corresponding to the list of flux amplitudes on or off
  const std::vector<double> _times;
};
