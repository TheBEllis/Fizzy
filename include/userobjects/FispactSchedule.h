#pragma once

#include "GeneralUserObject.h"

class FispactSchedule : public GeneralUserObject {
public:
  static InputParameters validParams();

  FispactSchedule(const InputParameters &params);

  virtual void initialize() {}
  virtual void finalize() {}
  virtual void execute() {}

  const std::vector<double> &getFluxAmplitude() const;
  std::vector<double> getFluxAmplitude();
  const std::vector<double> &getTimes() const;
  const std::vector<double> &getCumulativeTimes() const;
  const size_t &getNumInventories() const;
  const size_t &getNumSolutionInventories() const;

protected:
  /// Vector of flux amplitudes
  const std::vector<double> _flux_amplitude;

  /// List of times corresponding to the list of flux amplitudes on or off
  const std::vector<double> _times;

  /// Cumulative version of _times
  std::vector<double> _cumulative_times;

  /// Total number of solution inventories (ignoring initial)
  const size_t _n_solution_inventories;

  /// Total number of FISPACT inventories including initial inventory
  const size_t _n_inventories;
};
