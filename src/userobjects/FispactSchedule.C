#include "FispactSchedule.h"
#include "MooseError.h"

registerMooseObject("FizzyApp", FispactSchedule);

InputParameters FispactSchedule::validParams() {
  InputParameters params = GeneralUserObject::validParams();

  params.addRequiredParam<std::vector<double>>("times",
                                               "a list of times with each ");

  params.addRequiredParam<std::vector<double>>(
      "flux_amplitude", "Vector of flux amplitudes for each time.");

  return params;
}

FispactSchedule::FispactSchedule(const InputParameters &parameters)
    : GeneralUserObject(parameters),
      _flux_amplitude(getParam<std::vector<double>>("flux_amplitude")),
      _times(getParam<std::vector<double>>("times")),
      _n_solution_inventories(_times.size()),
      _n_inventories(_n_solution_inventories + 1) {
  // Check vectors passed in are the right size
  if (_times.size() != _flux_amplitude.size()) {
    mooseError("The length of the flux_schdule and times must be equal");
  }

  _cumulative_times.resize(_times.size());
  std::inclusive_scan(_times.begin(), _times.end(), _cumulative_times.begin());
}

const std::vector<double> &FispactSchedule::getFluxAmplitude() const {
  return _flux_amplitude;
}

std::vector<double> FispactSchedule::getFluxAmplitude() {
  return _flux_amplitude;
}

const std::vector<double> &FispactSchedule::getTimes() const { return _times; }

const std::vector<double> &FispactSchedule::getCumulativeTimes() const {
  return _cumulative_times;
}

const size_t &FispactSchedule::getNumInventories() const {
  return _n_inventories;
}

const size_t &FispactSchedule::getNumSolutionInventories() const {
  return _n_solution_inventories;
}
