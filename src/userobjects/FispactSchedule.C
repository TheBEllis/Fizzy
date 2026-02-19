#include "FispactSchedule.h"

registerMooseObject("FizzyApp", FispactSchedule);

InputParameters FispactSchedule::validParams() {
  InputParameters params = GeneralUserObject::validParams();

  params.addRequiredParam<std::vector<double>>("times",
                                               "a list of times with each ");

  params.addRequiredParam<std::vector<double>>(
      "flux_schedule",
      "List of 1s and 0s of length one less than the list of times. \
                               List represents whether heating or cooling should occur.");

  return params;
}

FispactSchedule::FispactSchedule(const InputParameters &parameters)
    : GeneralUserObject(parameters),
      _flux_schedule(getParam<std::vector<double>>("flux_schedule")),
      _times(getParam<std::vector<double>>("times")) {
  // Check vectors passed in are the right size
  mooseAssert(_times.size() == _flux_schedule.size(),
              "The length of the flux_schdule and times must be equal");

  _cumulative_times.resize(_times.size());
  std::inclusive_scan(_times.begin(), _times.end(), _cumulative_times.begin());
}

const std::vector<double> &FispactSchedule::getFluxSchedule() const {
  return _flux_schedule;
}

std::vector<double> FispactSchedule::getFluxSchedule() {
  return _flux_schedule;
}

const std::vector<double> &FispactSchedule::getTimes() const { return _times; }

const std::vector<double> &FispactSchedule::getCumulativeTimes() const {
  return _cumulative_times;
}
