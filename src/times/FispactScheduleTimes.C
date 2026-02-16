#include "FispactScheduleTimes.h"
#include "MooseError.h"
#include <numeric>

registerMooseObject("MooseApp", FispactScheduleTimes);

InputParameters FispactScheduleTimes::validParams() {
  InputParameters params = Times::validParams();

  params.addClassDescription(
      "Times set from a FispactSchedule UserObject defined in the input file");

  params.addRequiredParam<UserObjectName>(
      "FispactScheduleName",
      "Name of the FispactSchedule UserObject used to set times");

  params.addParam<std::vector<int>>(
      "FispactScheduleTimeIndices",
      "The indices of the values within the schedule times the user wishes to "
      "give to the executioner. This is useful if the user wishes to only "
      "store data and/or run multiapps at a subset of the times defined within "
      "the FispactSchedule");

  // Times are known for all processes already
  params.set<bool>("auto_broadcast") = false;
  params.set<bool>("dynamic_time_sequence") = false;

  return params;
}

FispactScheduleTimes::FispactScheduleTimes(const InputParameters &parameters)
    : Times(parameters),
      _schedule(getUserObject<FispactSchedule>("FispactScheduleName")),
      _fispact_times(_schedule.getCumulativeTimes()) {

  if (isParamValid("FispactScheduleTimeIndices")) {

    const std::vector<int> &indices =
        getParam<std::vector<int>>("FispactScheduleTimeIndices");

    std::vector<bool> mask(_fispact_times.size(), false);

    for (auto index : indices) {
      mask[index] = true;
    }

    auto mask_it = mask.begin();
    _fispact_times.erase(
        std::remove_if(_fispact_times.begin(), _fispact_times.end(),
                       [&](double const &) { return !*mask_it++; }),
        _fispact_times.end());

    for (auto &time : _fispact_times) {
      _console << time << std::endl;
    }

    _times = _fispact_times;
  } else {
    _times = _fispact_times;
  }
}
