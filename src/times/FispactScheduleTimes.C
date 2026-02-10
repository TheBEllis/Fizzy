#include "FispactSchedule.h"
#include "FispactScheduleTimes.h"

registerMooseObject("MooseApp", FispactScheduleTimes);

InputParameters FispactScheduleTimes::validParams() {
  InputParameters params = Times::validParams();
  params.addClassDescription(
      "Times set from a FispactSchedule UserObject defined in the input file");

  params.addRequiredParam<std::vector<Real>>(
      "times", "Times to store in the times vector");

  params.addRequiredParam<std::string>(
      "FispactSchedule",
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
    : Times(parameters) {

  const FispactSchedule &schedule =
      getUserObject<FispactSchedule>(getParam<std::string>("FispactSchedule"));

  if (isParamValid("FispactScheduleTimeIndices")) {
    std::vector<int> indices =
        getParam<std::vector<int>>("FispactScheduleTimeIndices");

    std::vector<double> schedule_times = schedule.getTimes();
    for (int index : indices) {
      _times.push_back(schedule_times[index]);
    }
    _times = schedule.getTimes();
  } else {

    _times = schedule.getTimes();
  }
}
