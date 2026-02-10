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

  // Times are known for all processes already
  params.set<bool>("auto_broadcast") = false;
  params.set<bool>("dynamic_time_sequence") = false;

  return params;
}

FispactScheduleTimes::FispactScheduleTimes(const InputParameters &parameters)
    : Times(parameters) {

  const FispactSchedule &schedule =
      getUserObject<FispactSchedule>(getParam<std::string>("FispactSchedule"));
  _times = schedule.getTimes();
}
