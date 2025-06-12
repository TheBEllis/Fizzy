#include "FispactSchedule.h"

registerMooseObject("fizzyApp", FispactSchedule);

InputParameters FispactSchedule::validParams() {
  InputParameters params = GeneralUserObject::validParams();

  params.addRequiredParam<std::vector<double>>("times",
                                               "a list of times with each ");

  params.addRequiredParam<std::vector<int>>(
      "schedule",
      "List of 1s and 0s of length one less than the list of times. \
                               List represents whether heating or cooling should occur.");

  return params;
}

FispactSchedule::FispactSchedule(const InputParameters &parameters)
    : GeneralUserObject(parameters),
      _schedule(getParam<std::vector<int>>("schedule")),
      _times(getParam<std::vector<double>>("times")) {
  // Check vectors passed in are the right size
  mooseAssert(
      _times.size() != _schedule.size() + 1,
      "The length of schedule must be one larger than the list of times");
}
