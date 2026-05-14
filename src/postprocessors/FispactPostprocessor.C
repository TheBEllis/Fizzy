#include "FispactPostprocessor.h"

FispactPostprocessor::FispactPostprocessor(const InputParameters &params)
    : GeneralPostprocessor(params),
      _fispact_problem(static_cast<FispactProblem &>(_fe_problem)) {}

size_t FispactPostprocessor::getFispactInventoryIdx() {

  const FispactSchedule &schedule = getFispactProblem().getSchedule();

  const std::vector<double> &schedule_times = schedule.getCumulativeTimes();

  double inventory_time;
  if (!_is_transient) {
    inventory_time = getFispactProblem().getOutputInventoryTime();

  } else {
    inventory_time = _t;
  }

  auto schedule_iterator =
      std::find(schedule_times.begin(), schedule_times.end(), inventory_time);

  if (schedule_iterator == schedule_times.end()) {
    mooseError("Current time " + std::to_string(inventory_time) +
               " does not match any entry in the FISPACT schedule. Cannot "
               "do interprocess communication");
  }
  size_t inventory_idx =
      std::distance(schedule_times.begin(), schedule_iterator);

  return inventory_idx;
}
