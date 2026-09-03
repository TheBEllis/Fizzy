#include "FispactPostprocessor.h"
#include "FispactProblem.h"

FispactPostprocessor::FispactPostprocessor(const InputParameters &params)
    : GeneralPostprocessor(params),
      _fispact_problem(static_cast<FispactProblem &>(_fe_problem)) {}

FispactProblem &FispactPostprocessor::getFispactProblem() {
  return _fispact_problem;
}
