#include "ElementPostprocessor.h"
#include "FispactElementPostprocessor.h"
#include "FispactProblem.h"

FispactElementPostprocessor::FispactElementPostprocessor(
    const InputParameters &params)
    : ElementPostprocessor(params),
      _fispact_problem(static_cast<FispactProblem &>(_fe_problem)) {}

FispactProblem &FispactElementPostprocessor::getFispactProblem() {
  return _fispact_problem;
}
