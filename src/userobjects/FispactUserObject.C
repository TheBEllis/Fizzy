#include "FispactProblem.h"
#include "FispactUserObject.h"

FispactUserObject::FispactUserObject(const InputParameters &params)
    : GeneralUserObject(params),
      _fispact_problem(static_cast<FispactProblem &>(_fe_problem)) {}
