#pragma once

#include "GeneralUserObject.h"
#include <unordered_map>

class FispactProblem;

class FispactUserObject : public GeneralUserObject {
public:
  static InputParameters validParams() {
    InputParameters params = GeneralUserObject::validParams();
    return params;
  };

  FispactUserObject(const InputParameters &params);

  virtual void initialize() {};
  virtual void finalize() {};
  virtual void execute() {};

  FispactProblem &getFispactProblem() const { return _fispact_problem; }

protected:
  FispactProblem &_fispact_problem;
};
