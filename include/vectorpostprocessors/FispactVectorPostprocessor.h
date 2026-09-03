#pragma once
#include "FispactProblem.h"
#include "GeneralVectorPostprocessor.h"

class FispactVectorPostprocessor : public GeneralVectorPostprocessor {

public:
  static InputParameters validParams() {
    InputParameters params = GeneralVectorPostprocessor::validParams();
    return params;
  };

  FispactVectorPostprocessor(const InputParameters &params);

  virtual void initialize() = 0;
  virtual void execute() = 0;
  virtual void finalize() = 0;

  FispactProblem &getFispactProblem() { return _fispact_problem; }

  size_t getFispactInventoryIdx();

  FispactProblem &_fispact_problem;
};
