#pragma once
#include "FispactProblem.h"
#include "GeneralPostprocessor.h"
#include "InputParameters.h"

class FispactPostprocessor : public GeneralPostprocessor {

public:
  static InputParameters validParams() {
    InputParameters params = GeneralPostprocessor::validParams();
    return params;
  };

  FispactPostprocessor(const InputParameters &params);

  virtual void initialize() = 0;
  virtual void execute() = 0;
  virtual void finalize() = 0;

  virtual PostprocessorValue getValue() const = 0;

  FispactProblem &getFispactProblem() { return _fispact_problem; }

  size_t getFispactInventoryIdx();

  FispactProblem &_fispact_problem;
};
