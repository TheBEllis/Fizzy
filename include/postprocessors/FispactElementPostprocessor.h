#pragma once
#include "ElementPostprocessor.h"
#include "InputParameters.h"
#include "UserObject.h"

// Forward dec
class FispactProblem;

class FispactElementPostprocessor : public ElementPostprocessor {

public:
  static InputParameters validParams() {
    InputParameters params = ElementPostprocessor::validParams();
    return params;
  };

  FispactElementPostprocessor(const InputParameters &params);

  virtual void initialize() = 0;
  virtual void execute() = 0;
  virtual void finalize() = 0;
  virtual void threadJoin(const UserObject &y) = 0;

  virtual PostprocessorValue getValue() const = 0;

  FispactProblem &getFispactProblem();

protected:
  FispactProblem &_fispact_problem;
};
