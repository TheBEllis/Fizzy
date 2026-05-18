#pragma once
#include "FispactElementPostprocessor.h"
#include "FizzyEnums.h"
#include "UserObject.h"

class NuclideMetric : public FispactElementPostprocessor {

public:
  static InputParameters validParams();

  NuclideMetric(const InputParameters &params);

  virtual void initialize();
  virtual void execute();
  virtual void finalize();

  virtual void threadJoin(const UserObject &y);

  virtual PostprocessorValue getValue() const;

  Real _sum;

  std::vector<std::string> _nuclides;

  nuclide_quantities::NuclideQuantitiesEnum _metric;
};
