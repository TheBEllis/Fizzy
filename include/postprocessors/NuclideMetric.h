#pragma once
#include "FispactPostprocessor.h"

class NuclideMetric : public FispactPostprocessor {

public:
  static InputParameters validParams();

  NuclideMetric(const InputParameters &params);

  virtual void initialize();
  virtual void execute();
  virtual void finalize();

  virtual PostprocessorValue getValue() const;

  std::vector<dof_id_type> _element_ids;

  Real _sum;
};
