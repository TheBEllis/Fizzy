#pragma once
#include "BlockRestrictable.h"
#include "FispactPostprocessor.h"

class SubdomainPhotonEmission : public FispactPostprocessor,
                                public BlockRestrictable {

public:
  static InputParameters validParams();

  SubdomainPhotonEmission(const InputParameters &params);

  virtual void initialize();
  virtual void execute();
  virtual void finalize();

  virtual PostprocessorValue getValue() const;

  Real _sum;
};
