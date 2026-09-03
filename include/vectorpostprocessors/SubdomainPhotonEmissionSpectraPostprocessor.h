#pragma once
#include "BlockRestrictable.h"
#include "FispactVectorPostprocessor.h"
#include "InputParameters.h"
#include "MooseTypes.h"

class SubdomainPhotonEmissionSpectraPostprocessor
    : public FispactVectorPostprocessor,
      public BlockRestrictable {
public:
  static InputParameters validParams();

  SubdomainPhotonEmissionSpectraPostprocessor(const InputParameters &params);

  virtual void initialize();
  virtual void execute();
  virtual void finalize();

protected:
  VectorPostprocessorValue &_photon_emission_spectra;
};
