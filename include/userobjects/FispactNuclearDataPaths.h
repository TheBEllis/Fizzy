#pragma once

#include "GeneralUserObject.h"
#include "fispactnucleardata.hpp"

class FispactNuclearDataPaths : public GeneralUserObject {
public:
  static InputParameters validParams();

  FispactNuclearDataPaths(const InputParameters &params);

  virtual void initialize() {}
  virtual void finalize() {}
  virtual void execute() {}

  void loadNuclearData(
      fispact::NuclearData &nd, fispact::FispactMonitor &fp_monitor,
      std::function<void(std::string, std::string, int, int)> callback);

protected:
  /// Nuclear data strings
  const std::string _base_path;

  const std::string _ND_IND_NUC_KEY;
  const std::string _ND_HAZARDS_KEY;
  const std::string _ND_ABSORP_KEY;
  const std::string _ND_CLEAR_KEY;
  const std::string _ND_A2DATA_KEY;
  const std::string _ND_ENBINS_KEY;
  const std::string _ND_DECAY_KEY;
  const std::string _ND_DK_ENDF_KEY;
  const std::string _ND_PROB_TAB_KEY;
  const std::string _ND_ASSCFY_KEY;
  const std::string _ND_FISSYLD_KEY;
  const std::string _ND_FY_ENDF_KEY;
  const std::string _ND_SF_ENDF_KEY;
  const std::string _ND_SP_ENDF_KEY;
  const std::string _ND_XS_EXTRA_KEY;
  const std::string _ND_CROSSEC_KEY;
  const std::string _ND_CROSSUNC_KEY;
  const std::string _ND_XS_ENDF_KEY;
  const std::string _ND_XS_ENDFB_KEY;
};
