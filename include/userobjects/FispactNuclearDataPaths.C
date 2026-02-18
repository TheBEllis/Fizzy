#include "FispactNuclearDataPaths.h"
#include "MooseError.h"

#include <algorithm>
#include <filesystem>
#include <numeric>
#include <stdexcept>
#include <string>

registerMooseObject("FizzyApp", FispactNuclearDataPaths);

InputParameters FispactNuclearDataPaths::validParams() {
  InputParameters params = GeneralUserObject::validParams();

  params.addParam<std::string>("base_path", "", "");
  params.addParam<std::string>("ND_IND_NUC_KEY", "");
  params.addParam<std::string>("ND_HAZARDS_KEY", "");
  params.addParam<std::string>("ND_ABSORP_KEY", "");
  params.addParam<std::string>("ND_CLEAR_KEY", "");
  params.addParam<std::string>("ND_A2DATA_KEY", "");
  params.addParam<std::string>("ND_ENBINS_KEY", "");
  params.addParam<std::string>("ND_DECAY_KEY", "");
  params.addParam<std::string>("ND_DK_ENDF_KEY", "");
  params.addParam<std::string>("ND_PROB_TAB_KEY", "");
  params.addParam<std::string>("ND_ASSCFY_KEY", "");
  params.addParam<std::string>("ND_FISSYLD_KEY", "");
  params.addParam<std::string>("ND_FY_ENDF_KEY", "");
  params.addParam<std::string>("ND_SF_ENDF_KEY", "");
  params.addParam<std::string>("ND_SP_ENDF_KEY", "");
  params.addParam<std::string>("ND_XS_EXTRA_KEY", "");
  params.addParam<std::string>("ND_CROSSEC_KEY", "");
  params.addParam<std::string>("ND_CROSSUNC_KEY", "");
  params.addParam<std::string>("ND_XS_ENDF_KEY", "");
  params.addParam<std::string>("ND_XS_ENDFB_KEY", "");

  return params;
}

FispactNuclearDataPaths::FispactNuclearDataPaths(
    const InputParameters &parameters)
    : GeneralUserObject(parameters) {}

void FispactNuclearDataPaths::loadNuclearData(
    fispact::NuclearData &nd, fispact::FispactMonitor &fp_monitor,
    std::function<void(std::string, std::string, int, int)> callback) {
  fispact::io::NuclearDataReader nd_reader(fp_monitor);

  std::string base_path = getParam<std::string>("base_path");

  if (isParamValid("ND_IND_NUC_KEY")) {
    /// check path exists
    std::string path = base_path + getParam<std::string>("ND_IND_NUC_KEY");

    if (!(std::filesystem::exists(path))) {
      paramError("ND_IND_NUC_KEY",
                 path + " is not a valid path! Could not set ");
    }

    nd_reader.setPath(FISPACT_ND_IND_NUC_KEY, path);
  }

  if (isParamValid("ND_HAZARDS_KEY")) {
    /// check path exists
    std::string path = base_path + getParam<std::string>("ND_HAZARDS_KEY");

    if (!(std::filesystem::exists(path))) {
      paramError("ND_HAZARDS_KEY",
                 path + " is not a valid path! Could not set ");
    }

    nd_reader.setPath(FISPACT_ND_HAZARDS_KEY, path);
  }

  if (isParamValid("ND_ABSORP_KEY")) {
    /// check path exists
    std::string path = base_path + getParam<std::string>("ND_ABSORP_KEY");

    if (!(std::filesystem::exists(path))) {
      paramError("ND_ABSORP_KEY",
                 path + " is not a valid path! Could not set ");
    }

    nd_reader.setPath(FISPACT_ND_ABSORP_KEY, path);
  }

  if (isParamValid("ND_CLEAR_KEY")) {
    /// check path exists
    std::string path = base_path + getParam<std::string>("ND_CLEAR_KEY");

    if (!(std::filesystem::exists(path))) {
      paramError("ND_CLEAR_KEY", path + " is not a valid path! Could not set ");
    }

    nd_reader.setPath(FISPACT_ND_CLEAR_KEY, path);
  }

  if (isParamValid("ND_A2DATA_KEY")) {
    /// check path exists
    std::string path = base_path + getParam<std::string>("ND_A2DATA_KEY");

    if (!(std::filesystem::exists(path))) {
      paramError("ND_A2DATA_KEY",
                 path + " is not a valid path! Could not set ");
    }

    nd_reader.setPath(FISPACT_ND_A2DATA_KEY, path);
  }

  if (isParamValid("ND_ENBINS_KEY")) {
    /// check path exists
    std::string path = base_path + getParam<std::string>("ND_ENBINS_KEY");

    if (!(std::filesystem::exists(path))) {
      paramError("ND_ENBINS_KEY",
                 path + " is not a valid path! Could not set ");
    }

    nd_reader.setPath(FISPACT_ND_ENBINS_KEY, path);
  }

  if (isParamValid("ND_DECAY_KEY")) {
    /// check path exists
    std::string path = base_path + getParam<std::string>("ND_DECAY_KEY");

    if (!(std::filesystem::exists(path))) {
      paramError("ND_DECAY_KEY", path + " is not a valid path! Could not set ");
    }

    nd_reader.setPath(FISPACT_ND_DECAY_KEY, path);
  }

  if (isParamValid("ND_DK_ENDF_KEY")) {
    /// check path exists
    std::string path = base_path + getParam<std::string>("ND_DK_ENDF_KEY");

    if (!(std::filesystem::exists(path))) {
      paramError("ND_DK_ENDF_KEY",
                 path + " is not a valid path! Could not set ");
    }

    nd_reader.setPath(FISPACT_ND_DK_ENDF_KEY, path);
  }

  if (isParamValid("ND_PROB_TAB_KEY")) {
    /// check path exists
    std::string path = base_path + getParam<std::string>("ND_PROB_TAB_KEY");

    if (!(std::filesystem::exists(path))) {
      paramError("ND_PROB_TAB_KEY",
                 path + " is not a valid path! Could not set ");
    }

    nd_reader.setPath(FISPACT_ND_PROB_TAB_KEY, path);
  }

  if (isParamValid("ND_ASSCFY_KEY")) {
    /// check path exists
    std::string path = base_path + getParam<std::string>("ND_ASSCFY_KEY");

    if (!(std::filesystem::exists(path))) {
      paramError("ND_ASSCFY_KEY",
                 path + " is not a valid path! Could not set ");
    }

    nd_reader.setPath(FISPACT_ND_ASSCFY_KEY, path);
  }

  if (isParamValid("ND_FISSYLD_KEY")) {
    /// check path exists
    std::string path = base_path + getParam<std::string>("ND_FISSYLD_KEY");

    if (!(std::filesystem::exists(path))) {
      paramError("ND_FISSYLD_KEY",
                 path + " is not a valid path! Could not set ");
    }

    nd_reader.setPath(FISPACT_ND_FISSYLD_KEY, path);
  }

  if (isParamValid("ND_FY_ENDF_KEY")) {
    /// check path exists
    std::string path = base_path + getParam<std::string>("ND_FY_ENDF_KEY");

    if (!(std::filesystem::exists(path))) {
      paramError("ND_FY_ENDF_KEY",
                 path + " is not a valid path! Could not set ");
    }

    nd_reader.setPath(FISPACT_ND_FY_ENDF_KEY, path);
  }

  if (isParamValid("ND_SF_ENDF_KEY")) {
    /// check path exists
    std::string path = base_path + getParam<std::string>("ND_SF_ENDF_KEY");

    if (!(std::filesystem::exists(path))) {
      paramError("ND_SF_ENDF_KEY",
                 path + " is not a valid path! Could not set ");
    }

    nd_reader.setPath(FISPACT_ND_SF_ENDF_KEY, path);
  }

  if (isParamValid("ND_SP_ENDF_KEY")) {
    /// check path exists
    std::string path = base_path + getParam<std::string>("ND_SP_ENDF_KEY");

    if (!(std::filesystem::exists(path))) {
      paramError("ND_SP_ENDF_KEY",
                 path + " is not a valid path! Could not set ");
    }

    nd_reader.setPath(FISPACT_ND_SP_ENDF_KEY, path);
  }

  if (isParamValid("ND_XS_EXTRA_KEY")) {
    /// check path exists
    std::string path = base_path + getParam<std::string>("ND_XS_EXTRA_KEY");

    if (!(std::filesystem::exists(path))) {
      paramError("ND_XS_EXTRA_KEY",
                 path + " is not a valid path! Could not set ");
    }

    nd_reader.setPath(FISPACT_ND_XS_EXTRA_KEY, path);
  }

  if (isParamValid("ND_CROSSEC_KEY")) {
    /// check path exists
    std::string path = base_path + getParam<std::string>("ND_CROSSEC_KEY");

    if (!(std::filesystem::exists(path))) {
      paramError("ND_CROSSEC_KEY",
                 path + " is not a valid path! Could not set ");
    }

    nd_reader.setPath(FISPACT_ND_CROSSEC_KEY, path);
  }

  if (isParamValid("ND_CROSSUNC_KEY")) {
    /// check path exists
    std::string path = base_path + getParam<std::string>("ND_CROSSUNC_KEY");

    if (!(std::filesystem::exists(path))) {
      paramError("ND_CROSSUNC_KEY",
                 path + " is not a valid path! Could not set ");
    }

    nd_reader.setPath(FISPACT_ND_CROSSUNC_KEY, path);
  }

  /// If the user has provided binary xs AND uncompressed xs, prefer the binary
  /// one
  if (isParamValid("ND_XS_ENDF_KEY")) {
    if (isParamValid("ND_XS_ENDFB_KEY")) {
      paramWarning("ND_XS_ENDF_KEY",
                   "User has provided both compressed xs and uncompressed xs, "
                   "preferring compressed option.");
    } else {
      /// check path exists
      std::string path = base_path + getParam<std::string>("ND_XS_ENDF_KEY");

      if (!(std::filesystem::exists(path))) {
        paramError("ND_XS_ENDF_KEY",
                   path + " is not a valid path! Could not set ");
      }

      nd_reader.setPath(FISPACT_ND_XS_ENDF_KEY, path);
    }
  }

  if (isParamValid("ND_XS_ENDFB_KEY")) {
    /// check path exists
    std::string path = base_path + getParam<std::string>("ND_XS_ENDFB_KEY");

    if (!(std::filesystem::exists(path))) {
      paramError("ND_XS_ENDFB_KEY",
                 path + " is not a valid path! Could not set ");
    }

    nd_reader.setUseXSBinary(true);

    nd_reader.setPath(FISPACT_ND_XS_ENDFB_KEY, path);
  }

  // Load the nuclear data
  nd_reader.load(nd, callback);
}
