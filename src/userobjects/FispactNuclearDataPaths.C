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
    : GeneralUserObject(parameters),
      _base_path(getParam<std::string>("base_path")) {
  loadNuclearDataPaths();
}

std::unordered_map<std::string, std::string> &getNuclearDataPathMap();

std::unordered_map<std::string, std::string> &
FispactNuclearDataPaths::getNuclearDataPathMap() {
  return _datatype_to_path_map;
}

void FispactNuclearDataPaths::loadNuclearDataPaths() {

  std::vector<std::string> nuclear_data{
      "ND_IND_NUC_KEY",  "ND_HAZARDS_KEY", "ND_ABSORP_KEY",   "ND_CLEAR_KEY",
      "ND_A2DATA_KEY",   "ND_ENBINS_KEY",  "ND_DECAY_KEY",    "ND_DK_ENDF_KEY",
      "ND_PROB_TAB_KEY", "ND_ASSCFY_KEY",  "ND_FISSYLD_KEY",  "ND_FY_ENDF_KEY",
      "ND_SF_ENDF_KEY",  "ND_SP_ENDF_KEY", "ND_XS_EXTRA_KEY", "ND_CROSSEC_KEY",
      "ND_CROSSUNC_KEY", "ND_XS_ENDF_KEY", "ND_XS_ENDFB_KEY"};
  for (auto &data : nuclear_data) {
    if (isParamValid(data)) {
      /// check path exists
      std::string path = _base_path + getParam<std::string>(data);

      if (!(std::filesystem::exists(path))) {
        paramError(data, path + " is not a valid path! Could not set ");
      }

      if (data == ("ND_XS_ENDF_KEY")) {
        if (isParamSetByUser("ND_XS_ENDFB_KEY")) {
          paramWarning(
              "ND_XS_ENDF_KEY",
              "User has provided both compressed xs and uncompressed xs, "
              "preferring compressed option.");

          /// If the user has set binary and non-binary cross sections, skip
          /// setting the non-binary ones

          continue;
        }
      }
      _datatype_to_path_map[data] = path;
    }
  }
}
