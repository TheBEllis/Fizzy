#pragma once

#include "GeneralUserObject.h"
#include <unordered_map>

class FispactNuclearDataPaths : public GeneralUserObject {
public:
  static InputParameters validParams();

  FispactNuclearDataPaths(const InputParameters &params);

  virtual void initialize() {}
  virtual void finalize() {}
  virtual void execute() {}

  void loadNuclearDataPaths();

  std::unordered_map<std::string, std::string> &getNuclearDataPathMap();

protected:
  /// Nuclear data strings
  const std::string _base_path;
  std::unordered_map<std::string, std::string> _datatype_to_path_map;
};
