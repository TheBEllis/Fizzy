#pragma once

#include "GeneralUserObject.h"
#include "InputParameters.h"
#include "UserObject.h"
#include "libmesh/id_types.h"

class FispactFluxInput : public GeneralUserObject {

public:
  static InputParameters validParams() {

    InputParameters params = GeneralUserObject::validParams();

    params.addRequiredParam<double>("wall_loading", "");

    return params;
  }

  FispactFluxInput(const InputParameters &parameters)
      : GeneralUserObject(parameters),
        _wall_loading(getParam<double>("wall_loading")) {}

  virtual void initialize() {}
  virtual void finalize() {}
  virtual void execute() {}

  virtual std::vector<double> getElemFlux(dof_id_type elem_id) = 0;

  virtual const std::vector<double> &getFluxEnergyGroups() {
    return _flux_energy_groups;
  };

  const size_t &getNumEnergyGroups() { return _n_bins; };

  const double &getWallLoading() { return _wall_loading; };

protected:
  std::vector<double> _flux_energy_groups;

  size_t _n_bins;

  double _wall_loading;
};
