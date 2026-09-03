#include "FispactFluxInput.h"
#include "InputParameters.h"
#include "libMeshReducedNamespace.h"
#include "libmesh/id_types.h"

class OpenMCFluxInput : public FispactFluxInput {

public:
  static InputParameters validParams();

  OpenMCFluxInput(const InputParameters &params);

  virtual std::vector<double> getElemFlux(dof_id_type elem_id);

  void readEnergyFilterBins();
  void readNumEnergyBins();

protected:
  size_t _flux_tally_id;

  size_t _energy_filter_id;

  std::string _statepoint_filename;
};
