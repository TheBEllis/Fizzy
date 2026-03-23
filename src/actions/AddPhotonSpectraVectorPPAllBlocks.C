#include "AddPhotonSpectraVectorPPAllBlocks.h"
#include "InputParameters.h"
#include "MooseTypes.h"
#include "Registry.h"
#include "SubdomainPhotonEmissionSpectraPostprocessor.h"

registerMooseAction("FizzyApp", AddPhotonSpectraVectorPPAllBlocks,
                    "add_vector_postprocessor");

InputParameters AddPhotonSpectraVectorPPAllBlocks::validParams() {
  InputParameters params = Action::validParams();
  return params;
}

AddPhotonSpectraVectorPPAllBlocks::AddPhotonSpectraVectorPPAllBlocks(
    const InputParameters &params)
    : Action(params) {}

void AddPhotonSpectraVectorPPAllBlocks::act() {
  for (SubdomainID block_id : _mesh->meshSubdomains()) {
    const std::string &subdomain_name = _mesh->getSubdomainName(block_id);

    InputParameters vpp_params =
        _factory.getValidParams("SubdomainPhotonEmissionSpectraPostprocessor");

    vpp_params.set<std::vector<SubdomainName>>("block") = {subdomain_name};

    std::string vpp_name = subdomain_name + "_photon_emission_vpp";

    _problem->addVectorPostprocessor(
        "SubdomainPhotonEmissionSpectraPostprocessor", vpp_name, vpp_params);
  }
}
