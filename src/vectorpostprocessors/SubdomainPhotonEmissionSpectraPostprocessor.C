#include "SubdomainPhotonEmissionSpectraPostprocessor.h"

registerMooseObject("FizzyApp", SubdomainPhotonEmissionSpectraPostprocessor);

InputParameters SubdomainPhotonEmissionSpectraPostprocessor::validParams() {

  InputParameters params = FispactVectorPostprocessor::validParams();
  params += BlockRestrictable::validParams();

  return params;
}

SubdomainPhotonEmissionSpectraPostprocessor::
    SubdomainPhotonEmissionSpectraPostprocessor(const InputParameters &params)
    : FispactVectorPostprocessor(params), BlockRestrictable(this),
      _photon_emission_spectra(declareVector("photon_emission_spectra")) {

  _photon_emission_spectra.resize(24);
}

void SubdomainPhotonEmissionSpectraPostprocessor::execute() {
  PhotonSpectra *spectra = getFispactProblem().getPhotonSpectra();

  std::unordered_map<uint64_t, uint64_t> &local_element_index =
      getFispactProblem().getLocalElemIndexMap();

  size_t inventory_index = getFispactInventoryIdx();
  for (auto &block_id : blockIDs()) {
    for (libMesh::Elem *elem :
         getFispactProblem()
             .mesh()
             .getMesh()
             .active_subdomain_elements_ptr_range(block_id)) {

      if (elem->processor_id() == processor_id()) {
        std::vector<double>::iterator begin = spectra->spectrum_begin(
            inventory_index, local_element_index[elem->id()]);

        std::vector<double>::iterator end = spectra->spectrum_end(
            inventory_index, local_element_index[elem->id()]);

        std::transform(begin, end, _photon_emission_spectra.begin(),
                       _photon_emission_spectra.begin(), std::plus<double>());
      }
    }
  }
}

void SubdomainPhotonEmissionSpectraPostprocessor::finalize() {

  _communicator.sum(_photon_emission_spectra);
}
