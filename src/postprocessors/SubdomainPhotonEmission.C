#include "SubdomainPhotonEmission.h"
#include <numeric>

registerMooseObject("FizzyApp", SubdomainPhotonEmission);

InputParameters SubdomainPhotonEmission::validParams() {

  InputParameters params = FispactPostprocessor::validParams();
  params += BlockRestrictable::validParams();

  return params;
}

SubdomainPhotonEmission::SubdomainPhotonEmission(const InputParameters &params)
    : FispactPostprocessor(params), BlockRestrictable(this) {}

void SubdomainPhotonEmission::initialize() { _sum = 0; }

void SubdomainPhotonEmission::execute() {

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

        _sum += std::accumulate(begin, end, 0.0);
      }
    }
  }
}

void SubdomainPhotonEmission::finalize() { comm().sum(_sum); }

PostprocessorValue SubdomainPhotonEmission::getValue() const { return _sum; }
