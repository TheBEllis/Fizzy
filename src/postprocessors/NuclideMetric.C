#include "NuclideMetric.h"
#include <numeric>

registerMooseObject("FizzyApp", NuclideMetric);

InputParameters NuclideMetric::validParams() {

  InputParameters params = GeneralPostprocessor::validParams();

  params.addRequiredParam<std::vector<dof_id_type>>(
      "element_ids",
      "Global ID's off all the elements to sum photon emission over.");
  return params;
}

NuclideMetric::NuclideMetric(const InputParameters &params)
    : FispactPostprocessor(params),
      _element_ids(getParam<std::vector<dof_id_type>>("element_ids")) {}

void NuclideMetric::initialize() { _sum = 0; }

void NuclideMetric::execute() {

  PhotonSpectra *spectra = getFispactProblem().getPhotonSpectra();

  std::unordered_map<uint64_t, uint64_t> &local_element_index =
      getFispactProblem().getLocalElemIndexMap();

  size_t inventory_index = getFispactInventoryIdx();

  for (auto &elem_id : _element_ids) {

    if (getFispactProblem().mesh().elemPtr(elem_id)->processor_id() ==
        processor_id()) {

      std::vector<double>::iterator begin = spectra->spectrum_begin(
          inventory_index, local_element_index[elem_id]);
      std::vector<double>::iterator end =
          spectra->spectrum_end(inventory_index, local_element_index[elem_id]);

      _sum += std::accumulate(begin, end, 0.0);
    }
  }
}

void NuclideMetric::finalize() { comm().sum(_sum); }

PostprocessorValue NuclideMetric::getValue() const { return _sum; }
