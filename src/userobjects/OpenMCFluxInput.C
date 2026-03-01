#include "HDF5Utils.h"
#include "InputParameters.h"
#include "OpenMCFluxInput.h"
#include "Problem.h"

registerMooseObject("FizzyApp", OpenMCFluxInput);

InputParameters OpenMCFluxInput::validParams() {
  InputParameters params = FispactFluxInput::validParams();

  params.addRequiredParam<size_t>(
      "flux_tally_id", "OpenMC tally ID of the tallies neutron flux");

  params.addRequiredParam<size_t>(
      "energy_filter_id",
      "OpenMC filter ID of the filter specifying the flux energy bins");

  params.addRequiredParam<FileName>(
      "statepoint_filename",
      "Filename of the statepoint file containing flux data");

  return params;
}

OpenMCFluxInput::OpenMCFluxInput(const InputParameters &parameters)
    : FispactFluxInput(parameters),
      _flux_tally_id(getParam<size_t>("flux_tally_id")),
      _energy_filter_id(getParam<size_t>("energy_filter_id")),
      _statepoint_filename(getParam<FileName>("statepoint_filename")) {

  readEnergyFilterBins();
  readNumEnergyBins();
}

std::vector<double> OpenMCFluxInput::getElemFlux(dof_id_type elem_id) {
  /**
   * Don't do this with HDF5 MPI driver, as we are going to call H5Dopen a
   * different number of times on each rank!
   */
  bool parallel = false;

  hid_t hdf5_file = hdf5_utils::file_open(_statepoint_filename.c_str(), 'r',
                                          parallel, comm().get());

  hid_t tallies_group = hdf5_utils::open_group(hdf5_file, "tallies");

  /// Open up the neutron flux tally group within the tallies group
  std::string tally_group_name = "tally " + std::to_string(_flux_tally_id);
  hid_t neutron_flux_tally_group =
      hdf5_utils::open_group(tallies_group, tally_group_name);

  /// Define size of selection of flux
  hsize_t dataset_dims[3];
  hsize_t hyperslab_dims[3]{_n_bins, 1, 1};
  hsize_t start[]{elem_id * _n_bins, 0, 0};
  hsize_t count[]{_n_bins, 1, 1};

  /// Number of dimensions of hyperslab to read
  hsize_t rank = 3;

  // Check this statepoint file is of
  hid_t flux_dataset =
      hdf5_utils::open_dataset(neutron_flux_tally_group, "results");
  hdf5_utils::get_shape(flux_dataset, dataset_dims);

  if (dataset_dims[0] != _n_bins * getSubProblem().mesh().getMesh().n_elem()) {
    mooseError("Provided tally dimensions do not conform to this mesh and "
               "neutron bin structure. Tally dimensions are (" +
               std::to_string(dataset_dims[0]) + "," +
               std::to_string(dataset_dims[1]) + "," +
               std::to_string(dataset_dims[2]) + ")");
  }

  std::vector<double> neutron_flux_results(_n_bins, 0);

  hdf5_utils::read_double_hyperslab(neutron_flux_tally_group, "results", rank,
                                    hyperslab_dims, start, count,
                                    neutron_flux_results.data(), parallel);

  /// Get number of tally realizations
  int n_realizations;
  hdf5_utils::read_int(neutron_flux_tally_group, "n_realizations",
                       &n_realizations, parallel, true);

  hdf5_utils::close_dataset(flux_dataset);
  hdf5_utils::close_group(tallies_group);
  hdf5_utils::close_group(neutron_flux_tally_group);

  hdf5_utils::file_close(hdf5_file);

  /**
   * Divide every value in neutron flux results vector by n_realizations to
   * get the mean value
   */
  for (auto &bin : neutron_flux_results) {
    bin /= n_realizations;
  }

  return neutron_flux_results;
}

void OpenMCFluxInput::readEnergyFilterBins() {
  /**
   * Don't do this with HDF5 MPI driver, as we are going to call H5Dopen a
   * different number of times on each rank!
   */
  bool parallel = false;

  hid_t hdf5_file = hdf5_utils::file_open(_statepoint_filename.c_str(), 'r',
                                          parallel, comm().get());

  hid_t filter_group = hdf5_utils::open_group(hdf5_file, "tallies/filters");

  /// Open up the neutron flux tally group within the tallies group
  std::string filter_name = "filter " + std::to_string(_energy_filter_id);
  hid_t energy_filter_group = hdf5_utils::open_group(filter_group, filter_name);

  // Check this statepoint file is of
  hid_t energy_bins_dataset =
      hdf5_utils::open_dataset(energy_filter_group, "bins");

  hsize_t energy_bins_dataset_dims[3];
  hdf5_utils::get_shape(energy_bins_dataset, energy_bins_dataset_dims);

  _flux_energy_groups.resize(energy_bins_dataset_dims[0]);

  hdf5_utils::read_double(energy_bins_dataset, nullptr,
                          _flux_energy_groups.data(), parallel);

  hdf5_utils::close_group(energy_filter_group);
  hdf5_utils::close_group(filter_group);
  hdf5_utils::close_dataset(energy_bins_dataset);

  hdf5_utils::file_close(hdf5_file);
}

void OpenMCFluxInput::readNumEnergyBins() {
  /**
   * Don't do this with HDF5 MPI driver, as we are going to call H5Dopen a
   * different number of times on each rank!
   */
  bool parallel = false;

  hid_t hdf5_file = hdf5_utils::file_open(_statepoint_filename.c_str(), 'r',
                                          parallel, comm().get());

  hid_t filter_group = hdf5_utils::open_group(hdf5_file, "tallies/filters");

  /// Open up the neutron flux tally group within the tallies group
  std::string filter_name = "filter " + std::to_string(_energy_filter_id);
  hid_t energy_filter_group = hdf5_utils::open_group(filter_group, filter_name);

  // Check this statepoint file is of
  int n_bins;
  hdf5_utils::read_int(energy_filter_group, "n_bins", &n_bins, parallel);

  hdf5_utils::close_group(energy_filter_group);
  hdf5_utils::close_group(filter_group);

  hdf5_utils::file_close(hdf5_file);

  _n_bins = n_bins;
}
