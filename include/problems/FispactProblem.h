#pragma once

#include "ExternalProblem.h"

#include "../../include/userobjects/FispactNuclearDataPaths.C"
#include "FISPACTMaterial.h"
#include "FispactFluxInput.h"
#include "FispactSchedule.h"

#include "HDF5Utils.h"

#include "MooseTypes.h"
#include "PhotonSpectra.h"

// Fispact includes
#include "fispactcompute.hpp"
#include "fispactmonitor.hpp"
#include "fispactutil.hpp"

#include "fispactelementaldata.hpp"
#include "fispactgroupconvert.hpp"
#include "fispactgroupstructures.hpp"
#include "fispactinputdata.hpp"
#include "fispactnucleardata.hpp"
#include "fispactoutputdata.hpp"
#include "fispactutil.hpp"

#include "fispactinputdata.hpp"
#include "fispactmonitor.hpp"
#include "fispactnucleardata.hpp"
#include "fispactoutputdata.hpp"
// Include for interprocess communication data structure
#include "PhotonSharingData.h"
#include "libmesh/id_types.h"
#include <string>
#include <unordered_map>

namespace fp = fispact;

#ifdef LIBMESH_HAVE_BOOST
namespace bi = boost::interprocess;
#endif

class FispactProblem : public ExternalProblem {

public:
  FispactProblem(const InputParameters &params);

  ~FispactProblem();

  static InputParameters validParams();

  virtual void initialSetup() override;

  virtual void externalSolve() override;
  virtual void syncSolutions(ExternalProblem::Direction direction) override;
  virtual bool converged(unsigned int) override { return true; }

  void timestepSetup() override;

private:
  /**
   *
   *
   *
   */
  static void load_callback(std::string key, std::string path, int i, int t) {
    std::cout << "\33[2K\r" << key << ": " << path << " [" << i << "/" << t
              << "]" << std::flush;
  }

  /**
   *
   *
   *
   */
  static void process_callback(std::string process_name, int i, int t) {
    std::cout << "\33[2K\r [" << i << "/" << t << "] " << process_name
              << std::flush;
  }

  /**
   * Set the nuclear data for FISPACT
   * @param[in] nd_base_path The directory containing the various nuclear data
   * resources required for fispact
   *
   */
  void setNuclearData();

  /**
   * Method to write calculated photon flux to a HDF5 file, primarily for
   * debugging purposes.
   * @param[in] filename the name of the resulting HDF5 file containing the
   * photon flux
   *
   */
  void writePhotonFlux(const std::string &filename,
                       const std::vector<double> &inventory_times);

  /**
   * Method used within @ref writePhotonFlux to write the photon energy bins
   * used by FISPACT to the same HDF5 file
   * @param[in] file_id the HDF5 file identifier for the file we are writing to
   * @param[in] photon_bins a vector of doubles containing the photon energy
   * bins used
   * @param[in] parallel indicates whether or not to use the parallel HDF5
   * driver
   */
  void writePhotonFluxBins(const hid_t &file_id, const bool parallel);

  /**
   *
   * Method used to set input parameters for a FISPACT activation calculation
   * @param[in] monitor the FISPACT monitor object required by all FISPACT
   * methods
   * @param[in] material the material definition for the FISPACT calculation
   * @param[in] neutron_flux the energy binned neutron flux
   * @param[in] bins the energy bin boundaries for neutron flux
   * @param[in] volume the volume in m^3 for the mesh element this FISPACT
   * calculation corresponds to
   * @param[out] input the now correctly setup FISPACT input object
   */
  void setFispactInputData(const fp::FispactMonitor &monitor,
                           const FISPACTMaterial &material,
                           const std::vector<double> &neutron_flux,
                           const double &volume, fp::InputData &input) const;

  /// Generate a log file name for the fispact logs
  std::string fispactLogName();

  /**
   *
   *
   */
  void read_material_xml_data();

  /**
   * Method to check that for all mesh subdomains, a corresponding
   * FispactMaterial exists. This allows for the program to error out at the
   * start as opposed to running into a subdomain lacking a material half way
   * through the solve.
   */
  void checkMaterialsExist() const;

  /**
   * Retrieve Fispact radiation schedule from FispactSchedule UserObject
   */
  void setFispactSchedule(fp::InputData &input, const double &volume,
                          const double &neutron_flux_sum) const;

  /**
   * Converts FISPACT gamma spectra outputs from MeV s^-1 to cm^-3 s^-1
   * @param[in] input FISPACT input
   * @param[in] photon_energy_spectra_ev Photon energy spectra as output by
   * FISPACT
   * @param[out] photon_energy_spectra_per_cc_s Photon energy spectra in
   * photons/cc-s
   */
  void
  convertGammaEvToCount(const fp::InputData &input,
                        const std::vector<double> &photon_energy_spectra_ev,
                        std::vector<double> &photon_energy_spectra_per_cc_s);
  /**
   * Calculates the total strength of one elements photon source term
   * @param[in] element a ptr to the libmesh element whose strength we are
   * calculating
   * @param[in] element_flux the vector representing the energy binned photon
   * flux for the relevent mesh element
   * @return the total strength of the photon source term for this element
   *
   */
  double calculateElementStrength(const libMesh::Elem *element,
                                  const std::vector<double> element_flux);
  /**
   *
   */
  void insertElementStrength(const int inv_index, const libMesh::Elem *element,
                             const std::vector<double> element_flux);

  /**
   *
   */
  void calculateLocalDomainStrength();

  /**
   *
   */
  void printInventoryByHeat(fp::OutputData &output, int timestep_index,
                            std::ostream &stream = std::cout);

  /**
   *
   */
  void printInventoryByMass(fp::OutputData &output, int timestep_index,
                            std::ostream &stream);

  /**
   *
   */
  void printInvData(fp::OutputData &output, std::ostream &stream);

  /**
   *
   */
  void getTotalDomainStrength();

  /**
   *
   */
  double extractHalflifeFromNuc(fp::NuclearData &nuclear_data, int zai);

  /**
   *
   */
  bool isFlux(const std::vector<double> &flux) const;

  void resolveFispactUserObjects();

  /**
   * Check if neutron group structure is consistent between the passed in flux
   * data and the read in nuclear data
   */
  void checkForEnergyGroupConsistency();

  /**
   *
   */
  std::vector<double> &getNeutronBins() {
    if (_flux_energy_groups.empty()) {
      mooseError("_flux_energy_groups not set, so could not be gotten");
    }
    return _flux_energy_groups;
  }

  /**
   * Return a vector populated with the number of atoms of each isotope in a
   * given material
   */
  std::vector<double> getMaterialAtoms(const FISPACTMaterial &material);

  /**
   * Load in molar mass data from HDF5 file set using input params
   */
  void loadMolarMasses();

  /**
   * Method to calculate the number of atoms in a given mass of substance
   */
  double getNumAtoms(const double &mass, const double &molar_mass,
                     const double &avogadro) const;

  /**
   *
   */
  std::vector<double> &getPhotonBins() {
    if (_photon_bins.empty()) {
      mooseError("_photon_bins not set, so could not be gotten");
    }
    return _photon_bins;
  }

  size_t photonEnergySpectraIdx(size_t inv_idx, dof_id_type elem_id) {
    return ((inv_idx * (_mesh.nActiveLocalElem() * _n_photon_bins)) +
            (_local_elem_index[elem_id] * _n_photon_bins));
  }

  /**
   *
   */
  int calculateMemorySize();

  void convertFluxEnergyGroups(std::vector<double> &input_flux);
  /**
   * Method to generate name of the interprocess memory segment generated by
   * this MPI rank
   * @param[out] ipc_name Reference to the string we want to populate
   */
  const std::string generateInterprocessName();

  /// -- Interprocess bits --
#ifdef LIBMESH_HAVE_BOOST
  bi::managed_shared_memory _segment;

  PhotonSharingData *_photon_sharing_instance;
#endif

  /// FISPACT monitor
  fp::FispactMonitor _fp_monitor;

  /// FISPACT nuclear data
  fp::NuclearData _fp_nuclear_data;

  /// FISPACT neutron flux
  std::unordered_map<int, std::vector<double>> _neutron_fluxes;

  /// FISPACT Photon fluxes
  /// Indexed by time major, element id minor
  // std::vector<double> _photon_energy_spectra;

  std::unique_ptr<PhotonSpectra> _photon_energy_spectra;

  /// Vector to store all local element strengths over all inventory times
  /// Indexed by time major, element id minor
  std::vector<double> _element_strengths;

  /// hdf5 filename for neutron flux
  std::string _neutron_flux_filename;

  /// path to neutron flux array in hdf5 file
  size_t _neutron_flux_tally_id;

  /// hdf5 filename for photon flux
  std::string _photon_flux_filename;

  bool _materials_from_xml;

  /// Filename of xml file to read materials from
  std::string _materials_xml_file;

  ///
  const FISPACTMaterial &getElementMaterial(dof_id_type &elem_id);

  ///
  // const size_t _input_neutron_bin_structure;

  ///
  std::vector<double> _flux_energy_groups;

  ///
  // std::vector<double> _input_neutron_bins;

  ///
  std::vector<double> _photon_bins;

  ///
  // uint64_t _n_input_neutron_bins;

  ///
  uint64_t _n_photon_bins;

  /// Number of FISPACT inventories
  uint64_t _n_inventories;

  ///
  UserObjectName _fp_schedule_uo_name;

  ///
  UserObjectName _fp_flux_uo_name;

  /// Fispact nuclear data path
  UserObjectName _fp_nuclear_data_uo_name;

  FispactSchedule *_fp_schedule_uo;

  FispactNuclearDataPaths *_fp_nuclear_data_uo;

  FispactFluxInput *_fp_flux_input_uo;

  std::vector<FISPACTMaterial *> _fp_fispact_materials;

  /// Vector to store local domain strength over all inventory times
  std::vector<double> _local_domain_strength;

  /// Vector to store local domain strength over all inventory times
  std::vector<double> _total_domain_strength;

  ///
  bool _write_photon_flux;

  ///
  bool _comm_photon_flux;

  ///
  bool _solved;

  ///
  bool _convert_energy_groups = false;

  /// name given to the interprocess memory segment generated by this rank
  const std::string _interprocess_segment_name;

  const std::string &_molar_mass_data_filename;

  /// Map from element ZAI to molar mass (g/mol)
  std::unordered_map<int, double> _molar_mass_map;

  /// Map from global element id to "local element id"
  std::unordered_map<uint64_t, uint64_t> _local_elem_index;

  std::unordered_map<size_t, std::vector<double>> _neutron_group_map = {
      {100, fp::groups::G100()},
      {709, fp::groups::G709()},
      {1102, fp::groups::G1102()}};
};
