#pragma once

#include "ExternalProblem.h"

#include "FISPACTMaterial.h"

#include "HDF5Utils.h"

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
#include <string>
#include <unordered_map>

namespace fp = fispact;

#ifdef LIBMESH_HAVE_BOOST
namespace bi = boost::interprocess;
#endif

class FispactProblem : public ExternalProblem {

  /**
   * Struct to store FISPACT Material definitions
   *
   */
  struct MaterialDefinition {
    std::string _mat_name;
    std::vector<std::pair<std::string, double>> _mat_atomic_composition;
    double _mat_density;
  };

public:
  FispactProblem(const InputParameters &params);

  static InputParameters validParams();

  // virtual void initialSetup() override;
  virtual void externalSolve() override;
  virtual void syncSolutions(ExternalProblem::Direction direction) override;
  virtual bool converged(unsigned int) override { return true; }

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
  void setNuclearData(const std::string &fp_nuclear_data_uo);

  /**
   *
   * Read neutron flux spectra from an OpenMC statepoint file
   * or a similarly organised HDF5 file.
   * @param[in] filename the name of the statepoint/HDF5 file to read neutron
   * flux from
   * @param[in] elem_id the mesh element ID whose neutron flux spectra we wish
   * to read
   * @param[in] tally_id the OpenMC tally ID containing the neutron flux spectra
   * @param[num_neutron_bins] the number of energy bins the neutron flux is
   * binned over
   * @return Vector of length num_neutron_bins containing the neutron flux
   * spectra from elemenet elem_id
   *
   */
  std::vector<double> readElementNeutronFlux(const std::string &filename,
                                             const dof_id_type &elem_id,
                                             const int &tally_id) const;

  /**
   * Method to write calculated photon flux to a HDF5 file, primarily for
   * debugging purposes.
   * @param[in] filename the name of the resulting HDF5 file containing the
   * photon flux
   *
   */
  void writePhotonFlux(const std::string &filename);

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
  void setFispactSchedule(fp::InputData &input) const;

  /**
   * Converts FISPACT gamma spectra outputs from MeV s^-1 to cm^-3 s^-1
   * @param[in] input FISPACT input
   * @param[in] photon_spectra
   * @param[in] photon_flux_bins
   * @param[out] photons_per_cc_per_s
   */
  void convertGammaEvToCount(const fp::InputData &input,
                             const std::vector<double> &photon_spectra,
                             std::vector<double> &photons_per_cc_per_s);
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
  void insertElementStrength(const libMesh::Elem *element,
                             const std::vector<double> element_flux);

  /**
   *
   */
  void updateLocalDomainStrength(const libMesh::Elem *element);

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

  /**
   *
   */
  void setNeutronBins();

  /**
   *
   */
  void setPhotonBins(const fp::OutputData &fispact_output);

  /**
   *
   */
  std::vector<double> &getNeutronBins() {
    if (_neutron_bins.empty()) {
      mooseError("_neutron_bins not set, so could not be gotten");
    }
    return _neutron_bins;
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
  /**
   *
   */
  int calculateMemorySize();

  /**
   * Method to generate name of the interprocess memory segment generated by
   * this MPI rank
   * @param[out] ipc_name Reference to the string we want to populate
   */
  const std::string generateInterprocessName();

  /// -- Interprocess bits --
#ifdef LIBMESH_HAVE_BOOST
  bi::managed_shared_memory _segment;
#endif

  /// FISPACT monitor
  fp::FispactMonitor _fp_monitor;

  /// FISPACT nuclear data
  fp::NuclearData _fp_nuclear_data;

  /// FISPACT neutron flux
  std::unordered_map<int, std::vector<double>> _neutron_fluxes;

  /// FISPACT Photon fluxes
  std::unordered_map<int, std::vector<double>> _photon_fluxes;

  ///
  std::unordered_map<int, double> _element_strengths;

  /// Fispact nuclear data path
  std::string _fp_nuclear_data_uo;

  /// hdf5 filename for neutron flux
  std::string _neutron_flux_filename;

  /// path to neutron flux array in hdf5 file
  int _neutron_flux_tally_id;

  /// hdf5 filename for photon flux
  std::string _photon_flux_filename;

  bool _materials_from_xml;

  /// Filename of xml file to read materials from
  std::string _materials_xml_file;

  /// Mappings from material name to material definitions
  std::unordered_map<std::string, MaterialDefinition> _mat_definitions;

  ///
  const FISPACTMaterial &getElementMaterial(dof_id_type &elem_id);

  ///
  std::string _neutron_bin_type;

  ///
  std::vector<double> _neutron_bins;

  ///
  std::vector<double> _photon_bins;

  ///
  int _num_neutron_bins;

  ///
  int _num_photon_bins;

  ///
  std::string _schedule_uo_name;

  ///
  double _local_domain_strength;

  ///
  double _total_domain_strength;

  ///
  bool _write_photon_flux;

  ///
  bool _comm_photon_flux;

  /// name given to the interprocess memory segment generated by this rank
  const std::string _interprocess_segment_name;

  const std::string &_molar_mass_data_filename;

  /// Map from element ZAI to molar mass (g/mol)
  std::unordered_map<int, double> _molar_mass_map;
};
