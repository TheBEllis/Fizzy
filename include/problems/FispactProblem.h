#pragma once
#include "ExternalProblem.h"

#include "fispactcompute.hpp"
#include "fispactelementaldata.hpp"
#include "fispactgroupconvert.hpp"
#include "fispactgroupstructures.hpp"
#include "fispactinputdata.hpp"
#include "fispactnucleardata.hpp"
#include "fispactoutputdata.hpp"
#include "fispactutil.hpp"
#include "libmesh/elem.h"
#include <hdf5/openmpi/H5Ipublic.h>
#include <string>
#include <unordered_map>

// Include for interprocess communication data structure
#include "PhotonSharingData.h"

// Use fp as short for fispact
namespace fp = fispact;

#ifdef LIBMESH_HAVE_BOOST
namespace bi = boost::interprocess;
#endif

class FispactProblem : public ExternalProblem {
  /// Struct to store Material definitions
  struct MaterialDefinition {
    std::string _mat_name;
    std::vector<std::pair<std::string, double>> _mat_atomic_composition;
    double _mat_density;
  };

public:
  // Constructor for FispactProblem
  FispactProblem(const InputParameters &params);

  static InputParameters validParams();

  // virtual void initialSetup() override;
  virtual void externalSolve() override;
  virtual void syncSolutions(ExternalProblem::Direction direction) override;
  virtual bool converged(unsigned int) override { return true; }

private:
  static void load_callback(std::string key, std::string path, int i, int t) {
    std::cout << "\33[2K\r" << key << ": " << path << " [" << i << "/" << t
              << "]" << std::flush;
  }

  static void process_callback(std::string process_name, int i, int t) {
    std::cout << "\33[2K\r [" << i << "/" << t << "] " << process_name
              << std::flush;
  }

  void setNuclearData(std::string nd_base_path);

  /// Read a neutron flux spectra from a hdf5 file
  std::vector<double> readNeutronFluxFromHDF5(const std::string &filename);

  /// Write output photon flux to HDF5
  void writePhotonFluxToHDF5(const std::string &filename,
                             fp::OutputData &fispact_output);

  void writePhotonFluxBins(hid_t file_id, fp::OutputData &fispact_output);

  void setFispactInputData(fp::FispactMonitor &monitor, fp::InputData &input,
                           MaterialDefinition &material,
                           std::vector<double> &neutron_flux,
                           const std::vector<double> &bins, double volume);

  /// Generate a log file name for the fispact logs
  std::string fispactLogName();

  // Initialise FISPACT monitor object
  void initFispactMonitor(std::string);

  void readNeutronFluxFromHDF5(std::string filename, std::string tally_dir);

  void read_material_xml_data();

  /**
   * Method to check that for all mesh subdomains, a corresponding
   * FispactMaterial exists. This allows for the program to error out at the
   * start as opposed to running into a subdomain lacking a material half way
   * through the solve.
   */
  void checkMaterialsExist();

  /**
   * Retrieve Fispact radiation schedule from FispactSchedule UserObject
   */
  void setFispactSchedule(fp::InputData &input);

  void convertGammaEvToCount(fp::InputData &input,
                             const std::vector<double> &photon_spectra,
                             const std::vector<double> &photon_flux_bins,
                             std::vector<double> &photons_per_cc_per_s);

  double calculateElementStrength(const libMesh::Elem *element,
                                  const std::vector<double> element_flux);

  void insertElementStrength(const libMesh::Elem *element,
                             const std::vector<double> element_flux);

  void updateLocalDomainStrength(const libMesh::Elem *element);

  int calculateMemorySize();

  void printInventoryByHeat(fp::OutputData &output, int timestep_index,
                            std::ostream &stream = std::cout);

  void printInventoryByMass(fp::OutputData &output, int timestep_index,
                            std::ostream &stream);

  void printInvData(fp::OutputData &output, std::ostream &stream);

  void getTotalDomainStrength();

  double extractHalflifeFromNuc(fp::NuclearData &nuclear_data, int zai);

  bool isFlux(int elem_id);

  // -- Interprocess bits --

#ifdef LIBMESH_HAVE_BOOST
  bi::managed_shared_memory _segment;
#endif

  /// FISPACT monitor
  fp::FispactMonitor _fp_monitor;

  /// FISPACT nuclear data
  fp::NuclearData _fp_nuclear_data;

  /// FISPACT neutron flux
  std::unordered_map<int, std::vector<double>> _neutron_fluxes;

  // FISPACT Photon fluxes
  std::unordered_map<int, std::vector<double>> _photon_fluxes;

  std::unordered_map<int, double> _element_strengths;

  /// hdf5 filename for neutron flux
  std::string _neutron_flux_filename;

  /// path to neutron flux array in hdf5 file
  std::string _neutron_flux_hdf5_path;

  /// hdf5 filename for photon flux
  std::string _photon_flux_filename;

  void setNeutronBins();
  ///
  bool _materials_from_xml;

  /// Boolean to determine whether a rank is still running FISPACT calculations
  bool _calculating = true;

  /// Filename of xml file to read materials from
  std::string _materials_xml_file;

  /// Mappings from material name to material definitions
  std::unordered_map<std::string, MaterialDefinition> _mat_definitions;

  MaterialDefinition &getElementMaterial(int &elem_id);

  std::string _neutron_bin_type;

  std::vector<double> _neutron_bins;

  std::vector<double> _photon_bins;

  int _num_neutron_bins;

  std::string _schedule_uo_name;

  double _local_domain_strength;

  double _total_domain_strength;

  bool _write_photon_flux;

  bool _comm_photon_flux;
};
