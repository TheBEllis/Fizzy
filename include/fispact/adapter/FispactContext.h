#pragma once
#include "FispactContextBase.h"

#include <fstream>
#include <memory>
/// Fispact includes
#include "fispactcompute.hpp"
#include "fispactelementaldata.hpp"
#include "fispactgroupconvert.hpp"
#include "fispactgroupstructures.hpp"
#include "fispactinputdata.hpp"
#include "fispactmonitor.hpp"
#include "fispactnucleardata.hpp"
#include "fispactoutputdata.hpp"
#include "fispactutil.hpp"

// IFispactInputData
class IFispactInputData : public IFispactInputDataBase {

public:
  IFispactInputData(fispact::FispactMonitor &monitor)
      : IFispactInputDataBase(), _monitor(monitor), _input(_monitor) {}

  virtual void setFlux(const std::vector<double> &flux_energy_groups,
                       const std::vector<double> &flux) {
    _input.setFlux(flux_energy_groups, flux);
  }

  virtual void setFluxWallLoading(double wall_loading) {
    _input.setFluxWallLoading(wall_loading);
  }

  virtual void setFluxName(std::string flux_name) {
    _input.setFluxName(flux_name);
  }

  virtual void setDensity(double density) { _input.setDensity(density); }

  virtual void setMassTotal(double total_mass) {
    _input.setMassTotal(total_mass);
  }

  virtual void setMass(const std::vector<int> &atomicnumbers,
                       const std::vector<double> &percentages) {
    _input.setMass(atomicnumbers, percentages);
  }

  virtual void setFuel(const std::vector<int> &zais,
                       const std::vector<double> &values) {
    _input.setFuel(zais, values);
  }

  virtual void setSchedule(const std::vector<double> &deltatime,
                           const std::vector<double> &fluxamp) {
    _input.setSchedule(deltatime, fluxamp);
  }

  virtual std::pair<std::vector<double>, std::vector<double>> getSchedule() {
    return _input.getSchedule();
  }

  virtual void setAtomsThreshold(double threshold) {
    _input.setAtomsThreshold(threshold);
  }

  fispact::InputData &getInput() { return _input; };

private:
  fispact::FispactMonitor &_monitor;
  fispact::InputData _input;
};
//~IFispactInputData

// IFispactOutputData
class IFispactOutputData : public IFispactOutputDataBase {

public:
  IFispactOutputData(fispact::FispactMonitor &monitor)
      : IFispactOutputDataBase(), _monitor(monitor), _output(_monitor) {}

  virtual std::vector<double> getGammaSpectrumBins(int inv_index) {
    return _output.getGammaSpectrumBins(inv_index);
  };

  virtual std::vector<double> getGammaSpectrumBoundaries(int inv_index) {
    return _output.getGammaSpectrumBoundaries(inv_index);
  };

  fispact::OutputData &getOutput() { return _output; };

private:
  fispact::FispactMonitor &_monitor;
  fispact::OutputData _output;
};
//~IFispactOutputData

class IFispactUtils : public IFispactUtilsBase {

public:
  IFispactUtils(fispact::FispactMonitor &monitor)
      : IFispactUtilsBase(), _monitor(monitor) {

    _neutron_energy_groups_map[709] = fispact::groups::G709();
    _neutron_energy_groups_map[1102] = fispact::groups::G1102();

    _photon_energy_groups_map[22] = {0.0,   1.0e4, 1.0e5, 2.0e5, 4.0e5, 1.0e6,
                                     1.5e6, 2.0e6, 2.5e6, 3.0e6, 3.5e6, 4.0e6,
                                     4.5e6, 5.0e6, 5.5e6, 6.0e6, 6.5e6, 7.0e6,
                                     7.5e6, 8.0e6, 1.0e7, 1.2e7, 1.4e7};
    _photon_energy_groups_map[24] = {
        1.000e-11, 1.000e+4, 2.000e+4, 5.000e+4, 1.000e+5, 2.000e+5, 3.000e+5,
        4.000e+5,  6.000e+5, 8.000e+5, 1.000e+6, 1.220e+6, 1.440e+6, 1.660e+6,
        2.000e+6,  2.500e+6, 3.000e+6, 4.000e+6, 5.000e+6, 6.500e+6, 8.000e+6,
        1.000e+7,  1.200e+7, 1.400e+7, 2.000e+7};
  }

  virtual int GetZai(std::string nuclidename) {
    return fispact::util::GetZai(_monitor, nuclidename);
  }

  virtual int GetAtomicNumberFromElementName(std::string elementname) {
    return fispact::util::GetAtomicNumberFromElementName(_monitor, elementname);
  }

  virtual std::vector<double> getNeutronEnergyBounds(size_t n_bins) {
    return _neutron_energy_groups_map[n_bins];
  }

  virtual std::vector<double> getPhotonEnergyBounds(size_t n_bins) {
    return _photon_energy_groups_map[n_bins];
  }

  virtual std::vector<double>
  GroupConvertByEnergy(const std::vector<double> &inbounds,
                       const std::vector<double> &invals,
                       const std::vector<double> &outbounds) {
    return fispact::groupconvert::GroupConvertByEnergy(_monitor, inbounds,
                                                       invals, outbounds);
  }

  virtual std::vector<double>
  GroupConvertByLethargy(const std::vector<double> &inbounds,
                         const std::vector<double> &invals,
                         const std::vector<double> &outbounds) {
    return fispact::groupconvert::GroupConvertByLethargy(_monitor, inbounds,
                                                         invals, outbounds);
  }

private:
  fispact::FispactMonitor &_monitor;
  std::unordered_map<size_t, std::vector<double>> _neutron_energy_groups_map;
  std::unordered_map<size_t, std::vector<double>> _photon_energy_groups_map;
};

class FispactContext : public FispactContextBase {
public:
  FispactContext() : FispactContextBase(), _monitor("log_name"), _nd(_monitor) {

    _i_input_data = std::make_unique<IFispactInputData>(_monitor);
    _i_output_data = std::make_unique<IFispactOutputData>(_monitor);
    _i_utils = std::make_unique<IFispactUtils>(_monitor);
  }

  virtual void globalInitialise() { fispact::GlobalInitialise(_monitor); }

  virtual void globalFinalise() { fispact::GlobalFinalise(_monitor); }

  virtual void process() {
    auto *real_input = dynamic_cast<IFispactInputData *>(_i_input_data.get());
    auto *real_output =
        dynamic_cast<IFispactOutputData *>(_i_output_data.get());
    fispact::InputData &input = real_input->getInput();
    fispact::OutputData &output = real_output->getOutput();
    fispact::Process(input, _nd, output, _monitor, process_callback);
  }

  virtual size_t getNuclearDataCrossSections() {
    return _nd.getReactionXS(0, 0).size();
  }
  virtual void setNuclearData(
      std::unordered_map<std::string, std::string> nuclear_data_paths) {

    fispact::io::NuclearDataReader nd_reader(_monitor);

    for (auto &[data, path] : nuclear_data_paths) {
      if (data == "ND_IND_NUC_KEY") {
        nd_reader.setPath(FISPACT_ND_IND_NUC_KEY, path);
        continue;
      }

      if (data == "ND_HAZARDS_KEY") {
        nd_reader.setPath(FISPACT_ND_HAZARDS_KEY, path);
        continue;
      }

      if (data == "ND_ABSORP_KEY") {

        nd_reader.setPath(FISPACT_ND_ABSORP_KEY, path);
        continue;
      }

      if (data == "ND_CLEAR_KEY") {

        nd_reader.setPath(FISPACT_ND_CLEAR_KEY, path);
        continue;
      }

      if (data == "ND_A2DATA_KEY") {
        nd_reader.setPath(FISPACT_ND_A2DATA_KEY, path);
        continue;
      }

      if (data == "ND_ENBINS_KEY") {
        nd_reader.setPath(FISPACT_ND_ENBINS_KEY, path);
        continue;
      }

      if (data == "ND_DECAY_KEY") {
        nd_reader.setPath(FISPACT_ND_DECAY_KEY, path);
        continue;
      }

      if (data == "ND_DK_ENDF_KEY") {
        nd_reader.setPath(FISPACT_ND_DK_ENDF_KEY, path);
        continue;
      }

      if (data == "ND_PROB_TAB_KEY") {
        nd_reader.setPath(FISPACT_ND_PROB_TAB_KEY, path);
        continue;
      }

      if (data == "ND_ASSCFY_KEY") {
        nd_reader.setPath(FISPACT_ND_ASSCFY_KEY, path);
        continue;
      }

      if (data == "ND_FISSYLD_KEY") {
        nd_reader.setPath(FISPACT_ND_FISSYLD_KEY, path);
        continue;
      }

      if (data == "ND_FY_ENDF_KEY") {
        nd_reader.setPath(FISPACT_ND_FY_ENDF_KEY, path);
        continue;
      }

      if (data == "ND_SF_ENDF_KEY") {
        nd_reader.setPath(FISPACT_ND_SF_ENDF_KEY, path);
        continue;
      }

      if (data == "ND_SP_ENDF_KEY") {
        nd_reader.setPath(FISPACT_ND_SP_ENDF_KEY, path);
        continue;
      }

      if (data == "ND_XS_EXTRA_KEY") {
        nd_reader.setPath(FISPACT_ND_XS_EXTRA_KEY, path);
        continue;
      }

      if (data == "ND_CROSSEC_KEY") {
        nd_reader.setPath(FISPACT_ND_CROSSEC_KEY, path);
        continue;
      }

      if (data == "ND_CROSSUNC_KEY") {
        nd_reader.setPath(FISPACT_ND_CROSSUNC_KEY, path);
        continue;
      }

      if (data == "ND_XS_ENDF_KEY") {
        std::cout << "HERE" << std::endl;
        nd_reader.setPath(FISPACT_ND_XS_ENDF_KEY, path);
        continue;
      }

      if (data == "ND_XS_ENDFB_KEY") {
        nd_reader.setUseXSBinary(true);
        nd_reader.setPath(FISPACT_ND_XS_ENDFB_KEY, path);
        continue;
      }
    }
    // Load the nuclear data
    nd_reader.load(_nd, load_callback);
  }

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

private:
  fispact::FispactMonitor _monitor;
  fispact::NuclearData _nd;
};
