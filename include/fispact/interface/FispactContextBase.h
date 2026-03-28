#pragma once
#include <memory>
#include <unordered_map>
#include <vector>

class IFispactInputDataBase {
public:
  virtual void setFlux(const std::vector<double> &flux_energy_groups,
                       const std::vector<double> &flux) = 0;

  virtual void setFluxWallLoading(double wall_loading) = 0;

  virtual void setFluxName(std::string flux_name) = 0;

  virtual void setDensity(double density) = 0;

  virtual void setMassTotal(double total_mass) = 0;

  virtual void setMass(const std::vector<int> &atomicnumbers,
                       const std::vector<double> &percentages) = 0;

  virtual void setFuel(const std::vector<int> &zais,
                       const std::vector<double> &values) = 0;

  virtual void setSchedule(const std::vector<double> &deltatime,
                           const std::vector<double> &fluxamp) = 0;

  virtual std::pair<std::vector<double>, std::vector<double>> getSchedule() = 0;

  virtual void setAtomsThreshold(double threshold) = 0;

  virtual void setSolverTolerance(double rtol, double atol) = 0;
};

class IFispactOutputDataBase {
public:
  virtual std::vector<double> getGammaSpectrumBins(int inv_index) = 0;

  virtual std::vector<double> getGammaSpectrumBoundaries(int inv_index) = 0;
};

class IFispactUtilsBase {
public:
  virtual int GetZai(std::string nuclidename) = 0;

  virtual int GetAtomicNumberFromElementName(std::string elementname) = 0;

  virtual std::vector<double> getNeutronEnergyBounds(size_t n_groups) = 0;

  virtual std::vector<double> getPhotonEnergyBounds(size_t n_groups) = 0;

  virtual std::vector<double>
  GroupConvertByEnergy(const std::vector<double> &inbounds,
                       const std::vector<double> &invals,
                       const std::vector<double> &outbounds) = 0;

  virtual std::vector<double>
  GroupConvertByLethargy(const std::vector<double> &inbounds,
                         const std::vector<double> &invals,
                         const std::vector<double> &outbounds) = 0;
};

class FispactContextBase {
public:
  virtual void globalInitialise() = 0;
  virtual void globalFinalise() = 0;
  virtual void process() = 0;

  virtual void setNuclearData(
      std::unordered_map<std::string, std::string> nuclear_data_paths) = 0;

  virtual size_t getNuclearDataCrossSections() = 0;
  virtual IFispactInputDataBase &getInput() { return *_i_input_data; }
  virtual IFispactOutputDataBase &getOutput() { return *_i_output_data; }
  virtual IFispactUtilsBase &getUtils() { return *_i_utils; }

protected:
  std::unique_ptr<IFispactInputDataBase> _i_input_data;
  std::unique_ptr<IFispactOutputDataBase> _i_output_data;
  std::unique_ptr<IFispactUtilsBase> _i_utils;
};
