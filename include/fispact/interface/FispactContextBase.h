#pragma once
#include <memory>
#include <unordered_map>
#include <vector>

class FispactOutputNuclideDataBase {
public:
  virtual ~FispactOutputNuclideDataBase() = default;

  virtual std::string getElement() const = 0;

  virtual std::string getState() const = 0;

  virtual int getIsotope() const = 0;

  virtual int getZAI() const = 0;

  // The half life (s)
  virtual double getHalfLife() const = 0;
  // The number of atoms
  virtual double getAtoms() const = 0;
  // The grams (g)
  virtual double getGrams() const = 0;
  // The activity (Bq)
  virtual double getActivity() const = 0;
  // The alpha fraction of the activity (Bq)
  virtual double getAlphaActivity() const = 0;
  // The beta fraction of the activity (Bq)
  virtual double getBetaActivity() const = 0;
  // The gamma fraction of the activity (Bq)
  virtual double getGammaActivity() const = 0;
  // The total heat (kW)
  virtual double getTotalHeat() const = 0;
  // The alpha heat (kW)
  virtual double getAlphaHeat() const = 0;
  // The beta heat (kW)
  virtual double getBetaHeat() const = 0;
  // The gamma heat (kW)
  virtual double getGammaHeat() const = 0;
  // The dose rate (Sv/hr)
  virtual double getDoseRate() const = 0;
  // The ingestion (Sv)
  virtual double getIngestion() const = 0;
  // The inhalation (Sv)
  virtual double getInhalation() const = 0;

private:
};

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

  virtual std::vector<std::unique_ptr<FispactOutputNuclideDataBase>>
  getInventoryNuclides(int inventory_index) = 0;

  enum FispactOutputs {
    FISPACT_OUTPUT_DATA_INVENTORY_IRRAD_TIME,
    FISPACT_OUTPUT_DATA_INVENTORY_COOL_TIME,
    FISPACT_OUTPUT_DATA_INVENTORY_TOTAL_ACTIVITY,
    FISPACT_OUTPUT_DATA_INVENTORY_ALPHA_ACTIVITY,
    FISPACT_OUTPUT_DATA_INVENTORY_BETA_ACTIVITY,
    FISPACT_OUTPUT_DATA_INVENTORY_GAMMA_ACTIVITY,
    FISPACT_OUTPUT_DATA_INVENTORY_TOTAL_HEAT,
    FISPACT_OUTPUT_DATA_INVENTORY_ALPHA_HEAT,
    FISPACT_OUTPUT_DATA_INVENTORY_BETA_HEAT,
    FISPACT_OUTPUT_DATA_INVENTORY_GAMMA_HEAT,
    FISPACT_OUTPUT_DATA_INVENTORY_TOTAL_MASS,
    FISPACT_OUTPUT_DATA_INVENTORY_TOTAL_ATOMS,
    FISPACT_OUTPUT_DATA_INVENTORY_INGESTION,
    FISPACT_OUTPUT_DATA_INVENTORY_INHALATION,
    FISPACT_OUTPUT_DATA_INVENTORY_FLUX_AMP
  };

  virtual std::pair<std::vector<int>, std::vector<double>>
  getSortedInventory(int inv_index, FispactOutputs key) const = 0;
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
