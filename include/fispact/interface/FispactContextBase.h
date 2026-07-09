#pragma once
#include "FizzyEnums.h"
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

  virtual double
  getQuantity(nuclide_quantities::NuclideQuantitiesEnum quantity) const {

    switch (quantity) {
    case (nuclide_quantities::ATOMS):
      return getAtoms();
    case (nuclide_quantities::GRAMS):
      return getGrams();
    case (nuclide_quantities::ACTIVITY):
      return getActivity();
    case (nuclide_quantities::ALPHA_ACTIVITY):
      return getAlphaActivity();
    case (nuclide_quantities::BETA_ACTIVITY):
      return getBetaActivity();
    case (nuclide_quantities::GAMMA_ACTIVITY):
      return getGammaActivity();
    case (nuclide_quantities::TOTAL_HEAT):
      return getTotalHeat();
    case (nuclide_quantities::ALPHA_HEAT):
      return getAlphaHeat();
    case (nuclide_quantities::BETA_HEAT):
      return getBetaHeat();
    case (nuclide_quantities::GAMMA_HEAT):
      return getGammaHeat();
    case (nuclide_quantities::DOSE):
      return getDoseRate();
    default:
      throw std::invalid_argument(
          "Invalid quantity requested from FISPACT nuclide inventory");
    }
  }

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

  virtual std::pair<std::vector<int>, std::vector<double>>
  getSortedInventory(int inv_index,
                     inventory_outputs::InventoryOutputsEnum key) const = 0;

  virtual double
  getInventoryValue(int inv_index,
                    inventory_outputs::InventoryOutputsEnum key) = 0;

  virtual int findInventoryIndex(int inv_index, int zai) = 0;

  virtual bool findInventoryExists(int inv_index, int zai) = 0;
};

class IFispactUtilsBase {
public:
  virtual int GetZai(std::string nuclidename) = 0;

  virtual std::string getNuclideName(int zai) = 0;

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
