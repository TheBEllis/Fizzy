
#pragma once
#include "FispactContextBase.h"

#include <fstream>
#include <memory>

class IMockFispactOutputNuclideData : public FispactOutputNuclideDataBase {

public:
  IMockFispactOutputNuclideData(
      const std::string &name, const std::string &state, const int &isotope,
      const int &zai, const double &hl, const double &atoms,
      const double &grams, const double &activity, const double &a_act,
      const double &b_act, const double &g_act, const double &total_heat,
      const double &a_heat, const double &b_heat, const double &c_heat,
      const double &g_heat, const double &ingestion, const double &inhalation)
      : FispactOutputNuclideDataBase(), _name(name), _state(state), _zai(zai),
        _isotope(isotope), _hl(hl), _atoms(atoms), _grams(grams),
        _activity(activity), _a_act(a_act), _b_act(b_act), _g_act(g_act),
        _total_heat(total_heat), _a_heat(a_heat), _b_heat(b_heat),
        _g_heat(g_heat), _ingestion(ingestion), _inhalation(inhalation) {}

  virtual std::string getElement() const { return _name; }

  virtual std::string getState() { return "state"; }

  virtual int getIsotope() const { return _isotope; }

  virtual int getZAI() const { return _zai; }

  // The half life (s)
  virtual double getHalfLife() const { return _hl; }
  // The number of atoms
  virtual double getAtoms() const { return _atoms; }
  // The grams (g)
  virtual double getGrams() const { return _grams; }
  // The activity (Bq)
  virtual double getActivity() const { return _activity; }
  // The alpha fraction of the activity (Bq)
  virtual double getAlphaActivity() const { return _a_act; }
  // The beta fraction of the activity (Bq)
  virtual double getBetaActivity() const { return _b_act; }
  // The gamma fraction of the activity (Bq)
  virtual double getGammaActivity() const { return _g_act; }
  // The total heat (kW)
  virtual double getTotalHeat() const { return _total_heat; }
  // The alpha heat (kW)
  virtual double getAlphaHeat() const { return _a_heat; }
  // The beta heat (kW)
  virtual double getBetaHeat() const { return _b_heat; }
  // The gamma heat (kW)
  virtual double getGammaHeat() const { return _g_heat; }
  // The dose rate (Sv/hr)
  virtual double getDoseRate() const { return _dose; }
  // The ingestion (Sv)
  virtual double getIngestion() const { return _ingestion; }
  // The inhalation (Sv)
  virtual double getInhalation() const { return _inhalation; }

protected:
  std::string _name;
  std::string _state;
  int _zai;
  int _isotope;
  double _hl, _atoms, _grams, _activity, _a_act, _b_act, _g_act, _total_heat,
      _a_heat, _b_heat, _g_heat, _dose, _ingestion, _inhalation;
};

// IFispactInputData
class IMockFispactInputData : public IFispactInputDataBase {

public:
  IMockFispactInputData() : IFispactInputDataBase() {}

  virtual void setFlux(const std::vector<double> &flux_energy_groups,
                       const std::vector<double> &flux) {}

  virtual void setFluxWallLoading(double wall_loading) {}

  virtual void setExcludeXrays(bool enable) {}

  virtual void setFluxName(std::string flux_name) {}

  virtual void setDensity(double density) {}

  virtual void setMassTotal(double total_mass) {}

  virtual void setMass(const std::vector<int> &atomicnumbers,
                       const std::vector<double> &percentages) {}

  virtual void setFuel(const std::vector<int> &zais,
                       const std::vector<double> &values) {}

  virtual void setSchedule(const std::vector<double> &deltatime,
                           const std::vector<double> &fluxamp) {}

  virtual std::pair<std::vector<double>, std::vector<double>> getSchedule() {
    return schedule;
  }

  virtual void setSolverTolerance(double rtol, double atol) {
    this->rtol = rtol;
    this->atol = atol;
  }

  virtual void setAtomsThreshold(double threshold) {}

private:
  double rtol, atol;
  std::pair<std::vector<double>, std::vector<double>> schedule;
};
//~IFispactInputData

// IFispactOutputData
class IMockFispactOutputData : public IFispactOutputDataBase {

public:
  IMockFispactOutputData() : IFispactOutputDataBase() {}

  virtual std::vector<double> getGammaSpectrumBins(int inv_index) {
    std::vector<double> vec = {1, 1, 1};
    return vec;
  };

  virtual std::vector<double> getGammaSpectrumBoundaries(int inv_index) {
    std::vector<double> vec = {1, 1, 1};
    return vec;
  };

  virtual std::vector<std::unique_ptr<FispactOutputNuclideDataBase>>
  getInventoryNuclides(int inventory_index) {

    std::vector<std::unique_ptr<FispactOutputNuclideDataBase>> nuclide_data;
    nuclide_data.reserve(1);

    // nuclide_data.push_back(std::make_unique<IMockFispactOutputNuclideData>("H1",
    // "liquid", 1, 1, 0.4, 0, 100, 1, 1, 1, 1,1,1));

    return nuclide_data;
  }

  virtual std::pair<std::vector<int>, std::vector<double>>
  getSortedInventory(int, inventory_outputs::InventoryOutputsEnum) const {
    std::vector<int> zais = {1, 2, 3};
    std::vector<double> value = {1, 2, 3};

    return std::pair<std::vector<int>, std::vector<double>>(zais, value);
  }

  virtual double getInventoryValue(int,
                                   inventory_outputs::InventoryOutputsEnum) {
    return 1.0;
  }

  virtual int findInventoryIndex(int, int) { return 1; };

  virtual bool findInventoryExists(int inv_index, int zai) { return true; }

private:
};
//~IFispactOutputData

class IMockFispactUtils : public IFispactUtilsBase {

public:
  IMockFispactUtils() : IFispactUtilsBase() {}

  virtual int GetZai(std::string nuclidename) {
    if (nuclidename == "H1") {
      return 4;
    }

    else if (nuclidename == "H2") {
      return 5;
    }

    else if (nuclidename == "He3") {
      return 6;
    }
    return -1;
  }

  virtual int GetAtomicNumberFromElementName(std::string elementname) const {
    return 1;
  }

  virtual std::vector<double>
  GroupConvertByEnergy(const std::vector<double> &inbounds,
                       const std::vector<double> &invals,
                       const std::vector<double> &outbounds) {
    return invals;
  }

  virtual std::vector<double>
  GroupConvertByLethargy(const std::vector<double> &inbounds,
                         const std::vector<double> &invals,
                         const std::vector<double> &outbounds) {
    return invals;
  }

  virtual std::string getNuclideName(int zai) { return "zai"; };

  virtual int GetAtomicNumberFromElementName(std::string elementname) {
    return 1;
  }

  virtual std::vector<double> getNeutronEnergyBounds(size_t n_groups) {
    return {1, 2, 3};
  }

  virtual std::vector<double> getPhotonEnergyBounds(size_t n_groups) {
    return {1, 2, 3};
  };

private:
};

class FispactContextMock : public FispactContextBase {
public:
  FispactContextMock() : FispactContextBase() {

    _i_input_data = std::make_unique<IMockFispactInputData>();
    _i_output_data = std::make_unique<IMockFispactOutputData>();
    _i_utils = std::make_unique<IMockFispactUtils>();
  }

  virtual void globalInitialise() {}

  virtual void globalFinalise() {}

  virtual void process() {}

  virtual size_t getNuclearDataCrossSections() { return 709; }

  virtual void setNuclearData(
      std::unordered_map<std::string, std::string> nuclear_data_paths) {}

private:
};
