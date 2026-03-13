
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
class IMockFispactInputData : public IFispactInputDataBase {

public:
  IMockFispactInputData() : IFispactInputDataBase() {}

  virtual void setFlux(const std::vector<double> &flux_energy_groups,
                       const std::vector<double> &flux) {}

  virtual void setFluxWallLoading(double wall_loading) {}

  virtual void setFluxName(std::string flux_name) {}

  virtual void setDensity(double density) {}

  virtual void setMassTotal(double total_mass) {}

  virtual void setMass(const std::vector<int> &atomicnumbers,
                       const std::vector<double> &percentages) {}

  virtual void setFuel(const std::vector<int> &zais,
                       const std::vector<double> &values) {}

  virtual void setSchedule(const std::vector<double> &deltatime,
                           const std::vector<double> &fluxamp) {}

  virtual void setAtomsThreshold(double threshold) {}

private:
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

private:
};
//~IFispactOutputData

class IMockFispactUtils : public IFispactUtilsBase {

public:
  IMockFispactUtils() : IFispactUtilsBase() {}

  virtual int GetZai(std::string nuclidename) const {
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
