#pragma once

#include "ExternalProblem.h"


#include "fispactnucleardata.hpp"
#include "fispactinputdata.hpp"
#include "fispactoutputdata.hpp"
#include "fispactgroupstructures.hpp"
#include "fispactgroupconvert.hpp"
#include "fispactelementaldata.hpp"
#include "fispactutil.hpp"
 
// Use fp as short for fispact
namespace fp = fispact;

class FispactProblem : public ExternalProblem
{
    public:
        FispactProblem(const InputParameters & params);

        static InputParameters validParams();

        
    private:

        /// Read a neutron flux spectra from a hdf5 file
        std::vector<double> readNeutronFluxFromH5(const std::string& filename);
        
        /// Write output photon flux to HDF5
        void writePhotonFluxToHDF5(const std::string& filename);

        /// Generate a log file name for the fispact logs
        std::string fispactLogName();
        
        // Initialise FISPACT monitor object
        void initFispactMonitor(std::string);
    
        /// FISPACT monitor
        fp::FispactMonitor _fp_monitor;       

        /// FISPACT nuclear data
        fp::NuclearData _fp_nuclear_data;      
};