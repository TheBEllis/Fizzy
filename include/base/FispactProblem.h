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
        // Constructor for FispactProblem
        FispactProblem(const InputParameters & params);

        static InputParameters validParams();

        // virtual void initialSetup() override;
        virtual void externalSolve() override;
        // virtual void syncSolutions(ExternalProblem::Direction direction) override;

        
    private:

        static void load_callback(std::string key, std::string path, int i, int t){
            std::cout << "\33[2K\r" << key << ": " << path << " [" << i << "/" << t << "]" << std::flush;
        }
        
        static void process_callback(std::string process_name, int i, int t){
            std::cout << "\33[2K\r [" << i << "/" << t << "] " << process_name << std::flush;
        }

        void setNuclearData(std::string nd_base_path);

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