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

        void setFispactInputData(fp::Monitor& monitor, fp::InputData& ip, fp::OutputData& ip, std::vector<double> neutron_flux);        

        /// Generate a log file name for the fispact logs
        std::string fispactLogName();
        
        // Initialise FISPACT monitor object
        void initFispactMonitor(std::string);

        void readNeutronFluxFromHDF5(std::string filename);

        void read_material_xml_data();
    
        /// FISPACT monitor
        fp::FispactMonitor _fp_monitor;       

        /// FISPACT nuclear data
        fp::NuclearData _fp_nuclear_data;      

        /// FISPACT neutron flux
        std::vector<std::vector<double>> _neutron_flux;

        /// hdf5 filename for neutron flux
        std::string _neutron_flux_filename;
        
        /// path to neutron flux array in hdf5 file
        std::string _neutron_flux_hdf5_path

        ///
        bool materials_from_xml;
        
        /// Filename of xml file to read materials from
        std::string materials_xml_file;

        /// Struct to store Material definitions
        struct MaterialDefinition
        {
            std::string _mat_name;
            double _mat_density;
            std::vector<std::pair<std::string, double>> _m_atomic_composition;
        };

        /// Mappings from material name to material definitions
        std::unorderedmap<std::string, MaterialDefinition>;
        
};