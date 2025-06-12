#pragma once

#include "GeneralUserObject.h"

class FISPACTMaterial : public GeneralUserObject
{   
public:
    static InputParameters validParams();

    FISPACTMaterial(const InputParameters & params);

    virtual void initialize() {}
    virtual void finalize() {}
    virtual void execute() {}
  
protected:
  
    // Nuclide names
    std::vector<std::string> _nuclides;
  
    // List of nuclide proportions corresponding to each nuclide in the nuclide list 
    std::vector<double> _nuclide_proportion;

    // Density of the fispact material
    double _density;
};
