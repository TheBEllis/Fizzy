#pragma once

#include "GeneralUserObject.h"


class FISPACTSchedule : public GeneralUserObject
{   
public:
    static InputParameters validParams();

    FISPACTSchedule(const InputParameters & params);

    virtual void initialize() {}
    virtual void finalize() {}
    virtual void execute() {}
  
protected:
  
    /// Nuclide names
    const std::vector<int32_t> & _schedule;
  
    /// List of times corresponding to the binary on or off
    const std::vector<double> & _times;

};