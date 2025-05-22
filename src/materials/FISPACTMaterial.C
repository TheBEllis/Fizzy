#include "FISPACTMaterial.h"

registerMooseObject("fizzyApp", FISPACTMaterial);

InputParameters
FISPACTMaterial::validParams()
{
    InputParameters params = GeneralUserObject::validParams();

    params.addRequiredParam<std::vector<std::string>>("nuclides",
        "a list of nuclides present in the material");

    params.addRequiredParam<std::vector<double>>("nuclide_proportions",
                               "List of percentages dictating the proportion of the nuclides in the material.");
    params.addRequiredParam<double>("density",
    "Total material density in g/cm^3");

    return params;
}

FISPACTMaterial::FISPACTMaterial(const InputParameters & parameters)
    : GeneralUserObject(parameters),
        _schedule(getParam<std::vector<int>>("schedule")),
        _times(getParam<std::vector<double>>("times"))
{
    // Check vectors passed in are the right size
    mooseAssert(_times.size() != _schedule.size() + 1, "The length of schedule must be one larger than the list of times");
}