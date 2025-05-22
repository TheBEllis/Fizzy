#include "FizzyApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
FizzyApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

FizzyApp::FizzyApp(InputParameters parameters) : MooseApp(parameters)
{
  FizzyApp::registerAll(_factory, _action_factory, _syntax);
}

FizzyApp::~FizzyApp() {}

void
FizzyApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<FizzyApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"FizzyApp"});
  Registry::registerActionsTo(af, {"FizzyApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
FizzyApp::registerApps()
{
  registerApp(FizzyApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
FizzyApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  FizzyApp::registerAll(f, af, s);
}
extern "C" void
FizzyApp__registerApps()
{
  FizzyApp::registerApps();
}
