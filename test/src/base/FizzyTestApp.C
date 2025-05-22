//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "FizzyTestApp.h"
#include "FizzyApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
FizzyTestApp::validParams()
{
  InputParameters params = FizzyApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

FizzyTestApp::FizzyTestApp(InputParameters parameters) : MooseApp(parameters)
{
  FizzyTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

FizzyTestApp::~FizzyTestApp() {}

void
FizzyTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  FizzyApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"FizzyTestApp"});
    Registry::registerActionsTo(af, {"FizzyTestApp"});
  }
}

void
FizzyTestApp::registerApps()
{
  registerApp(FizzyApp);
  registerApp(FizzyTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
FizzyTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  FizzyTestApp::registerAll(f, af, s);
}
extern "C" void
FizzyTestApp__registerApps()
{
  FizzyTestApp::registerApps();
}
