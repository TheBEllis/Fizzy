#pragma once
#include "FispactContextBase.h"
#include <memory>

std::unique_ptr<FispactContextBase> createFispactContext(bool make_mock);
