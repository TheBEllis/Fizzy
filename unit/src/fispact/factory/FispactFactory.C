#include "FispactContext.h"
#include "FispactContextBase.h"
#include "FispactContextMock.h"
#include "FispactFactory.h"
#include <memory>

std::unique_ptr<FispactContextBase> createFispactContext(bool make_mock) {

  return std::make_unique<FispactContext>();
}
