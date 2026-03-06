#include "FispactContext.h"
#include "FispactContextBase.h"
#include "FispactFactory.h"
#include <memory>

std::unique_ptr<FispactContextBase> createFispactContext(bool make_mock) {
  if (make_mock) {
    return 0;
  } else {

    return std::make_unique<FispactContext>();
  }
}
