#ifdef FIZZY_UNIT_TEST
#include "FispactContextMock.h"
#else
#include "FispactContext.h"
#endif

#include "FispactFactory.h"

std::unique_ptr<FispactContextBase> createFispactContext() {
#ifdef FIZZY_UNIT_TEST
  return std::make_unique<FispactContextMock>();
#else
  return std::make_unique<FispactContext>();
#endif
}
