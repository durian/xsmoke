#pragma once

#include <chrono>

#include "XPLMDisplay.h"
#include "XPCWidget.h"
#include "XPStandardWidgets.h"
#include "XPLMDataAccess.h"
#include "XPLMGraphics.h"
#include "XPWidgetUtils.h"

#include "Log.h"

namespace DROP {

  class InfoWindow;
  void closeInfoWindow();
  void closeInfoWindows();
  void showInfoWindow(const std::string&, const std::string&, const std::string&);
  InfoWindow * showInfoWindow(const std::string&, const std::string&, const std::string&, int, int, float tmr = -1.0);
  
}
