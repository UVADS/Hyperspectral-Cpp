#pragma once
#include <string>

struct settingsValues {
  int verbose;
  std::string loggerPath;
  std::string outDir;
  bool displayRX;
  bool saveRX;
  bool displayACE;
  bool saveACE;
};

// Create the settings (verbose value, logger filename)
settingsValues createSettings(int, const std::string);