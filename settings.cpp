#include "settings.h"
#include <string>
#include <filesystem>


// Function definition
settingsValues createSettings(int verboseValue, const std::string hsi_header_path) {
  settingsValues settings;
  settings.verbose = verboseValue; 
  settings.loggerPath = std::filesystem::path(hsi_header_path).replace_filename("output").string()+"\\log.txt";
  settings.outDir = std::filesystem::path(hsi_header_path).replace_filename("output").string();
  settings.displayRX = true;
  settings.saveRX = true;
  settings.displayACE = true;
  settings.saveACE = true;

  std::filesystem::create_directory(settings.outDir);

  return settings;
}