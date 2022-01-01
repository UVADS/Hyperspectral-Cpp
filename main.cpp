#include "hsi_data_reader.h"
#include "spectral.h"
#include "settings.h"
#include "write_log_file.h"
#include <iostream>
#include <Eigen/Dense>
#include <filesystem>

using namespace hsi;
using hsi::HSIData;
using hsi::HSIDataOptions;
using hsi::HSIDataReader;


int main(int argc, char** argv) {
  // Start timeing
  auto t_start = std::chrono::high_resolution_clock::now();

  // Set parameters
  int verbose = 1;

  // Get the file names for the image and header
  // from user input or use a default
  const std::string hsi_data_path = getDataFname(argc, argv);
  const std::string hsi_header_path = getHeaderFname(argc, argv);

  // Initialize the settings and log file
  settingsValues settings = createSettings(verbose, hsi_data_path);
  Logger_initialize(settings);

  // Set parameters for Eigen computations
  Logger("Number of Eigen Threads: "+std::to_string(Eigen::nbThreads()), settings);
  Eigen::initParallel();

  // Determine the metadata from the header file
  HSIDataOptions data_options(hsi_data_path);
  data_options.ReadHeaderFromFile(hsi_header_path);
  hsi::HSIDataRange data_range = getDataRange(data_options);

  // Read the image data.
  Logger("Reading image.", settings);
  HSIDataReader reader(data_options);  
  reader.ReadData(data_range);
  const hsi::HSIData& hsi_data = reader.GetData();
  
  // Create eigen array of image from the data. 
  Logger("Converting image data to Eigen matrix.", settings);
  ImData Im = HSIData2EigenData(hsi_data, data_options, settings);

  // Computing the covariance (and outputting computation time)
  Logger("Computing statistics using Eigen.", settings);
  ImStats stats = computeAnomalyCleanedImageStats(Im, settings);
  //ImStats stats = computeImageStats(Im, settings);

  // Whiten the image
  Logger("Whitening the image.", settings);
  Im.white2d = computeWhitenedImage(Im, stats, settings);

  // Compute the RX anomaly image
  Logger("Compute the RX anomaly image.", settings);
  Im.RX = computeRXimage(Im, settings);
  // Display the image if desired
  showRXimage(Im, settings);

  // Get the file names for the target library and header
  // from user input or use a default
  const std::string tgt_lib_data_path = getTgtLibDataFname(hsi_data_path);
  const std::string tgt_lib_header_path = getTgtLibHeaderFname(hsi_data_path);
  std::cout << tgt_lib_data_path << std::endl;
  std::cout << tgt_lib_header_path << std::endl;

  
  // Determine the metadata from the header file
  HSIDataOptions tgt_lib_data_options(tgt_lib_data_path);
  tgt_lib_data_options.ReadHeaderFromFile(tgt_lib_header_path);
  hsi::HSIDataRange tgt_lib_data_range = getDataRange(tgt_lib_data_options);

  // Read the image data.
  Logger("Reading Spectral library.", settings);
  HSIDataReader tgt_lib_reader(tgt_lib_data_options);  
  tgt_lib_reader.ReadData(tgt_lib_data_range);
  const hsi::HSIData& tgt_lib_data = tgt_lib_reader.GetData();

  // Create eigen array of library from the data. 
  Logger("Converting library data to Eigen data.", settings);
  LibData lib_tgt_fullres = SpecLibData2EigenData(tgt_lib_data, tgt_lib_data_options);

  // Create eigen array of library from the data. 
  Logger("Resampling library to image wavelengths.", settings);
  LibData lib_tgt = resampleLibrary(lib_tgt_fullres, Im);

  lib_tgt.spectraWhite = computeWhitenedLibrary(lib_tgt, stats, settings);

  // Compute the ACE target detection image
  Eigen::MatrixXf ACE = computeACEimage(Im, lib_tgt, settings);

  // Output results to a subdirectory called \output
  showACEResults(Im, lib_tgt, ACE, settings);
  
  auto t_end = std::chrono::high_resolution_clock::now();
  auto secondsElapsedWhite = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);
  Logger("Total Runtim: "+std::to_string(secondsElapsedWhite.count()/1000.)+"s", settings);
  return 0;
}
