#pragma once
#include "settings.h"
#include <Eigen/Dense>
#include <unordered_map>

// Structure to hold the image data and metadata
struct ImData {
  Eigen::MatrixXf im2d;
  Eigen::MatrixXf white2d;
  Eigen::MatrixXf RX;
  Eigen::VectorXf fwhm;
  Eigen::VectorXf wl;
  int wlScale; // 1 if wl is in nm, 1000 if wl is in micrometers
  int rows;
  int cols;
  int bands;
  int bandRed = 0;
  int bandGreen = 0;
  int bandBlue = 0;
};

struct RGBwavelengthsStruct {
  int bandRed = 0;
  int bandGreen = 0;
  int bandBlue = 0;
};

// Structure to hold the spectral library data and metadata
struct LibData {
  Eigen::MatrixXf spectra;
  Eigen::MatrixXf spectraWhite;
  Eigen::VectorXf wl;
  std::vector<std::string> specNames;
  int nSpectra;
  int bands;
};

// Structure to hold the image statistics
struct ImStats {
  Eigen::MatrixXf cov;
  Eigen::VectorXf mean;
  Eigen::MatrixXf evals;
  Eigen::MatrixXf evecs;
  Eigen::MatrixXf W;
};

// Compute the pdf for the normal distribution
float normalPDF(float, float, float);

// convert a comma seperate sting list of the form
// {first, second, third, ... , last}
// to a Eigen::VectorXf 
std::vector<std::string> stringCS2Vector(std::string);

// convert a comma seperate sting list of the form
// {first, second, third, ... , last}
// to a vector of floats
Eigen::VectorXf stringCS2VectorFloats(std::string, int);

// Compute the incides for the wavelengths for the RGB colors
RGBwavelengthsStruct computeRGBwavelengths(Eigen::VectorXf);

// Get the name of the data file based on user input
std::string getDataFname(int, char**);

// Get the name of the header file based on user input
std::string getHeaderFname(int, char**);

// Get the range of data (lines, columns, bands) to read
// Currently read the whole image, but subsets could be read here
hsi::HSIDataRange getDataRange(hsi::HSIDataOptions);

// Compute a covarnaince matrix for the input image
// Input image matrix is pixels x bands
ImData HSIData2EigenData(const hsi::HSIData&, hsi::HSIDataOptions, settingsValues);

// Show the RX image
void showRGBimage(ImData, settingsValues);

// Compute a covarnaince matrix and other stats for the input image
// Input image matrix is pixels x bands
ImStats computeImageStats(ImData, settingsValues);

// Compute the whitened image
// The whitened image is the mean-subtracted image
// projected on the eignevectors scaled by eigenvalues
Eigen::MatrixXf computeWhitenedImage(ImData, ImStats, settingsValues);

// Compute the RX anomaly image
Eigen::MatrixXf computeRXimage(ImData, settingsValues);

// Show the RX image
void showRXimage(ImData, settingsValues);

// Get the name of the spectral target library
const std::string getTgtLibDataFname(const std::string);

// Get the name of the spectral target library header
const std::string getTgtLibHeaderFname(const std::string);

// Compute a structure to hold the spectral library information
LibData SpecLibData2EigenData(const hsi::HSIData&, hsi::HSIDataOptions);

// Resample library to a given set os wavelengths
LibData resampleLibrary(LibData, ImData);

// Compute the whitened spectral library
// The whitened library spectra are the mean-subtracted spectra
// projected on the eignevectors scaled by eigenvalues
Eigen::MatrixXf computeWhitenedLibrary(LibData, ImStats, settingsValues);

// Compute the ACE target detection image
Eigen::MatrixXf computeACEimage(ImData, LibData, settingsValues);

//  Outputs results to a folder
void showACEResults(ImData, LibData, Eigen::MatrixXf, settingsValues);

// Compute a covarnaince matrix and other stats for the input image
// using a 3 step process: compute stats and RX anomaly detection,
// remove the anomalies, and recompute stats from non-anomalous pixels.
// Input image matrix is pixels x bands
ImStats computeAnomalyCleanedImageStats(ImData, settingsValues);