#include "hsi_data_reader.h"
#include "spectral.h"
#include "settings.h"
#include "write_log_file.h"
#include <cmath>
#include <Eigen/Dense>
#include <chrono>
#include <filesystem>

#include "opencv2/highgui.hpp"
#include <vector>
#include <iostream>
#include <sstream>
#include <ctype.h>
#include <cstdlib>
#include <unordered_map>

using hsi::HSIData;
using hsi::HSIDataOptions;
using hsi::HSIDataReader;

float normalPDF(float s, float x, float m){
  return (1/(s*2.50663))*std::pow(2.71828,-0.5*std::pow((x-m)/s,2.0));
}

// convert a comma seperate sting list of the form
// [first, second, third, ... , last]
// to a vector of strings
std::vector<std::string> stringCS2Vector(std::string str){
  str.erase(std::remove(str.begin(), str.end(), '}'), str.end());
  str.erase(std::remove(str.begin(), str.end(), '{'), str.end());
  str.erase(std::remove(str.begin(), str.end(), ']'), str.end());
  str.erase(std::remove(str.begin(), str.end(), '['), str.end());
  std::vector<std::string> result;
  std::stringstream ss(str);
  while (ss.good()){
    std::string substr;
    std::getline(ss, substr, ',');
    result.push_back(substr);
  }
  return result;
}


// convert a comma seperate sting list of the form
// {first, second, third, ... , last}
// to a Eigen::VectorXf 
Eigen::VectorXf stringCS2VectorFloats(std::string str, int len){
  Eigen::VectorXf result(len);
  str.erase(std::remove(str.begin(), str.end(), '}'), str.end());
  str.erase(std::remove(str.begin(), str.end(), '{'), str.end());
  str.erase(std::remove(str.begin(), str.end(), ']'), str.end());
  str.erase(std::remove(str.begin(), str.end(), '['), str.end());
  std::stringstream ss(str);
  int i = 0;
  while (ss.good()){
    std::string substr;
    std::getline(ss, substr, ',');
    float num_float = std::stof(substr);
    result(i) = num_float;
    i++;
  }
  return result;
}

RGBwavelengthsStruct computeRGBwavelengths(Eigen::VectorXf wl){
  RGBwavelengthsStruct RGBwavelengths;
  float redDiffWL = abs(wl(0)-650);
  float greenDiffWL = abs(wl(0)-550);
  float blueDiffWL = abs(wl(0)-450);
  for (int band_idx = 0; band_idx < wl.size(); ++band_idx){
    if (abs(wl(band_idx)-650) < redDiffWL){
      RGBwavelengths.bandRed = band_idx;
    }
    if (abs(wl(band_idx)-550) < redDiffWL){
      RGBwavelengths.bandGreen = band_idx;
    }
    if (abs(wl(band_idx)-450) < redDiffWL){
      RGBwavelengths.bandBlue = band_idx;
    }
  }
  return RGBwavelengths;
}

std::string getDataFname(int argc, char** argv){
  // Set paths to image files.
  if (argc == 1) {
    // if there are no input file names then use the default image
    // AVIRS which should be located in the same directory as the executable file.
    return std::filesystem::path(argv[0]).replace_filename("AVIRIS").string();
  }else {
    // input image and header file names
    return argv[1];
  }
}

std::string getHeaderFname(int argc, char** argv){
  // Set paths to image files.
  if (argc == 1) {
    // if there are no input file names then use the default
    return std::filesystem::path(argv[0]).replace_filename("AVIRIS.hdr").string();;
  }else if (argc == 2) {
    // case where file name is image and header is just the image with .hdr added
    return strcat(argv[1],".hdr");
  } else {
    // input image and header file names
    return argv[2];
  }
}

// Get the range of data (lines, columns, bands) to read
// Currently read the whole image, but subsets could be read here
hsi::HSIDataRange getDataRange(hsi::HSIDataOptions data_options){
  // Set range of data we want to read.
  hsi::HSIDataRange data_range;
  data_range.start_row = 0;
  data_range.end_row = data_options.num_data_rows;
  data_range.start_col = 0;
  data_range.end_col = data_options.num_data_cols;
  data_range.start_band = 0;
  data_range.end_band = data_options.num_data_bands;
  return data_range;
}

ImData HSIData2EigenData(const hsi::HSIData& hsi_data, hsi::HSIDataOptions data_options, settingsValues settings){
  // Start timing
  auto t_start_HSIData2EigenData = std::chrono::high_resolution_clock::now();
  
  // Create the structure to hold the image and metadata
  ImData Im;

  Im.rows = hsi_data.num_rows;
  Im.cols = hsi_data.num_cols;
  Im.bands = hsi_data.num_bands;
  Eigen::MatrixXf im2d(Im.rows * Im.cols, Im.bands);
  for (int band = 0; band < Im.bands; ++band) {
    for (int row = 0; row < Im.rows; ++row) {
      for (int col = 0; col < Im.cols; ++col) {
        im2d(row * Im.cols + col, band) = hsi_data.GetValue(row, col, band).value_as_float;
      }
    }
  }
  Im.im2d = im2d;

  // This creates the Eigen::VectorXf Im.wl from the wavelength string
  Im.wl = stringCS2VectorFloats(data_options.wavelength, Im.bands);

  // Im.wlScale is a scale factor on the wavelengths
  // multiplying by this factor puts the wavelengths
  // into units of nanometers
  Im.wlScale = 1;
  if (Im.wl.mean() < 10){
    Im.wlScale = 1000;
  }

  // Compute the indices for the wavelengths associated with 
  // Reg, Green, Blue colors (650nm, 550nm, 450nm)
  RGBwavelengthsStruct RGBwavelengths = computeRGBwavelengths(Im.wl);
  int bandRed = RGBwavelengths.bandRed;
  int bandGreen = RGBwavelengths.bandGreen;
  int bandBlue = RGBwavelengths.bandBlue;
  
  // This creates the Eigen::VectorXf Im.fwhm from the fwhm string
  Im.fwhm = stringCS2VectorFloats(data_options.fwhm, Im.bands);

  auto t_end_HSIData2EigenData = std::chrono::high_resolution_clock::now();
  auto secondsElapsedWhite = std::chrono::duration_cast<std::chrono::milliseconds>(t_end_HSIData2EigenData - t_start_HSIData2EigenData);
  Logger("Converting HSIDATA to Eigen Matrices completed in : "+std::to_string(secondsElapsedWhite.count()/1000.)+"s", settings);
  
  std::cout << "some data from the image: \n" << im2d(Eigen::seq(0,5),Eigen::seq(0,5)) << std::endl;


  return Im;
}

void showRGBimage(ImData Im, settingsValues settings){
  //  create grayscale image
  Eigen::MatrixXf RX2d = Im.RX.reshaped(Im.rows, Im.cols);
  const cv::Size RX_image_size(RX2d.cols(), RX2d.rows());
  float min_value = 0;
  float max_value = 0;
  cv::Mat RX_image(RX_image_size, CV_64FC1);
  for (int row = 0; row < RX2d.rows(); ++row) {
    for (int col = 0; col < RX2d.cols(); ++col) {
      const float pixel_value = RX2d(row, col);
      RX_image.at<double>(row, col) = pixel_value;
      min_value = std::min(min_value, pixel_value);
      max_value = std::max(max_value, pixel_value);
    }
  }
  RX_image = 255*(RX_image-min_value)/(max_value-min_value);
  cv::namedWindow("test", cv::WINDOW_AUTOSIZE);
  cv::imshow("test", RX_image);
  cv::waitKey(0);
}

// Computes statistics for an image
ImStats computeImageStats(ImData Im, settingsValues settings) {
  // Create the structure to hold all the stats
  ImStats stats; 

  // Compute the image mean
  stats.mean = Im.im2d.colwise().mean();

  // Compute the covariance matrix
  auto t_start_cov_Eigin = std::chrono::high_resolution_clock::now();
  Eigen::MatrixXd centered = (Im.im2d.rowwise() - stats.mean.transpose()).template cast<double>();
  Eigen::MatrixXd cov = (centered.transpose() * centered) / float(Im.im2d.rows() - 1);
  auto t_end_cov_Eigin = std::chrono::high_resolution_clock::now();
  auto secondsElapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t_end_cov_Eigin - t_start_cov_Eigin);
  Logger("Covariance computation completed in "+std::to_string(secondsElapsed.count()/1000.)+"s", 
    settings);

  // Compute the eigenvalues and eigenvectors
  auto t_start_evalsEvecs_Eigin = std::chrono::high_resolution_clock::now();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(cov);
  auto t_end_evalsEvecs_Eigin = std::chrono::high_resolution_clock::now();
  Eigen::MatrixXd evals;
  Eigen::MatrixXd evecs;
  if (eigensolver.info() == Eigen::Success){
    auto secondsElapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t_end_evalsEvecs_Eigin - t_start_evalsEvecs_Eigin);
    Logger("Eigenvalue and eignvector computation completed in "+std::to_string(secondsElapsed.count()/1000.)+"s", 
      settings);
    evals = eigensolver.eigenvalues();
    evecs = eigensolver.eigenvectors();
  }
  else {
    Logger("WARNING: eigenvalues and eigenvectors failed to compute.", settings);
  };

  
  std::cout << "Eigenvalues:\n";
  for (int i = 0; i < 3; ++i) {
    std::cout << evals(i) << " | " << evals(i)/evals(evals.rows()-1) << std::endl;
  }
  std::cout << "...\n";
  for (int i = evals.rows()-3; i < evals.rows(); ++i) {
    std::cout << evals(i) << " | " << evals(i)/evals(evals.rows()-1) << std::endl;
  }

  // Computing the whitening matrix
  auto t_start_W = std::chrono::high_resolution_clock::now();
  Eigen::MatrixXd D = Eigen::MatrixXd::Identity(evals.rows(),evals.rows());
  double maxEval = evals.maxCoeff();
  for (int i = 0; i < evals.rows(); ++i) {
    // Optional regularization on eigenvalues
    //stats.evals(i) = std::max(stats.evals(i),maxEval*std::pow(10,-8));
    D(i,i) = 1/sqrt(evals(i));
  }
  Eigen::MatrixXd W =  (evecs*D);
  std::cout << "some data from the whitening matrix: \n" << W(Eigen::seq(0,5),Eigen::seq(0,5)) << std::endl;
  
  auto t_end_W = std::chrono::high_resolution_clock::now();
  auto secondsElapsedW = std::chrono::duration_cast<std::chrono::milliseconds>(t_end_W - t_start_W);
  Logger("Whitening matrix computation completed in "+std::to_string(secondsElapsedW.count()/1000.)+"s", 
    settings);

  // Recast as float
  stats.cov = cov.template cast<float>();
  stats.evals = evals.template cast<float>();
  stats.evecs = evecs.template cast<float>();
  stats.W = W.template cast<float>();

  Logger("im2d Matrix is "+std::to_string(Im.im2d.rows())+"x"+std::to_string(Im.im2d.cols()),
    settings);
  Logger("Covariance is "+std::to_string(stats.cov.rows())+"x"+std::to_string(stats.cov.cols()),
    settings);  
  Logger("Whitening Matrix is "+std::to_string(stats.W.rows())+"x"+std::to_string(stats.W.cols()),
    settings);

  // VALIDATION: If we want to validate the Whitening by comparison to the inverse
  //Eigen::MatrixXd Diff = cov.inverse() - (W)*(W.transpose());
  //std::cout << "some data from the Diff matrix: \n" << Diff(Eigen::seq(0,5),Eigen::seq(0,5)) << std::endl;

  return stats;
}

// Whitens an image using provided statistics
Eigen::MatrixXf computeWhitenedImage(ImData Im, ImStats stats, settingsValues settings) {

  auto t_start_White = std::chrono::high_resolution_clock::now();
  // VALIDATION: If we want to validat the data be viewing the image mean
  //std::cout << "image mean: \n" << stats.mean;
  Eigen::MatrixXf centered = Im.im2d.rowwise() - stats.mean.transpose();
  Eigen::MatrixXf whitened = centered*stats.W;
  auto t_end_White = std::chrono::high_resolution_clock::now();
  auto secondsElapsedWhite = std::chrono::duration_cast<std::chrono::milliseconds>(t_end_White - t_start_White);
  Logger("Whitening image completed in "+std::to_string(secondsElapsedWhite.count()/1000.)+"s", 
    settings); 
  Logger("Whitening image is "+std::to_string(whitened.rows())+"x"+std::to_string(whitened.cols()),
    settings);

  return whitened;
}

// Computes the RX anomaly detetion image from the whitened image
Eigen::MatrixXf computeRXimage(ImData Im, settingsValues settings){
  
  auto t_start_rx = std::chrono::high_resolution_clock::now();
  Eigen::MatrixXf  RX = Im.white2d.rowwise().norm();
  Logger("RX image has size: "+std::to_string(RX.size()), settings);
  auto t_end_rx = std::chrono::high_resolution_clock::now();
  auto secondsElapsedRX = std::chrono::duration_cast<std::chrono::milliseconds>(t_end_rx - t_start_rx);
  Logger("Creating RX image completed in "+std::to_string(secondsElapsedRX.count()/1000.)+"s", settings);

  return RX;
}

// Displays the RX anomaly detection image
void showRXimage(ImData Im, settingsValues settings){

  if (settings.saveRX || settings.displayRX){
    //  create grayscale image
    Eigen::MatrixXf RX2d = Im.RX.reshaped(Im.rows, Im.cols);
    const cv::Size RX_image_size(RX2d.cols(), RX2d.rows());
    float min_value = 0;
    float max_value = 0;
    cv::Mat RX_image(RX_image_size, CV_64FC1);
    for (int row = 0; row < RX2d.rows(); ++row) {
      for (int col = 0; col < RX2d.cols(); ++col) {
        const float pixel_value = RX2d(row, col);
        RX_image.at<double>(row, col) = pixel_value;
        min_value = std::min(min_value, pixel_value);
        max_value = std::max(max_value, pixel_value);
      }
    }
    RX_image = 5*(RX_image-min_value)/(max_value-min_value);
    
    if (settings.displayRX){
      cv::namedWindow("test", cv::WINDOW_AUTOSIZE);
      cv::imshow("test", RX_image);
      cv::waitKey(0);
    }

    if (settings.saveRX){
      bool check = cv::imwrite(settings.outDir+"\\RX.jpg", 255*RX_image);
      if (check == false) {
        Logger("WARNING Saving RX image failed.", settings);
      } else {
        Logger("RX Image saves as "+settings.outDir+"\\RX.jpg", settings);
      }
    }
  }
}

const std::string getTgtLibDataFname(const std::string hsi_data_path){
  std::filesystem::path p(hsi_data_path);
  return p.parent_path().string()+"\\lib_detect_fullresolution.sli";
}

const std::string getTgtLibHeaderFname(const std::string hsi_data_path){
  std::filesystem::path p(hsi_data_path);
  return p.parent_path().string()+"\\lib_detect_fullresolution.hdr";
}

// Converts data and metadata for a spectral library into 
// a LibData structure using Eigen variable types.
LibData SpecLibData2EigenData(const hsi::HSIData& tgt_lib_data, hsi::HSIDataOptions tgt_lib_data_options){
  // Create the structure to hold the image and metadata
  LibData lib_tgt;
  lib_tgt.nSpectra = tgt_lib_data.num_cols;
  lib_tgt.bands = tgt_lib_data.num_rows;

  // This creates the vector Im.wl from the wavelength string
  lib_tgt.wl = stringCS2VectorFloats(tgt_lib_data_options.wavelength, lib_tgt.bands);

  // This creates the vector lib_tgt.specNames from the specNames string
  lib_tgt.specNames = stringCS2Vector(tgt_lib_data_options.specNames);

  // Make a matrix that will hold the spectra
  Eigen::MatrixXf spectra(lib_tgt.nSpectra, lib_tgt.bands);
  for (int i = 0; i < lib_tgt.nSpectra; ++i) {
    for (int j = 0; j < lib_tgt.bands; ++j) {
      spectra(i,j) = 10.0*tgt_lib_data.GetValueLib(j, i, 0).value_as_float;
    }
  }
  lib_tgt.spectra = spectra;
  std::cout << "some data from the library file spectra: \n" << spectra(Eigen::seq(0,5),Eigen::seq(0,5)) << std::endl;

  return lib_tgt;
}

LibData resampleLibrary(LibData lib_tgt_fullres, ImData Im) {
  // Create the structure to hold the image and metadata
  LibData lib_tgt;
  lib_tgt.specNames = lib_tgt_fullres.specNames;
  lib_tgt.nSpectra = lib_tgt_fullres.nSpectra;
  lib_tgt.bands = Im.bands;
  lib_tgt.wl = Im.wl;
  Eigen::MatrixXf spectra = Eigen::MatrixXf::Zero(lib_tgt_fullres.nSpectra, Im.bands);

  Eigen::VectorXf sigma(Im.bands);
  for (int i = 0; i < Im.bands; ++i) {
    sigma(i) = Im.fwhm(i)/2.355;
  }

  float sum = 0;
  float weight = 0;
  for (int im_bnd_idx = 0; im_bnd_idx < Im.bands; ++im_bnd_idx) {  //  looping over bands in the image
    for (int lb_bnd_idx = 0; lb_bnd_idx < lib_tgt_fullres.bands; ++lb_bnd_idx) { // looping over all bands in the library 
      //  compute the difference between the image (target) wavelength and library wavelength so that we 
      //  only compute contributions of library bands with wavelength that are withing 3 sigma of the image wavelength
      float diff = (lib_tgt_fullres.wl[lb_bnd_idx] - Im.wl[im_bnd_idx]);
      float standard_deviations_diff = diff/sigma(im_bnd_idx);
      
      /*  VALIDATION OUTPUT:  uncomment this to see how the resampling computation works
      std::cout << lib_tgt_fullres.wl[lb_bnd_idx] << " - " << Im.wl[im_bnd_idx] << std::endl;
      std::cout << diff << std::endl;
      std::cout << standard_deviations_diff << std::endl;
      std::cout << "Weight: " << weight << std::endl;
      std::cout << "Sum: " << sum << std::endl;
      */

      if (standard_deviations_diff > -3){
        if (standard_deviations_diff > 3){
          // stop the  iteration through library bands becasue the library band wavelength has
          // passed all wavelengths within 3 standard deviations of the image wavelength
          lb_bnd_idx = lib_tgt_fullres.bands;
        } else {
          weight = normalPDF(sigma(im_bnd_idx), 
                            Im.wl[im_bnd_idx],
                            lib_tgt_fullres.wl[lb_bnd_idx]);
          for (int spec_idx = 0; spec_idx < lib_tgt.nSpectra; ++spec_idx) { // looping over all spectra in the library
            spectra(spec_idx, im_bnd_idx) = spectra(spec_idx, im_bnd_idx) + 
                                        lib_tgt_fullres.spectra(spec_idx, lb_bnd_idx)*
                                        weight;
            sum = sum + weight;
          }
        }
      }
    }
    for (int spec_idx = 0; spec_idx < lib_tgt.nSpectra; ++spec_idx) { // looping over all spectra in the library
      spectra(spec_idx, im_bnd_idx) = spectra(spec_idx, im_bnd_idx)/sum;
    }
    sum = 0;    
  }
  lib_tgt.spectra = spectra;

  /* VALIDATION: If we want to validate data by checking the first resampled spectrum from the library
  std::cout << "First spectrum:\n";
  for (int idx = 0; idx < Im.bands; ++idx) {
    std::cout << spectra(0,idx) << std::endl;
  }
  */
  
  std::cout << "some data from the resampled spectra: \n";
  std::cout << spectra(Eigen::seq(0,5),Eigen::seq(0,5)) << std::endl;

  return lib_tgt;
}

Eigen::MatrixXf computeWhitenedLibrary(LibData lib_tgt, ImStats stats, settingsValues settings) {

  auto t_start_White = std::chrono::high_resolution_clock::now();
  Eigen::MatrixXf centered = lib_tgt.spectra.rowwise() - stats.mean.transpose();
  std::cout << "some data from the centered spectra: \n" << centered(Eigen::seq(0,5),Eigen::seq(0,5)) << std::endl;
  Eigen::MatrixXf spectraWhite = centered*stats.W;
  std::cout << "some data from the whitened spectra: \n" << spectraWhite(Eigen::seq(0,5),Eigen::seq(0,5)) << std::endl;
  auto t_end_White = std::chrono::high_resolution_clock::now();
  auto secondsElapsedWhite = std::chrono::duration_cast<std::chrono::milliseconds>(t_end_White - t_start_White);
  Logger("Whitening library completed in "+std::to_string(secondsElapsedWhite.count()/1000.)+"s", 
    settings); 
  Logger("Whitening library is "+std::to_string(spectraWhite.rows())+"x"+std::to_string(spectraWhite.cols()),
    settings);

  return spectraWhite;
}


Eigen::MatrixXf computeACEimage(ImData Im, LibData lib_tgt, settingsValues settings){
  
  auto t_start_ACE = std::chrono::high_resolution_clock::now();
  std::cout << "some data from the whitened image: \n" << Im.white2d(Eigen::seq(0,5),Eigen::seq(0,5)) << std::endl;
  std::cout << "some data from the whitened library: \n" << lib_tgt.spectraWhite(Eigen::seq(0,5),Eigen::seq(0,5)) << std::endl;
  Eigen::MatrixXf  ACE_num = Im.white2d * lib_tgt.spectraWhite.transpose();
  std::cout << "some data from the ACE numerator: \n" << ACE_num(Eigen::seq(0,5),Eigen::seq(0,5)) << std::endl;
  Eigen::MatrixXf  ACE_denom = (Im.RX.replicate(1,lib_tgt.nSpectra).array() * lib_tgt.spectraWhite.transpose().colwise().norm().replicate(Im.rows*Im.cols, 1).array()).matrix();
  std::cout << "some data from the ACE denominator: \n" << ACE_denom(Eigen::seq(0,5),Eigen::seq(0,5)) << std::endl;
  Eigen::MatrixXf  ACE = (ACE_num.array() / ACE_denom.array()).matrix();
  std::cout << "some data from the ACE image: \n" << ACE(Eigen::seq(0,5),Eigen::seq(0,5)) << std::endl;
  Logger("ACE image is "+std::to_string(ACE_num.rows())+"x"+std::to_string(ACE_num.cols()),
    settings);
  auto t_end_ACE = std::chrono::high_resolution_clock::now();
  auto secondsElapsedACE = std::chrono::duration_cast<std::chrono::milliseconds>(t_end_ACE - t_start_ACE);
  Logger("Creating ACE image completed in "+std::to_string(secondsElapsedACE.count()/1000.)+"s", settings);

  return ACE;
}

void showACEResults(ImData Im, LibData lib_tgt, Eigen::MatrixXf ACE, settingsValues settings){
    
  if (settings.saveACE || settings.displayACE){
    // Create the openCV ACE image
    const cv::Size ACEimSize(Im.cols, Im.rows);
    cv::Mat ACEim(ACEimSize, CV_64FC1);

    for (int spec_idx = 0; spec_idx < lib_tgt.nSpectra; ++spec_idx) {
      Eigen::MatrixXf ACEslice2d = ACE(Eigen::all, spec_idx);
      Eigen::MatrixXf ACE2d = ACEslice2d.reshaped(Im.rows, Im.cols);
      //  create grayscale image
      float min_value = 100;
      float max_value = -100;    
      for (int row = 0; row < Im.rows; ++row) {
        for (int col = 0; col < Im.cols; ++col) {
          float pixel_value = ACE2d(row, col);
          //pixel_value = ((pixel_value>0)-(pixel_value<0))*std::sqrt(std::abs(pixel_value));
          float min_thresh = 0.25;
          pixel_value = std::max(min_thresh,pixel_value);
          ACEim.at<double>(row, col) = pixel_value;
          min_value = std::min(min_value, pixel_value);
          max_value = std::max(max_value, pixel_value);
        }
      }
      ACEim = (ACEim-min_value)/(max_value-min_value);
      std::cout << lib_tgt.specNames[spec_idx] << std::endl;
      std::cout << "Min Value: " << min_value << std::endl;
      std::cout << "Max Value: " << max_value << std::endl;
 
      if (settings.displayACE){
        cv::namedWindow(lib_tgt.specNames[spec_idx], cv::WINDOW_AUTOSIZE);
        cv::imshow(lib_tgt.specNames[spec_idx], ACEim);
        cv::waitKey(0);
      }
      if (settings.saveACE){
        std::string specName = lib_tgt.specNames[spec_idx];
        specName.erase(std::remove(specName.begin(), specName.end(), ' '), specName.end());
        bool check = cv::imwrite(settings.outDir+"\\ACE_"+specName+".jpg", 255*ACEim);
        if (check == false) {
          Logger("WARNING Saving ACE image failed.", settings);
        } else {
          Logger("ACE Image saved as "+settings.outDir+"\\ACE_"+specName+".jpg", settings);
        }
      }
    }
  }

}
  
// Compute a covarnaince matrix and other stats for the input image
// using a 3 step process: compute stats and RX anomaly detection,
// remove the anomalies, and recompute stats from non-anomalous pixels.
// Input image matrix is pixels x bands
ImStats computeAnomalyCleanedImageStats(ImData Im, settingsValues settings) {
  
  // build a subseted compy of the image structure containing
  // only about 10,000 pixel spectra for quick approximate
  // first covariance computation
  ImData Im_sampled;
  Im_sampled.im2d = Im.im2d(Eigen::seq(0,Im.im2d.rows(),std::floor(Im.im2d.rows()/10000.)), Eigen::all);
  Im_sampled.rows = Im.rows;
  Im_sampled.cols = Im.cols;

  // compute the covariance for the subset of pixel spectra
  Logger("Computing statistics using Eigen (first pass)", settings);
  ImStats stats_sampled = computeImageStats(Im_sampled, settings);

  // Whiten the full image using the subset-computed covariance
  Logger("Whitening the image. (first pass)", settings);
  Im_sampled.white2d = computeWhitenedImage(Im, stats_sampled, settings);

  // Compute the RX anomaly image using the subset-computed covariance
  // This will be used to remove anomalies for the 
  // primary covariance computation
  Logger("Compute the RX anomaly image. (first pass)", settings);
  Eigen::MatrixXf RX = computeRXimage(Im_sampled, settings);
  
  // compute a threshold for anomaly removeal
  // removing the most anomalous 5% of pixels
  auto t_start_rx_pctile = std::chrono::high_resolution_clock::now();  
  std::vector<float> RXv;
  RXv.resize(RX.size());
  Eigen::VectorXf::Map(&RXv[0], RX.size()) = RX;  
  std::sort(RXv.begin(), RXv.end());
  int idx_95pct = std::round(RX.size()*0.95);
  for (int i = 0; i < 100; ++i) { 
    std::cout << RXv[i] << std::endl;
  }
  float thresh_95pct = RXv[idx_95pct];
  std::cout << "Thresh: " << thresh_95pct << std::endl;
  std::cout << "Max: " << RX.maxCoeff() << std::endl;
  auto t_end_rx_pctile = std::chrono::high_resolution_clock::now();
  auto secondsElapsedRX_pctile = std::chrono::duration_cast<std::chrono::milliseconds>(t_end_rx_pctile - t_start_rx_pctile);
  Logger("Computing RX threshold completed in "+std::to_string(secondsElapsedRX_pctile.count()/1000.)+"s", settings);
  
  auto t_subsample_start = std::chrono::high_resolution_clock::now();  
  Eigen::VectorXi is_selected = (RX.array() < thresh_95pct).cast<int>(); 
  Eigen::VectorXi is_selectedindices(is_selected.sum());
  int idx = 0;
  for (int i = 0; i < Im.rows; ++i) { 
    if (is_selected(i) == 1){
      is_selectedindices(idx) = i;
      idx++;
    }
  }

  Im_sampled.im2d = Im.im2d(is_selectedindices,Eigen::all);
  std::cout << Im_sampled.im2d.rows() << std::endl;
  std::cout << Im_sampled.im2d.cols() << std::endl;
  auto t_subsample_end = std::chrono::high_resolution_clock::now();
  auto secondsElapsed_subsample = std::chrono::duration_cast<std::chrono::milliseconds>(t_subsample_end - t_subsample_start);
  Logger("time for subset computation: "+std::to_string(secondsElapsed_subsample.count()/1000.)+"s", settings);

  // HAVING TROUBLE SUBSETTING THE IMAGE TO COMPUTE STATS FROM THE NON-ANOMALUS PIXELS
  // I AM MULTIPLYING ANOMALIES BY ZER FOR COV COMPUTATION BUT ITS NOT WORKING!


  // Now we compute the stats for the data with anomalies removed
  // Create the structure to hold all the stats
  ImStats stats; 

  // Compute the image mean
  stats.mean = Im_sampled.im2d.colwise().mean();

  // Compute the covariance matrix
  auto t_start_cov_Eigin = std::chrono::high_resolution_clock::now();
  Eigen::MatrixXd centered = (Im_sampled.im2d.rowwise() - stats.mean.transpose()).template cast<double>();
  Eigen::MatrixXd cov = (centered.transpose() * centered) / float(is_selected.sum() - 1);
  auto t_end_cov_Eigin = std::chrono::high_resolution_clock::now();
  auto secondsElapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t_end_cov_Eigin - t_start_cov_Eigin);
  Logger("Covariance computation completed in "+std::to_string(secondsElapsed.count()/1000.)+"s", 
    settings);

  // Compute the eigenvalues and eigenvectors
  auto t_start_evalsEvecs_Eigin = std::chrono::high_resolution_clock::now();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(cov);
  auto t_end_evalsEvecs_Eigin = std::chrono::high_resolution_clock::now();
  Eigen::MatrixXd evals;
  Eigen::MatrixXd evecs;
  if (eigensolver.info() == Eigen::Success){
    auto secondsElapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t_end_evalsEvecs_Eigin - t_start_evalsEvecs_Eigin);
    Logger("Eigenvalue and eignvector computation completed in "+std::to_string(secondsElapsed.count()/1000.)+"s", 
      settings);
    evals = eigensolver.eigenvalues();
    evecs = eigensolver.eigenvectors();
  }
  else {
    Logger("WARNING: eigenvalues and eigenvectors failed to compute.", settings);
  };

  
  std::cout << "Eigenvalues:\n";
  for (int i = 0; i < 3; ++i) {
    std::cout << evals(i) << " | " << evals(i)/evals(evals.rows()-1) << std::endl;
  }
  std::cout << "...\n";
  for (int i = evals.rows()-3; i < evals.rows(); ++i) {
    std::cout << evals(i) << " | " << evals(i)/evals(evals.rows()-1) << std::endl;
  }

  // Computing the whitening matrix
  auto t_start_W = std::chrono::high_resolution_clock::now();
  Eigen::MatrixXd D = Eigen::MatrixXd::Identity(evals.rows(),evals.rows());
  double maxEval = evals.maxCoeff();
  for (int i = 0; i < evals.rows(); ++i) {
    // Optional regularization on eigenvalues
    //evals(i) = std::max(evals(i),maxEval*std::pow(10,-6));
    D(i,i) = 1/sqrt(evals(i));
  }
  Eigen::MatrixXd W =  (evecs*D);
  std::cout << "some data from the whitening matrix: \n" << W(Eigen::seq(0,8),Eigen::seq(0,8)) << std::endl;
  
  auto t_end_W = std::chrono::high_resolution_clock::now();
  auto secondsElapsedW = std::chrono::duration_cast<std::chrono::milliseconds>(t_end_W - t_start_W);
  Logger("Whitening matrix computation completed in "+std::to_string(secondsElapsedW.count()/1000.)+"s", 
    settings);

  // Recast as float
  stats.cov = cov.template cast<float>();
  stats.evals = evals.template cast<float>();
  stats.evecs = evecs.template cast<float>();
  stats.W = W.template cast<float>();

  Logger("im2d Matrix is "+std::to_string(Im.im2d.rows())+"x"+std::to_string(Im.im2d.cols()),
    settings);
  Logger("Covariance is "+std::to_string(stats.cov.rows())+"x"+std::to_string(stats.cov.cols()),
    settings);  
  Logger("Whitening Matrix is "+std::to_string(stats.W.rows())+"x"+std::to_string(stats.W.cols()),
    settings);

  // VALIDATION: If we want to validate the Whitening by comparison to the inverse
  //Eigen::MatrixXd Diff = cov.inverse() - (W)*(W.transpose());
  //std::cout << "some data from the Diff matrix: \n" << Diff(Eigen::seq(0,5),Eigen::seq(0,5)) << std::endl;

  return stats;
}
