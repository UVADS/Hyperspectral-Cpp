# Hyperspectral-Cpp
Hyperspectral Image and Spectral Library Reading, Resampling, Anomaly Detection, and Target Detection in C++.

This is a commented code to read hypersptectral files in ENVI format using the reader C++ code (with modification) from from https://github.com/rteammco/hsi-data-reader written by Richard Teammco.  Algorithm doescriptions for whitening and target detection can be found in: 
<blockquote>
William F. Basener, Brian Allen, Kristen Bretney, "Geometry of statistical target detection," J. Appl. Rem. Sens. 11(1) 015012 (21 February 2017) https://doi.org/10.1117/1.JRS.11.015012
</blockquote>

The spectral.cpp file contains the spectral computations.  The main.cpp file reads an image and executes the following operations:
<ul>
  <li>Read a hyperspectral imag in ENVI format</li>
  <li>Compute covariance and related statistics</li>
  <li>Whiten the image</li>
  <li>Create RX anomaly detection image (option: view or save as jpg)</li>
  <li>Read a spectral library from ENVI format</li>
  <li>Resample the spectral library to image wavelengths using image fwhm when provided</li>
  <li>Compute ACE target detection images (option: view or save as jpg)</li>
</ul>

<b>This code is intended to be funtional for target detection, but also useful as commented core code for building further hyperspectral projects in C++, especially for bulk processing applications.</b>

# Running 

The code assumes you provide a hyperspectral image file named <code>AVIRIS</code>.  A spectral library <code>lib_detect_fullresolution.sli</code> of polymer spectra is provided.  If you put all files in a single directory you can run the code by <code>main.exe</code> or if you want to use a different image execute with <code> main image_filename</code> (which will look for a header with the same name as the image with .hdr appended) or <code>main image_filename image_header_filename</code>.  The code will look for a spectral library with the file name <code>lib_detect_fullresolution.sli</code> and header <code>lib_detect_fullresolution.hdr</code> in the same directory as the image.  It will create an output folder with a log file and out jpgs in the same directory as the image.


# Dependencies:

The Eigen C++ library(https://eigen.tuxfamily.org/index.php?title=Main_Page) is used for matrix operaitons because this is an easy to use and easy to install library that is also fast.

The openCV library is used to view and save output as .jpg images because this is a common library.
