#pragma once
#include "settings.h"
#include <iostream>

// Get the current date and time in nice format
std::string getCurrentDateTime(std::string);

// Output a log message for a given image being processed
void Logger(std::string, settingsValues);

// Create empty log file for a given image being processed
void Logger_initialize(settingsValues);