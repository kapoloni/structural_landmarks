/**
 * @file   __main.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup filterbank
 * @ingroup    filterbank
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

// #define BIP_DEBUG_MODE
// #define BIP_VERBOSE_MODE

#include <cstdio>
#include <cstdlib>
#include <string>
#include "CmdLine.h"
#include "triple.hpp"
#include "loggaborfilterbank.hpp"

// Default argument values.
#define DEFAULT_IMAGE_SIZE_X    "128"
#define DEFAULT_IMAGE_SIZE_Y    "128"
#define DEFAULT_IMAGE_SIZE_Z    "128"
#define DEFAULT_NUM_SCALES      "4"
#define DEFAULT_NUM_AZIMUTHS    "6"
#define DEFAULT_NUM_ELEVATIONS  "3"
#define DEFAULT_MAX_FREQUENCY   "0.3"
#define DEFAULT_MULT_FACTOR     "2.1"
#define DEFAULT_FREQUENCY_RATIO "0.55"
#define DEFAULT_ANGULAR_RATIO   "1.2"
#define DEFAULT_LOWPASS_ORDER   "15"
#define DEFAULT_LOWPASS_CUTOFF  "0.45"


void show_usage(char *argv[])
{
    printf("---------------------------------------------------------------------------------\n"
           "                              Log-Gabor Filter Bank                              \n"
           "---------------------------------------------------------------------------------\n"
           " Usage: %s                                                                       \n"
           "        -o output-filename-prefix                                                \n"
           "---------------------------------------------------------------------------------\n"
           " Options [description = default values]                                          \n"
           "---------------------------------------------------------------------------------\n"
           " Image parameters:                                                               \n"
           "    -sz   [Image sizes (cols rows slices) = %s %s %s]                            \n"
           " Filter bank parameters:                                                         \n"
           "    -ns   [Number of scales = %s]                                                \n"
           "    -na   [Number of azimuth angles = %s]                                        \n"
           "    -ne   [Number of elevation angles = %s]                                      \n"
           "    -maxf [Maximum central frequency = %s]                                       \n"
           "    -mulf [Multiplicative factor = %s]                                           \n"
           "    -fr   [Frequency bandwidth ratio = %s]                                       \n"
           "    -ar   [Angular spread ratio = %s]                                            \n"
           " Butterworth lowpass filter parameters:                                          \n"
           "    -lpo  [Lowpass filter order = %s]                                            \n"
           "    -lpc  [Lowpass filter cutoff = %s]                                           \n"
           " Other:                                                                          \n"
           "    -us   [Uniform sampling = true if given]                                     \n"
           "                                                                                 \n"
           " NOTE: for 2D images, use -sz cols rows 1 and -ne 1                              \n"
           "---------------------------------------------------------------------------------\n\n",
        argv[0],
        DEFAULT_IMAGE_SIZE_X,
        DEFAULT_IMAGE_SIZE_Y,
        DEFAULT_IMAGE_SIZE_Z,
        DEFAULT_NUM_SCALES,
        DEFAULT_NUM_AZIMUTHS,
        DEFAULT_NUM_ELEVATIONS,
        DEFAULT_MAX_FREQUENCY,
        DEFAULT_MULT_FACTOR,
        DEFAULT_FREQUENCY_RATIO,
        DEFAULT_ANGULAR_RATIO,
        DEFAULT_LOWPASS_ORDER,
        DEFAULT_LOWPASS_CUTOFF
    );
}


int main(int argc, char *argv[])
{
    CCmdLine cmdline;

    // Check for the help flag or for the lack of minimum required arguments.
    // Show the program usage if needed (and exit).
    if (cmdline.HasSwitch("-h") || cmdline.SplitLine(argc, argv) < 1) {
        show_usage(argv);
        return EXIT_FAILURE;
    }

    // --------------------------------------------------------------------------------------------
    // Step 1 : getting the input data.
    // --------------------------------------------------------------------------------------------

    // Try to get the minimum required arguments.
    // Show the program usage if needed (and exit).
    std::string filename_prefix;
    try
    {
        filename_prefix = cmdline.GetArgument("-o", 0);
    }
    catch (...) // catch any exception
    {
        show_usage(argv);
        return EXIT_FAILURE;
    }

    // Get the optional arguments.
    // Set default values for those not given by the user. 
    size_t size_x = static_cast<size_t>(
        atol(cmdline.GetSafeArgument("-sz", 0, DEFAULT_IMAGE_SIZE_X).c_str()));
    size_t size_y = static_cast<size_t>(
        atol(cmdline.GetSafeArgument("-sz", 1, DEFAULT_IMAGE_SIZE_Y).c_str()));
    size_t size_z = static_cast<size_t>(
        atol(cmdline.GetSafeArgument("-sz", 2, DEFAULT_IMAGE_SIZE_Z).c_str()));
    size_t num_scales = static_cast<size_t>(
        atol(cmdline.GetSafeArgument("-ns", 0, DEFAULT_NUM_SCALES).c_str()));
    size_t num_azimuths = static_cast<size_t>(
        atol(cmdline.GetSafeArgument("-na", 0, DEFAULT_NUM_AZIMUTHS).c_str()));
    size_t num_elevations = static_cast<size_t>(
        atol(cmdline.GetSafeArgument("-ne", 0, DEFAULT_NUM_ELEVATIONS).c_str()));
    float max_frequency = static_cast<float>(
        atof(cmdline.GetSafeArgument("-maxf", 0, DEFAULT_MAX_FREQUENCY).c_str()));
    float mult_factor = static_cast<float>(
        atof(cmdline.GetSafeArgument("-mulf", 0, DEFAULT_MULT_FACTOR).c_str()));
    float frequency_ratio = static_cast<float>(
        atof(cmdline.GetSafeArgument("-fr", 0, DEFAULT_FREQUENCY_RATIO).c_str()));
    float angular_ratio = static_cast<float>(
        atof(cmdline.GetSafeArgument("-ar", 0, DEFAULT_ANGULAR_RATIO).c_str()));
    float lowpass_order = static_cast<float>(
        atof(cmdline.GetSafeArgument("-lpo", 0, DEFAULT_LOWPASS_ORDER).c_str()));
    float lowpass_cutoff = static_cast<float>(
        atof(cmdline.GetSafeArgument("-lpc", 0, DEFAULT_LOWPASS_CUTOFF).c_str()));
    bool uniform_sampling = cmdline.HasSwitch("-us");

    bip::triple<size_t> sizes(size_x, size_y, size_z);

    // --------------------------------------------------------------------------------------------
    // Step 2: computing and writing the bank of log-gabor filters.
    // --------------------------------------------------------------------------------------------

    // Create the bank of filters.
    bip::loggabor_filter_bank lgab(filename_prefix, sizes, num_scales, num_azimuths,
                                   num_elevations, max_frequency, mult_factor, frequency_ratio,
                                   angular_ratio, lowpass_order, lowpass_cutoff, uniform_sampling);
    
    // Write the parameters of the bank of filters.
    bip::loggabor_filter_bank::write_parameters(lgab);
    
    // Compute and write all filters.
    lgab.compute();

    return EXIT_SUCCESS;
}
