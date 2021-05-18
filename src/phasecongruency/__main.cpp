/**
 * @file   __main.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup phasecongruency
 * @ingroup    phasecongruency
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
#include <vector>
#include <itkImage.h>
#include "CmdLine.h"
#include "triple.hpp"
#include "loggaborfilterbank.hpp"
#include "phasecongruency.hpp"

// Default argument values.
#define DEFAULT_NOISE_THRESHOLD "-1.0"
#define DEFAULT_NOISE_STD       "3.0"
#define DEFAULT_SIGMOID_GAIN    "10.0"
#define DEFAULT_SIGMOID_CUTOFF  "0.5"


void show_usage(char *argv[])
{
    printf("---------------------------------------------------------------------------------\n"
           "                                Phase Congruency                                 \n"
           "---------------------------------------------------------------------------------\n"
           " Usage: %s                                                                       \n"
           "        -i   input-image                                                         \n"
           "        -bof bank-of-filters                                                     \n"
           "        -o   output-filename-prefix                                              \n"
           "---------------------------------------------------------------------------------\n"
           " Options [description = default values]                                          \n"
           "---------------------------------------------------------------------------------\n"
           " Phase congruency parameters:                                                    \n"
           "    -nt   [Noise energy threshold = %s]                                          \n"
           "    -nstd [Noise standard deviation = %s]                                        \n"
           "    -sig  [Sigmoid weight parameters (gain cutoff) = %s %s]                      \n"
           " Other:                                                                          \n"
           "    -mask [Binary mask of the region of interest = NULL]                         \n"
           "    -hann [Apply a Hanning window to the input image = true if given]            \n"
           "                                                                                 \n"
           " NOTE: for automatic noise threshold estimation, use nt -1.0 (or any value < 0)  \n"
           "---------------------------------------------------------------------------------\n\n",
        argv[0],
        DEFAULT_NOISE_THRESHOLD,
        DEFAULT_NOISE_STD,
        DEFAULT_SIGMOID_GAIN,
        DEFAULT_SIGMOID_CUTOFF
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
    // Step 1: getting the input data.
    // --------------------------------------------------------------------------------------------

    // Try to get the minimum required arguments.
    // Show the program usage if needed (and exit).
    std::string filename_input;
    std::string filename_bof;
    std::string filename_output_prefix;
    try
    {
        filename_input         = cmdline.GetArgument("-i", 0);
        filename_bof           = cmdline.GetArgument("-bof", 0);
        filename_output_prefix = cmdline.GetArgument("-o", 0);
    }
    catch (...) // catch any exception
    {
        show_usage(argv);
        return EXIT_FAILURE;
    }

    // Get the optional arguments.
    // Set default values for those not given by the user.
    float noise_threshold = static_cast<float>(
        atof(cmdline.GetSafeArgument("-nt", 0, DEFAULT_NOISE_THRESHOLD).c_str()));
    float noise_std = static_cast<float>(
        atof(cmdline.GetSafeArgument("-nstd", 0, DEFAULT_NOISE_STD).c_str()));
    float sigmoid_gain = static_cast<float>(
        atof(cmdline.GetSafeArgument("-sig", 0, DEFAULT_SIGMOID_GAIN).c_str()));
    float sigmoid_cutoff = static_cast<float>(
        atof(cmdline.GetSafeArgument("-sig", 1, DEFAULT_SIGMOID_CUTOFF).c_str()));
    std::string filename_mask = cmdline.GetSafeArgument("-mask", 0, "");
    bool hanning_flag = cmdline.HasSwitch("-hann");

    typedef itk::Image<float, 3> itkImage;

    // Read the input image and its data array.
    typename itkImage::Pointer itk_input_img = bip::read_image<float, 3>(filename_input);
    float *input_img = bip::image2array<float, 3>(itk_input_img);

    // Get the image sizes.
    bip::triple<size_t> sizes(
        itk_input_img->GetLargestPossibleRegion().GetSize()[0],
        itk_input_img->GetLargestPossibleRegion().GetSize()[1],
        itk_input_img->GetLargestPossibleRegion().GetSize()[2]
    );

    // Read the mask (if needed).
    unsigned char *input_mask = nullptr;
    if (!filename_mask.empty()) {
        input_mask = bip::image2array<unsigned char, 3>(
                     bip::read_image<unsigned char, 3>(filename_mask));
    }

    // --------------------------------------------------------------------------------------------
    // Step 2: computing the phase congruency.
    // --------------------------------------------------------------------------------------------

    // If required, apply a hanning window to the image data, in order to avoid high frequency
    // artifacts due to spectral leakage.
    if (hanning_flag) {
        for (size_t i = 0, z = 0; z < sizes[2]; ++z)
            for (size_t y = 0; y < sizes[1]; ++y)
                for (size_t x = 0; x < sizes[0]; ++x, ++i)
                    input_img[i] *= bip::hanning(bip::triple<size_t>(x, y, z), sizes);
    }

    // Read the bank of log-Gabor filters (i.e. its parameters).
    bip::loggabor_filter_bank *lgab = bip::loggabor_filter_bank::read_parameters(filename_bof);

    // Compute the phase congruency.
    bip::phase_congruency pc(filename_output_prefix, sizes, lgab, input_img, input_mask,
                             itk_input_img, noise_threshold, noise_std, sigmoid_gain,
                             sigmoid_cutoff);
    pc.compute();

    // Free memory.
    if (input_mask != nullptr)
        delete[] input_mask;
    delete[] input_img;
    delete lgab;


    return EXIT_SUCCESS;
}
