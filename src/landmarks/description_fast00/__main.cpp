/**
 * @file   __main.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup landmark-detection
 * @ingroup    landmark-detection
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

// #define BIP_DEBUG_MODE
// #define BIP_VERBOSE_MODE

#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <itkImage.h>
#include "CmdLine.h"
#include "assert.hpp"
#include "triple.hpp"
#include "landmark.hpp"
#include "mathfunctions.hpp"
#include "imageio.hpp"
#include "landmarkdetector.hpp"

// Default argument values.
#define DEFAULT_SALIENCY_THRESHOLD          "0.5"
#define DEFAULT_NUM_MAX_LANDMARKS           "10000"
#define DEFAULT_LOCAL_DECRIPTOR_REGION_SIZE "1" //"16"
#define DEFAULT_GLOBAL_DECRIPTOR_MAX_RADIUS "128"


void show_usage(char *argv[])
{
    printf("---------------------------------------------------------------------------------\n"
           "                                Landmark Detector                                \n"
           "---------------------------------------------------------------------------------\n"
           " Usage: %s                                                                       \n"
           "        -io input-output-filename-prefix                                         \n"
           "        -r  reference-image                                                      \n"
           "        -m  map-name                                                             \n"
           "---------------------------------------------------------------------------------\n"
           " Options [description = default values]                                          \n"
           "---------------------------------------------------------------------------------\n"
           " Landmark detector parameters:                                                   \n"
           "    -st   [Saliency threshold = %s]                                              \n"
           "    -max  [Maximum number of landmarks = %s]                                     \n"
           "    -gr   [Global descriptor maximum radius = %s]                                \n"
           " Runtime options:                                                                \n"
           "    -getg [Get global landmark descriptors = true if given]                      \n"
           " Other:                                                                          \n"
           "    -mask [Binary mask with region of interest = NULL]                           \n"
           "---------------------------------------------------------------------------------\n\n",
        argv[0],
        DEFAULT_SALIENCY_THRESHOLD,
        DEFAULT_NUM_MAX_LANDMARKS,
        DEFAULT_LOCAL_DECRIPTOR_REGION_SIZE,
        DEFAULT_GLOBAL_DECRIPTOR_MAX_RADIUS
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
    std::string filename_io_prefix;
    std::string filename_ref_img;
    std::string map_name;
    try
    {
        filename_io_prefix = cmdline.GetArgument("-io", 0);
        filename_ref_img   = cmdline.GetArgument("-r", 0);
        map_name           = cmdline.GetArgument("-m", 0);
    }
    catch (...) // catch any exception
    {
        show_usage(argv);
        return EXIT_FAILURE;
    }

    // Get the optional arguments.
    // Set default values for those not given by the user.
    float saliency_threshold = static_cast<float>(
        atof(cmdline.GetSafeArgument("-st", 0, DEFAULT_SALIENCY_THRESHOLD).c_str()));
    size_t num_max_landmarks = static_cast<size_t>(
        atol(cmdline.GetSafeArgument("-max", 0, DEFAULT_NUM_MAX_LANDMARKS).c_str()));
    size_t local_descr_region_size = static_cast<size_t>(
        atol(cmdline.GetSafeArgument("-lr", 0, DEFAULT_LOCAL_DECRIPTOR_REGION_SIZE).c_str()));
    size_t global_descr_max_radius = static_cast<size_t>(
        atol(cmdline.GetSafeArgument("-gr", 0, DEFAULT_GLOBAL_DECRIPTOR_MAX_RADIUS).c_str()));
    std::string filename_mask = cmdline.GetSafeArgument("-mask", 0, "");
    
    bool flag_saliencies  = cmdline.HasSwitch("-gets");
    bool flag_local_desc  = cmdline.HasSwitch("-getl");
    bool flag_global_desc = cmdline.HasSwitch("-getg");

    unsigned char options = flag_saliencies *
                                bip::landmark_detector::OPTION_GET_SALIENCIES_MAP    +
                            flag_local_desc *
                                bip::landmark_detector::OPTION_GET_LOCAL_DESCRIPTORS +
                            flag_global_desc *
                                bip::landmark_detector::OPTION_GET_GLOBAL_DESCRIPTORS;

    // Read the reference image.
    typedef itk::Image<float, 3>  itkImage;
    typename itkImage::Pointer itk_ref_img = bip::read_image<float, 3>(filename_ref_img);

    // Get the image sizes.
    bip::triple<size_t> sizes(
        itk_ref_img->GetLargestPossibleRegion().GetSize()[0],
        itk_ref_img->GetLargestPossibleRegion().GetSize()[1],
        itk_ref_img->GetLargestPossibleRegion().GetSize()[2]
    );

    // Read the mask (if needed).
    unsigned char *input_mask = NULL;
    if (!filename_mask.empty()) {
        input_mask = bip::image2array<unsigned char, 3>(
                     bip::read_image<unsigned char, 3>(filename_mask));
    }

    // --------------------------------------------------------------------------------------------
    // Step 2: detecting the landmarks and computing their descriptors.
    // --------------------------------------------------------------------------------------------

    bip::landmark_detector detector(filename_io_prefix, map_name, sizes, itk_ref_img, input_mask,
                                    saliency_threshold, num_max_landmarks, local_descr_region_size,
                                    global_descr_max_radius, options);
    detector.compute();

    std::vector<bip::landmark> landmarks = detector.get_landmarks();
    std::cout << "Number of detected landmarks: " << landmarks.size() << std::endl;

    // Free memory.
    if (input_mask != NULL)
        delete[] input_mask;

    // --------------------------------------------------------------------------------------------
    // Step 3 writing the results.
    // --------------------------------------------------------------------------------------------

    bip::write_landmarks(filename_io_prefix + "_landmarks.txt", landmarks);

    return EXIT_SUCCESS;
}
