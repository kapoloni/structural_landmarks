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
#define DEFAULT_GLOBAL_DECRIPTOR_MAX_RADIUS "64"
#define DEFAULT_LOCAL_DECRIPTOR_MAX_RADIUS  "16"
#define DEFAULT_NR "3"
#define DEFAULT_NA "8"



void show_usage(char *argv[])
{
    printf("---------------------------------------------------------------------------------\n"
           "                                Landmark Detector                                \n"
           "---------------------------------------------------------------------------------\n"
           " Usage: %s                                                                       \n"
           "        -io input-output-filename-prefix                                         \n"
           "        -r  reference-image                                                      \n"
           "        -m  map-name                                                             \n"
           "        -o  output-file                                                          \n"
           "---------------------------------------------------------------------------------\n"
           " Options [description = default values]                                          \n"
           "---------------------------------------------------------------------------------\n"
           " Landmark detector parameters:                                                   \n"
           "    -st   [Saliency threshold = %s]                                              \n"
           "    -max  [Maximum number of landmarks = %s]                                     \n"
           "    -gr   [Global descriptor maximum radius = %s]                                \n"
           "    -gl   [Local descriptor maximum radius = %s]                                 \n"
           "    -nr   [Number of radius radius = %s]                                         \n"
           "    -na   [Number of azimuths radius = %s]                                       \n"
           " Runtime options:                                                                \n"
           "    -geto [Get orientations = true if given]                                     \n"
           "    -getp [Get landmarks locations]                                              \n"
           "    -getg [Get global landmark descriptors = true if given]                      \n"
           "    -getl [Get local landmark descriptors = true if given]                       \n"
           " Other:                                                                          \n"
           "    -mask [Binary mask with region of interest = NULL]                           \n"
           "---------------------------------------------------------------------------------\n\n",
           argv[0],
           DEFAULT_SALIENCY_THRESHOLD,
           DEFAULT_NUM_MAX_LANDMARKS,
           DEFAULT_GLOBAL_DECRIPTOR_MAX_RADIUS,
           DEFAULT_LOCAL_DECRIPTOR_MAX_RADIUS,
           DEFAULT_NR,
           DEFAULT_NA);
}


int main(int argc, char *argv[])
{
    CCmdLine cmdline;

    // Check for the help flag or for the lack of minimum required arguments.
    // Show the program usage if needed (and exit).
    if (cmdline.HasSwitch("-h") || cmdline.SplitLine(argc, argv) < 3) {
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
        filename_ref_img = cmdline.GetArgument("-r", 0);
        map_name = cmdline.GetArgument("-m", 0);
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
    size_t global_descr_max_radius = static_cast<size_t>(
        atol(cmdline.GetSafeArgument("-gr", 0, DEFAULT_GLOBAL_DECRIPTOR_MAX_RADIUS).c_str()));
    size_t local_descr_max_radius = static_cast<size_t>(
        atol(cmdline.GetSafeArgument("-gl", 0, DEFAULT_LOCAL_DECRIPTOR_MAX_RADIUS).c_str()));
    std::string filename_mask = cmdline.GetSafeArgument("-mask", 0, "");

    int nr = static_cast<int>(atoi(cmdline.GetSafeArgument("-nr", 0, DEFAULT_NR).c_str()));
    int na = static_cast<int>(atoi(cmdline.GetSafeArgument("-na", 0, DEFAULT_NA).c_str()));

    bool flag_calc_land  = cmdline.HasSwitch("-getp");
    bool flag_calc_ori  = cmdline.HasSwitch("-geto");
    bool flag_global_desc = cmdline.HasSwitch("-getg");
    bool flag_local_desc = cmdline.HasSwitch("-getl");

    unsigned char options =
                            flag_calc_ori *
                                bip::landmark_detector::OPTION_CALCULATE_ORIENTATION    +
                            flag_calc_land *
                                bip::landmark_detector::OPTION_CALCULATE_LANDMARK +
                            flag_global_desc *
                                bip::landmark_detector::OPTION_GET_GLOBAL_DESCRIPTORS +
                            flag_local_desc *
                                bip::landmark_detector::OPTION_GET_LOCAL_DESCRIPTORS;

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
                                    saliency_threshold, num_max_landmarks, local_descr_max_radius,
                                    global_descr_max_radius, nr, na, options);

    char filename[512];

    std::string filename_detected_points = "None";
    sprintf(filename, "%s", filename_detected_points.c_str());

    detector.compute(filename);

    std::vector<bip::landmark> landmarks = detector.get_landmarks();
    std::cout << "Number of detected landmarks: " << landmarks.size() << std::endl;

    // Free memory.
    if (input_mask != NULL)
        delete[] input_mask;

    // --------------------------------------------------------------------------------------------
    // Step 3 writing the results.
    // --------------------------------------------------------------------------------------------
    string file_output;
    try
    {
        file_output  = cmdline.GetArgument("-o", 0);
    }
    catch (...) // catch any exception
    {
        show_usage(argv);
        return EXIT_FAILURE;
    }
    std::stringstream filename_output;
    filename_output << file_output << "_landmarks.txt";

    bip::write_landmarks( filename_output.str(), landmarks);

    const size_t total_size = sizes[0] * sizes[1] * sizes[2];
    // Create a map of landmarks. It is simply an image with bright pixels indicating the
    // landmark locations.
    float *landmarks_map = new float[total_size]();

    for (size_t i = 0; i < landmarks.size(); ++i) {
        const bip::triple<size_t> location = landmarks[i].get_location();
        const size_t j = location[0] + sizes[0] * (location[1] + sizes[1] * location[2]);

        landmarks_map[j] = 1;
    }

    // std::stringstream filename_output2;
    // filename_output2 << file_output << "_landmarks.nii";

    // // And then we write such map.
    // bip::write_image<float, 3>(filename_output2.str(),
    //      bip::array2image<float, 3>(landmarks_map, sizes, itk_ref_img));

    delete[] landmarks_map;


    return EXIT_SUCCESS;
}
