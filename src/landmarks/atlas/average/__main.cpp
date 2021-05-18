/**
 * @file   __main.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup landmark-atlas-average
 * @ingroup    landmark-atlas-average
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

// #define BIP_DEBUG_MODE
// #define BIP_VERBOSE_MODE

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <deque>
#include <itkImage.h>
#include "CmdLine.h"
#include "assert.hpp"
#include "triple.hpp"
#include "landmark.hpp"
#include "atlassample.hpp"
#include "averagelandmarkatlas.hpp"

// Default argument values.
#define DEFAULT_RELATIVE_FREQUENCY_THRESHOLD "0.25"
#define DEFAULT_MAX_COUNTDOWN_VALUE          "1000"
#define DEFAULT_DESCRIPTORS_TRADEOFF         "0.5"
#define DEFAULT_MAX_DESCRIPTOR_DISTANCE      "0.0"
#define DEFAULT_MAX_LOCATION_DISTANCE        "0.0"


void show_usage(char *argv[])
{
    printf("---------------------------------------------------------------------------------\n"
           "                             Average Landmark Atlas                              \n"
           "---------------------------------------------------------------------------------\n"
           " Usage: %s                                                                       \n"
           "        -i samples-filenames-prefix                                              \n"
           "        -n num-samples                                                           \n"
           "        -r reference-image                                                       \n"
           "        -o output-filename-prefix                                                \n"
           "---------------------------------------------------------------------------------\n"
           " Options [description = default values]                                          \n"
           "---------------------------------------------------------------------------------\n"
           " Landmark averaging parameters:                                                  \n"
           "    -rft  [Relative frequency threshold = %s]                                    \n"
           "    -maxc [Maximum landmark countdown value = %s]                                \n"
           " Landmark matching parameters:                                                   \n"
           "    -to   [Trade-off between local and global descriptors = %s]                  \n"
           "    -maxd [Maximum distance of weighted landmark descriptors = %s]               \n"
           "    -maxl [Maximum distance of landmark locations = %s]                          \n"
           "                                                                                 \n"
           " NOTE: for infinite maximum distances, use -maxd (-maxl) <= 0.0                  \n"
           "---------------------------------------------------------------------------------\n\n",
        argv[0],
        DEFAULT_RELATIVE_FREQUENCY_THRESHOLD,
        DEFAULT_MAX_COUNTDOWN_VALUE,
        DEFAULT_DESCRIPTORS_TRADEOFF,
        DEFAULT_MAX_DESCRIPTOR_DISTANCE,
        DEFAULT_MAX_LOCATION_DISTANCE);
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
    std::string filename_samples_prefix;
    std::string filename_ref_img;
    std::string filename_output_prefix;
    size_t      num_samples;
    try
    {
        filename_samples_prefix = cmdline.GetArgument("-i", 0);
        filename_ref_img        = cmdline.GetArgument("-r", 0);
        filename_output_prefix  = cmdline.GetArgument("-o", 0);
        num_samples             = static_cast<size_t>(atol(cmdline.GetArgument("-n", 0).c_str()));
    }
    catch (...) // catch any exception
    {
        show_usage(argv);
        return EXIT_FAILURE;
    }

    // Get the optional arguments.
    // Set default values for those not given by the user.
    float rel_frequency_threshold = static_cast<float>(
        atof(cmdline.GetSafeArgument("-rft", 0, DEFAULT_RELATIVE_FREQUENCY_THRESHOLD).c_str()));
    size_t max_countdown_value = static_cast<size_t>(
        atol(cmdline.GetSafeArgument("-maxc", 0, DEFAULT_MAX_COUNTDOWN_VALUE).c_str()));
    float descriptors_tradeoff = static_cast<float>(
        atof(cmdline.GetSafeArgument("-to", 0, DEFAULT_DESCRIPTORS_TRADEOFF).c_str()));
    float max_descriptor_dist = static_cast<float>(
        atof(cmdline.GetSafeArgument("-maxd", 0, DEFAULT_MAX_DESCRIPTOR_DISTANCE).c_str()));
    float max_location_dist = static_cast<float>(
        atof(cmdline.GetSafeArgument("-maxl", 0, DEFAULT_MAX_LOCATION_DISTANCE).c_str()));

    // Read the reference image.
    typedef itk::Image<float, 3> itkImage;
    itkImage::Pointer itk_ref_img = bip::read_image<float, 3>(filename_ref_img);

    // Get the atlas sizes.
    bip::triple<size_t> sizes(
        itk_ref_img->GetBufferedRegion().GetSize()[0],
        itk_ref_img->GetBufferedRegion().GetSize()[1],
        itk_ref_img->GetBufferedRegion().GetSize()[2]);

    // --------------------------------------------------------------------------------------------
    // Step 2: computing the probabilistic atlas from samples.
    // --------------------------------------------------------------------------------------------

    bip::average_landmark_atlas atlas(filename_samples_prefix, sizes, num_samples, itk_ref_img,
                                      rel_frequency_threshold, max_countdown_value,
                                      descriptors_tradeoff, max_descriptor_dist,
                                      max_location_dist);
    try
    {
        atlas.compute();
    }
    catch (const char *exception)
    {
        std::cerr << "An exception has occurred!\n  " << exception << std::endl;
    }

    // --------------------------------------------------------------------------------------------
    // Step 3: writing the results.
    // --------------------------------------------------------------------------------------------

    char filename_output[512];

    const size_t total_size = sizes[0] * sizes[1] * sizes[2];

    const bip::atlas_sample *average_sample = atlas.get_average_sample();
    const std::vector<bip::landmark> landmarks(average_sample->get_landmarks().begin(),
                                               average_sample->get_landmarks().end());
    const std::deque<size_t> frequencies(average_sample->get_frequencies());

    // Write the averaged landmarks.
    sprintf(filename_output, "%s_averaged_landmarks.txt", filename_output_prefix.c_str());
    bip::write_landmarks(filename_output, landmarks);

    // Create a map of averaged landmarks. It is simply an image whose pixels indicate the
    // landmark locations and frequencies.
    float *average_landmarks_map = new float[total_size]();

    for (size_t i = 0; i < landmarks.size(); ++i) {
        const bip::triple<size_t> location = landmarks[i].get_location();
        const size_t j = location[0] + sizes[0] * (location[1] + sizes[1] * location[2]);

        average_landmarks_map[j] = frequencies[i];
    }

    // Write the landmarks map.
    bip::write_image<float, 3>(filename_output_prefix + "_averaged.nii",
         bip::array2image<float, 3>(average_landmarks_map, sizes, itk_ref_img));

    delete[] average_landmarks_map;


    return EXIT_SUCCESS;
}
