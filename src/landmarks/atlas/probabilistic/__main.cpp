/**
 * @file   __main.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup landmark-atlas-probabilistic
 * @ingroup    landmark-atlas-probabilistic
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
#include <itkImage.h>
#include "CmdLine.h"
#include "assert.hpp"
#include "triple.hpp"
#include "landmark.hpp"
#include "probabilisticlandmarkatlas.hpp"

// Default argument values.
#define DEFAULT_KERNEL_SIGMA              "1.0"
#define DEFAULT_REGULARIZATION_SIGMA      "1.0"
#define DEFAULT_EVALUATION_K_NEIGHBORHOOD "50"
#define DEFAULT_COVARIANCE_K_NEIGHBORHOOD "5"
#define DEFAULT_DENSITY_THRESHOLD         "0.25"
#define DEFAULT_AVERAGING_WINDOW_RADIUS   "5"


void show_usage(char *argv[])
{
    printf("---------------------------------------------------------------------------------\n"
           "                          Probabilistic Landmark Atlas                           \n"
           "---------------------------------------------------------------------------------\n"
           " Usage: %s                                                                       \n"
           "        -i samples-filenames-prefix                                              \n"
           "        -n num-samples                                                           \n"
           "        -r reference-image                                                       \n"
           "        -o output-filename-prefix                                                \n"
           "---------------------------------------------------------------------------------\n"
           " Options [description = default values]                                          \n"
           "---------------------------------------------------------------------------------\n"
           " Manifold Parzen Windows' parameters:                                            \n"
           "    -kstd [Standard deviation of Gaussian kernel = %s]                           \n"
           "    -rstd [Standard deviation of Gaussian regularization term = %s]              \n"
           "    -ek   [K-neighborhood value for evaluation = %s]                             \n"
           "    -ck   [K-neighborhood value for Gaussian covariances = %s]                   \n"
           " Landmarks averaging parameters:                                                 \n"
           "    -dt   [Relative density threshold = %s]                                      \n"
           "    -aw   [Averaging window size = %s]                                           \n"
           "---------------------------------------------------------------------------------\n\n",
        argv[0],
        DEFAULT_KERNEL_SIGMA,
        DEFAULT_REGULARIZATION_SIGMA,
        DEFAULT_EVALUATION_K_NEIGHBORHOOD,
        DEFAULT_COVARIANCE_K_NEIGHBORHOOD,
        DEFAULT_DENSITY_THRESHOLD,
        DEFAULT_AVERAGING_WINDOW_RADIUS);
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
    float kernel_sigma = static_cast<float>(
        atof(cmdline.GetSafeArgument("-kstd", 0, DEFAULT_KERNEL_SIGMA).c_str()));
    float regularization_sigma = static_cast<float>(
        atof(cmdline.GetSafeArgument("-rstd", 0, DEFAULT_REGULARIZATION_SIGMA).c_str()));
    size_t evaluation_k = static_cast<size_t>(
        atol(cmdline.GetSafeArgument("-ek", 0, DEFAULT_EVALUATION_K_NEIGHBORHOOD).c_str()));
    size_t covariance_k = static_cast<size_t>(
        atol(cmdline.GetSafeArgument("-ck", 0, DEFAULT_COVARIANCE_K_NEIGHBORHOOD).c_str()));
    float density_threshold = static_cast<float>(
        atof(cmdline.GetSafeArgument("-dt", 0, DEFAULT_DENSITY_THRESHOLD).c_str()));
    size_t averaging_window_radius = static_cast<size_t>(
        atol(cmdline.GetSafeArgument("-aw", 0, DEFAULT_AVERAGING_WINDOW_RADIUS).c_str()));

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

    bip::probabilistic_landmark_atlas atlas(filename_samples_prefix, sizes, num_samples,
                                            itk_ref_img, kernel_sigma, regularization_sigma,
                                            evaluation_k, covariance_k, density_threshold,
                                            averaging_window_radius);
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

    // Write the density image.
    sprintf(filename_output, "%s_density.nii", filename_output_prefix.c_str());
    bip::write_image<float, 3>(filename_output, atlas.get_density());

    const std::vector<bip::landmark> average_landmarks(atlas.get_average_landmarks());

    // Write the averaged landmarks.
    sprintf(filename_output, "%s_averaged_landmarks.txt", filename_output_prefix.c_str());
    bip::write_landmarks(filename_output, average_landmarks);

    // Create a map of averaged landmarks. It is simply an image with bright pixels indicating the
    // landmark locations.
    float *average_landmarks_map = new float[total_size]();

    for (size_t i = 0; i < average_landmarks.size(); ++i) {
        const bip::triple<size_t> location = average_landmarks[i].get_location();
        const size_t j = location[0] + sizes[0] * (location[1] + sizes[1] * location[2]);

        average_landmarks_map[j] = 1;
    }

    // And then we write such map.
    bip::write_image<float, 3>(filename_output_prefix + "_averaged.nii",
         bip::array2image<float, 3>(average_landmarks_map, sizes, itk_ref_img));

    delete[] average_landmarks_map;


    return EXIT_SUCCESS;
}
