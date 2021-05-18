/**
 * @file   __main.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup landmark-registration
 * @ingroup    landmark-registration
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
#include <utility>
#include <itkImage.h>
#include <itkVector.h>
#include <itkMesh.h>
#include "CmdLine.h"
#include "assert.hpp"
#include "landmark.hpp"
#include "imageio.hpp"
#include "landmarkio.hpp"
#include "landmarkguidedregistration.hpp"

// Default argument values.
#define DEFAULT_DESCRIPTORS_TRADEOFF    "0.5"
#define DEFAULT_MAX_DESCRIPTOR_DISTANCE "0.0"
#define DEFAULT_MAX_LOCATION_DISTANCE   "0.0"
#define DEFAULT_SPLINE_ORDER            "3"
#define DEFAULT_NUM_CONTROL_POINTS      "5"


void show_usage(char *argv[])
{
    printf("---------------------------------------------------------------------------------\n"
           "                      Landmark-Guided Nonrigid Registration                      \n"
           "---------------------------------------------------------------------------------\n"
           " Usage: %s                                                                       \n"
           "        -io fixed-image-filename-prefix moving-image-filename-prefix             \n"
           "---------------------------------------------------------------------------------\n"
           " Options [description = default values]                                          \n"
           "---------------------------------------------------------------------------------\n"
           " Registration parameters:                                                        \n"
           "    -so   [Spĺine order = %s (CURRENTLY UNCHANGEABLE)]                           \n"
           "    -cp   [Number of control points per dimension at the first level = %s]       \n"
           " Landmark matching parameters:                                                   \n"
           "    -to   [Trade-off between local and global descriptors = %s]                  \n"
           "    -maxd [Maximum distance of weighted landmark descriptors = %s]               \n"
           "    -maxl [Maximum distance of landmark locations = %s]                          \n"
           " Other:                                                                          \n"
           "    -wm   [Landmark weights map = NULL]                                          \n"
           "    -mesh [Meshes to be registered along with the moving image = NULL]           \n"
           "                                                                                 \n"
           " NOTE: for infinite maximum distances, use -maxd (-maxl) <= 0.0                  \n"
           "---------------------------------------------------------------------------------\n\n",
        argv[0],
        DEFAULT_SPLINE_ORDER,
        DEFAULT_NUM_CONTROL_POINTS,
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
    std::string filename_fixed_img_prefix;
    std::string filename_moving_img_prefix;
    try
    {
        filename_fixed_img_prefix  = cmdline.GetArgument("-io", 0);
        filename_moving_img_prefix = cmdline.GetArgument("-io", 1);
    }
    catch (...) // catch any exception
    {
        show_usage(argv);
        return EXIT_FAILURE;
    }

    // Get the optional arguments.
    // Set default values for those not given by the user.
    size_t spline_order = static_cast<size_t>(
        atol(cmdline.GetSafeArgument("-so", 0, DEFAULT_SPLINE_ORDER).c_str()));
    size_t num_control_points = static_cast<size_t>(
        atol(cmdline.GetSafeArgument("-cp", 0, DEFAULT_NUM_CONTROL_POINTS).c_str()));
    float descriptors_tradeoff = static_cast<float>(
        atof(cmdline.GetSafeArgument("-to", 0, DEFAULT_DESCRIPTORS_TRADEOFF).c_str()));
    float max_descriptor_dist = static_cast<float>(
        atof(cmdline.GetSafeArgument("-maxd", 0, DEFAULT_MAX_DESCRIPTOR_DISTANCE).c_str()));
    float max_location_dist = static_cast<float>(
        atof(cmdline.GetSafeArgument("-maxl", 0, DEFAULT_MAX_LOCATION_DISTANCE).c_str()));
    std::string filename_landmark_weights_map = cmdline.GetSafeArgument("-wm", 0, "");

    size_t num_meshes = std::max(0, cmdline.GetArgumentCount("-mesh"));
    std::vector<std::string> filenames_moving_meshes;
    for (size_t i = 0; i < num_meshes; ++i)
        filenames_moving_meshes.push_back(cmdline.GetSafeArgument("-mesh", i, ""));

    // Get the filenames of the landmarks files.
    char filename_fixed_img_landmarks[512];
    char filename_moving_img_landmarks[512];
    sprintf(filename_fixed_img_landmarks,  "%s_landmarks.txt", filename_fixed_img_prefix.c_str());
    sprintf(filename_moving_img_landmarks, "%s_landmarks.txt", filename_moving_img_prefix.c_str());

    std::vector<bip::landmark> fixed_img_landmarks;
    std::vector<bip::landmark> moving_img_landmarks;
    try
    {
        // Read the landmarks from both files.
        fixed_img_landmarks  = bip::read_landmarks(filename_fixed_img_landmarks);
        moving_img_landmarks = bip::read_landmarks(filename_moving_img_landmarks);
    }
    catch (const char *exception)
    {
        std::cerr << "An exception has occurred!\n  " << exception << std::endl;
    }

    typedef itk::Image<float, 3>          itkImage;
    typedef itk::Vector<float, 3>         itkVector;
    typedef itk::Image<itkVector, 3>      itkDisplField;
    typedef itk::Mesh<float, 3>           itkMesh;
    typedef std::vector<itkMesh::Pointer> itkMeshContainer;

    // Read the fixed and moving images.
    itkImage::Pointer itk_fixed_img = bip::read_image<float, 3>(
        filename_fixed_img_prefix + ".nii");
    itkImage::Pointer itk_moving_img = bip::read_image<float, 3>(
        filename_moving_img_prefix + ".nii");

    // Read the moving meshes (if needed).
    itkMeshContainer itk_moving_meshes;
    for (size_t i = 0; i < num_meshes; ++i)
        itk_moving_meshes.push_back(bip::read_mesh<float, 3>(filenames_moving_meshes[i]));

    // Read the weights (densities) map (if needed).
    itkImage::Pointer itk_landmark_weights_map = nullptr;

    if (!filename_landmark_weights_map.empty())
        itk_landmark_weights_map = bip::read_image<float, 3>(filename_landmark_weights_map);

    // --------------------------------------------------------------------------------------------
    // Step 2: performing the registration guided by pairs of matched landmarks.
    // --------------------------------------------------------------------------------------------

    bip::landmark_guided_registration reg(itk_fixed_img, itk_moving_img, itk_moving_meshes,
                                          fixed_img_landmarks, moving_img_landmarks,
                                          num_control_points, descriptors_tradeoff,
                                          max_descriptor_dist, max_location_dist,
                                          itk_landmark_weights_map);
    reg.compute();

    // --------------------------------------------------------------------------------------------
    // Step 3: writing the results.
    // --------------------------------------------------------------------------------------------

    char filename_output[512];

    // Get the results.
    itkImage::Pointer itk_output_img               = reg.get_output_img();
    itkMeshContainer  itk_output_meshes            = reg.get_output_meshes();
    itkDisplField::Pointer itk_displ_field_img     = reg.get_displacement_field();
    itkDisplField::Pointer itk_inv_displ_field_img = reg.get_inv_displacement_field();

    // Write the registered image.
    sprintf(filename_output, "%s_registration.nii", filename_moving_img_prefix.c_str());
    bip::write_image<float, 3>(filename_output, itk_output_img);

    // Write the registered meshes.
    for (int i = 0; i < num_meshes; ++i) {
        sprintf(filename_output, "%s_registration.vtk",
                filenames_moving_meshes[i].substr(0, filenames_moving_meshes[i].size()-4).c_str());
        bip::write_mesh<float, 3>(filename_output, itk_output_meshes[i]);
    }

    // Write the displacement field (applied to the moving meshes).
    sprintf(filename_output, "%s_displfield.nii", filename_moving_img_prefix.c_str());
    bip::write_image<itkVector, 3>(filename_output, itk_displ_field_img);

    // Write the inverse displacement field (applied to the moving image).
    sprintf(filename_output, "%s_inv_displfield.nii", filename_moving_img_prefix.c_str());
    bip::write_image<itkVector, 3>(filename_output, itk_inv_displ_field_img);


    return EXIT_SUCCESS;
}
