/**
 * @file   __main.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup landmark-matching
 * @ingroup    landmark-matching
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
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include "CmdLine.h"
#include "assert.hpp"
#include "triple.hpp"
#include "landmark.hpp"
#include "imageio.hpp"
#include "landmarkio.hpp"
#include "matchlandmarks.hpp"

// Default argument values.
#define DEFAULT_DESCRIPTORS_TRADEOFF    "0.5"
#define DEFAULT_MAX_DESCRIPTOR_DISTANCE "0.0"
#define DEFAULT_MAX_LOCATION_DISTANCE   "0.0"
#define DEFAULT_MIN_LOCATION_DISTANCE   "0.0"

void show_usage(char *argv[])
{
    printf("---------------------------------------------------------------------------------\n"
           "                                Landmark Matching                                \n"
           "---------------------------------------------------------------------------------\n"
           " Usage: %s                                                                       \n"
           "        -io input-1-filename-prefix input-2-filename-prefix                      \n"
        //    "        -r  reference-image-1 reference-image-2                                  \n"
           "---------------------------------------------------------------------------------\n"
           " Options [description = default values]                                          \n"
           "---------------------------------------------------------------------------------\n"
           "    -to   [Trade-off between local and global descriptors = %s]                  \n"
           "    -maxd [Maximum distance of weighted landmark descriptors = %s]               \n"
           "    -maxl [Maximum distance of landmark locations = %s]                          \n"
           "    -minl [Minimum distance of landmark locations = %s]                          \n"
           " NOTE: for infinite maximum distances, use -maxd (-maxl) <= 0.0                  \n"
           "---------------------------------------------------------------------------------\n\n",
        argv[0],
        DEFAULT_DESCRIPTORS_TRADEOFF,
        DEFAULT_MAX_DESCRIPTOR_DISTANCE,
        DEFAULT_MAX_LOCATION_DISTANCE,
        DEFAULT_MIN_LOCATION_DISTANCE);
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
    std::string filename_input1_prefix;
    std::string filename_input2_prefix;
    // std::string filename_reference_img_1;
    // std::string filename_reference_img_2;
    try
    {
        filename_input1_prefix   = cmdline.GetArgument("-io", 0);
        filename_input2_prefix   = cmdline.GetArgument("-io", 1);
        // filename_reference_img_1 = cmdline.GetArgument("-r", 0);
        // filename_reference_img_2 = cmdline.GetArgument("-r", 1);
    }
    catch (...) // catch any exception
    {
        show_usage(argv);
        return EXIT_FAILURE;
    }

    // Get the optional arguments.
    // Set default values for those not given by the user.
    float descriptors_tradeoff = static_cast<float>(
        atof(cmdline.GetSafeArgument("-to", 0, DEFAULT_DESCRIPTORS_TRADEOFF).c_str()));
    float max_descriptor_dist = static_cast<float>(
        atof(cmdline.GetSafeArgument("-maxd", 0, DEFAULT_MAX_DESCRIPTOR_DISTANCE).c_str()));
    float max_location_dist = static_cast<float>(
        atof(cmdline.GetSafeArgument("-maxl", 0, DEFAULT_MAX_LOCATION_DISTANCE).c_str()));
    float min_location_dist = static_cast<float>(
        atof(cmdline.GetSafeArgument("-minl", 0, DEFAULT_MIN_LOCATION_DISTANCE).c_str()));

    // Get the filenames of the landmarks files.
    char filename_input1[512];
    char filename_input2[512];
    sprintf(filename_input1, "%s_landmarks.txt", filename_input1_prefix.c_str());
    sprintf(filename_input2, "%s_landmarks.txt", filename_input2_prefix.c_str());

    std::vector<bip::landmark> landmarks1;
    std::vector<bip::landmark> landmarks2;
    try
    {
        // Read the landmarks from both files.
        // std::cout << filename_input1 << filename_input2;
        landmarks1 = bip::read_landmarks(filename_input1);
        // std::cout << landmarks1 << filename_input2;
        landmarks2 = bip::read_landmarks(filename_input2);
        // std::cout << landmarks2;
    }
    catch (const char *exception)
    {
        std::cerr << "An exception has occurred!\n  " << exception << std::endl;
    }

    // --------------------------------------------------------------------------------------------
    // Step 2: finding matches between pairs of landmarks.
    // --------------------------------------------------------------------------------------------

    std::vector<std::pair<size_t, size_t> > matches = match_landmarks(landmarks1,
                                                                      landmarks2,
                                                                      descriptors_tradeoff,
                                                                      max_descriptor_dist,
                                                                      max_location_dist,
                                                                      min_location_dist);

    // --------------------------------------------------------------------------------------------
    // Step 3: write the results.
    // --------------------------------------------------------------------------------------------

    typedef itk::Image<float, 3> itkImage;

    // // Read the reference images.
    // itkImage::Pointer itk_reference_img_1 = bip::read_image<float, 3>(filename_reference_img_1);
    // itkImage::Pointer itk_reference_img_2 = bip::read_image<float, 3>(filename_reference_img_2);

    // // Create the matches images with the help of the reference images.
    // itkImage::Pointer itk_matches_img_1 = itkImage::New();
    // itk_matches_img_1->SetRegions(itk_reference_img_1->GetLargestPossibleRegion());
    // itk_matches_img_1->SetOrigin(itk_reference_img_1->GetOrigin());
    // itk_matches_img_1->SetSpacing(itk_reference_img_1->GetSpacing());
    // itk_matches_img_1->SetDirection(itk_reference_img_1->GetDirection());
    // itk_matches_img_1->Allocate();
    // itk_matches_img_1->FillBuffer(0);

    // itkImage::Pointer itk_matches_img_2 = itkImage::New();
    // itk_matches_img_2->SetRegions(itk_reference_img_2->GetLargestPossibleRegion());
    // itk_matches_img_2->SetOrigin(itk_reference_img_2->GetOrigin());
    // itk_matches_img_2->SetSpacing(itk_reference_img_2->GetSpacing());
    // itk_matches_img_2->SetDirection(itk_reference_img_2->GetDirection());
    // itk_matches_img_2->Allocate();
    // itk_matches_img_2->FillBuffer(0);

    char filename_output[1024];

    // Open an output file stream to write the landmark matches in a text file too.
    sprintf(filename_output, "%s_&_%s_matches.txt", filename_input1_prefix.c_str(),
                                                    filename_input2_prefix.c_str());
    std::ofstream ofs(filename_output);

    // float label = 1;
    // std::vector<std::pair<size_t, size_t> >::iterator it;

    // // For each pair of matched landmarks...
    // for (it = matches.begin(); it != matches.end(); ++it) {
    //     // Get the index of both matched landmarks.
    //     size_t m1 = it->first;
    //     size_t m2 = it->second;

    //     // Write the matched landmark locations in the text file.
    //     ofs << landmarks1[m1].get_location() << " <==> "
    //         << landmarks2[m2].get_location() << "\n";

    //     // Define the positions of the landmarks in the image.
    //     itkImage::IndexType index1;
    //     index1[0] = landmarks1[m1].get_location()[0];
    //     index1[1] = landmarks1[m1].get_location()[1];
    //     index1[2] = landmarks1[m1].get_location()[2];

    //     itkImage::IndexType index2;
    //     index2[0] = landmarks2[m2].get_location()[0];
    //     index2[1] = landmarks2[m2].get_location()[1];
    //     index2[2] = landmarks2[m2].get_location()[2];

    //     // Set in both images the voxel value as the label value that identifies the match.
    //     // itk_matches_img_1->SetPixel(index1, label);
    //     // itk_matches_img_2->SetPixel(index2, label);

    //     // The next match (if any) will receive another identification label.
    //     label += 1;
    // }
    // ofs.close();

    // Write the result as a pair of images whose matching landmarks are represented by the same
    // intensity values.
    // sprintf(filename_output, "%s_matches.nii", filename_input1_prefix.c_str());
    // bip::write_image<float, 3>(filename_output, itk_matches_img_1);

    // sprintf(filename_output, "%s_matches.nii", filename_input2_prefix.c_str());
    // bip::write_image<float, 3>(filename_output, itk_matches_img_2);

    return EXIT_SUCCESS;
}
