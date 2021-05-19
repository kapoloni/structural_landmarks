/**
 * @file   __test.cpp
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

#include <cstdlib>
#include <iostream>
#include <string>
#include "triple.hpp"
#include "imageio.hpp"
#include "loggaborfilterbank.hpp"
#include "phasecongruency.hpp"
#include "landmark.hpp"
#include "landmarkio.hpp"
#include "landmarkdetector.hpp"


int main(int argc, char *argv[])
{
    typedef float                 TPixel;
    typedef unsigned char         TLabel;
    typedef itk::Image<TPixel, 3> TImage;

    // --------------------------------------------------------------------------------------------

    bip::triple<size_t> sizes2(256, 256, 1);

    TImage::Pointer itkimg2 = bip::read_image<TPixel, 3>("data/sample2.tif");

    bip::landmark_detector detector2("sample2", sizes2, itkimg2, NULL, 0.6, 1E4, 16, 256, 6);
    detector2.compute();

    std::vector<bip::landmark> landmarks2 = detector2.get_landmarks();
    std::cout << "Number of detected landmarks: " << landmarks2.size() << std::endl;

    bip::write_landmarks("sample2_landmarks.txt", landmarks2);

    TPixel *arr2 = new TPixel[sizes2[0] * sizes2[1] * sizes2[2]]();

    for (size_t i = 0; i < landmarks2.size(); ++i) {
        const bip::triple<size_t> location = landmarks2[i].get_location();
        const size_t j = location[0] + sizes2[0] * (location[1] + sizes2[1] * location[2]);

        arr2[j] = 1;
    }

    bip::write_image<TPixel, 3>("sample2_landmarks.nii",
         bip::array2image<TPixel, 3>(arr2, sizes2, itkimg2));

    delete[] arr2;

    // --------------------------------------------------------------------------------------------

    bip::triple<size_t> sizes3(33, 33, 33);

    TImage::Pointer itkimg3 = bip::read_image<TPixel, 3>("data/tipstruct.nii");

    bip::landmark_detector detector3("tipstruct", sizes3, itkimg3, NULL, 0.35, 1E4, 16, 33, 6);
    detector3.compute();

    std::vector<bip::landmark> landmarks3 = detector3.get_landmarks();
    std::cout << "Number of detected landmarks: " << landmarks3.size() << std::endl;

    bip::write_landmarks("tipstruct_landmarks.txt", landmarks3);

    TPixel *arr3 = new TPixel[sizes3[0] * sizes3[1] * sizes3[2]]();

    for (size_t i = 0; i < landmarks3.size(); ++i) {
        const bip::triple<size_t> location = landmarks3[i].get_location();
        const size_t j = location[0] + sizes3[0] * (location[1] + sizes3[1] * location[2]);

        arr3[j] = 1;
    }

    bip::write_image<TPixel, 3>("tipstruct_landmarks.nii",
         bip::array2image<TPixel, 3>(arr3, sizes3, itkimg3));

    delete[] arr3;


    return EXIT_SUCCESS;
}
