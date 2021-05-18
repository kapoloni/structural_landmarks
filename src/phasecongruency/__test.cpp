/**
 * @file   __test.cpp
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

#include <cstdlib>
#include <iostream>
#include <string>
#include "triple.hpp"
#include "imageio.hpp"
#include "loggaborfilterbank.hpp"
#include "phasecongruency.hpp"


int main(int argc, char *argv[])
{
    typedef float                 TPixel;
    typedef itk::Image<TPixel, 3> TImage;

    // --------------------------------------------------------------------------------------------

    bip::triple<size_t> sizes2(256, 256, 1);

    bip::loggabor_filter_bank bof2("loggabor2", sizes2, 4, 6, 1, 0.33, 2.1, 0.55, 1.2, 15, 0.45, false);
    bip::loggabor_filter_bank::write_parameters(bof2);
    bof2.compute();
    
    TImage::Pointer itkimg2 = bip::read_image<TPixel, 3>("data/sample2.tif");
    TPixel *img2 = bip::image2array<TPixel, 3>(itkimg2);

    bip::phase_congruency pc2("sample2", sizes2, &bof2, img2, NULL, itkimg2, -1.0, 2.0, 10, 0.5);
    pc2.compute();

    delete[] img2;

    // --------------------------------------------------------------------------------------------

    bip::triple<size_t> sizes3(33, 33, 33);

    bip::loggabor_filter_bank bof3("loggabor3", sizes3, 3, 6, 3, 0.33, 2.1, 0.55, 1.2, 15, 0.45, false);
    bip::loggabor_filter_bank::write_parameters(bof3);
    bof3.compute();
    
    TImage::Pointer itkimg3 = bip::read_image<TPixel, 3>("data/tipstruct.nii");
    TPixel *img3 = bip::image2array<TPixel, 3>(itkimg3);

    for (size_t i = 0, z = 0; z < sizes3[2]; ++z)
        for (size_t y = 0; y < sizes3[1]; ++y)
            for (size_t x = 0; x < sizes3[0]; ++x, ++i)
                img3[i] *= bip::hanning(bip::triple<size_t>(x, y, z), sizes3);

    bip::phase_congruency pc3("tipstruct", sizes3, &bof3, img3, NULL, itkimg3, -1.0, 2.0, 10, 0.5);
    pc3.compute();

    delete[] img3;


    return EXIT_SUCCESS;
}
