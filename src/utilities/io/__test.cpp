/**
 * @file   __test.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup io
 * @ingroup    io
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <itkImage.h>
#include <itkMesh.h>
#include "triple.hpp"
#include "imageio.hpp"
#include "meshio.hpp"


int main(int argc, char *argv[])
{
    const unsigned int Dimension = 3;

    typedef float                         TPixel;
    typedef float                         TPoint;
    typedef itk::Image<TPixel, Dimension> TImage;
    typedef itk::Mesh<TPixel, Dimension>  TMesh;

    // --------------------------------------------------------------------------------------------

    // Reading image from file.
    TImage::Pointer img = bip::read_image<TPixel, Dimension>("data/tipstruct.nii");

    // Getting array of pixels from image.
    TPixel *arr = bip::image2array<TPixel, Dimension>(img);

    bip::triple<size_t> sizes;
    sizes[0] = img->GetBufferedRegion().GetSize()[0];
    sizes[1] = img->GetBufferedRegion().GetSize()[1];
    sizes[2] = img->GetBufferedRegion().GetSize()[2];

    // Getting image from array of pixels.
    TImage::Pointer img2 = bip::array2image<TPixel, Dimension>(arr, sizes);

    // Writing image to file.
    bip::write_image<TPixel, Dimension>("tipstruct.nii", img2);

    delete[] arr;

    // --------------------------------------------------------------------------------------------

    // Reading mesh from file.
    TMesh::Pointer mesh = bip::read_mesh<TPoint, Dimension>("data/corpuscallosum.vtk");

    // Writing mesh to file.
    bip::write_mesh<TPoint, Dimension>("corpuscallosum.vtk", mesh);


    return EXIT_SUCCESS;
}
