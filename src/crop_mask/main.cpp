#include <itkImage.h>
#include <itkRescaleIntensityImageFilter.h>

#include <itkScalarImageToTextureFeaturesFilter.h>
#include <itkScalarImageToCooccurrenceMatrixFilter.h>
#include <itkRegionOfInterestImageFilter.h>
#include <itkImageMomentsCalculator.h>
#include <itkExtractImageFilter.h>
#include <itkChangeInformationImageFilter.h>

#include "CmdLine.h"
#include "bipUtils.h"
#include "bipMiscellaneous.h"

#define DEFAULT_SIZE "32"

void show_usage(char *argv[])
{
    printf("---------------------------------------------------------------------------------\n"
           "                                Crop Image Mask                                  \n"
           "---------------------------------------------------------------------------------\n"
           " Usage: %s                                                                       \n"
           "        -i input image                                                           \n"
           "        -m input mask                                                            \n"
           "        -o output name                                                           \n"
           "        -s [size/2, default =  %s -> final-size = 2*%s]                          \n"
           "---------------------------------------------------------------------------------\n\n",
        argv[0],
        DEFAULT_SIZE,
        DEFAULT_SIZE
    );
}

int main( int argc, char* argv[] )
{
    CCmdLine cmdLine;

    if (cmdLine.HasSwitch("-h") || cmdLine.SplitLine(argc, argv) < 3) {
        show_usage(argv);
        return EXIT_FAILURE;
    }

    // Setup types
    using TImageType = itk::Image< double, 3 >;
    using TMaskType = itk::Image< unsigned char, 3 >;

    // Read input image
    TImageType::Pointer input_img = bip::utils::ReadImage< TImageType >(cmdLine.GetArgument("-i", 0));

    // Read input mask
    TMaskType::Pointer input_mask = bip::utils::ReadImage< TMaskType >(cmdLine.GetArgument("-m", 0));

    using FixedImageCalculatorType = itk::ImageMomentsCalculator< TMaskType >;
    FixedImageCalculatorType::Pointer fixedCalculator = FixedImageCalculatorType::New();
    fixedCalculator->SetImage( input_mask );
    fixedCalculator->Compute();
    FixedImageCalculatorType::VectorType fixedCenter = fixedCalculator->GetCenterOfGravity();

    TImageType::IndexType pixelIndex;
    TImageType::PointType point = fixedCenter;
    input_mask->TransformPhysicalPointToIndex(point, pixelIndex);

    TImageType::IndexType desiredStart;
    int dim = atoi(cmdLine.GetSafeArgument("-s", 0, DEFAULT_SIZE).c_str());
    desiredStart[0] = pixelIndex[0] - dim;
    desiredStart[1] = pixelIndex[1] - dim;
    desiredStart[2] = pixelIndex[2] - dim;

    TImageType::SizeType desiredSize;
    desiredSize.Fill(dim*2);

    TImageType::RegionType desiredRegion(desiredStart, desiredSize);

    using FilterType = itk::ExtractImageFilter< TImageType, TImageType >;
    FilterType::Pointer filter = FilterType::New();
    filter->SetExtractionRegion( desiredRegion );
    filter->SetInput( input_img );
    filter->SetDirectionCollapseToIdentity();
    filter->Update();

    TImageType::Pointer output = filter->GetOutput();

    // Rescale image
    //std::string *intensity = bip::misc::get_min_max_intensity< TImageType >(output);
    //float value = atof(intensity[1].c_str()) - atof(intensity[0].c_str());
    //output = bip::misc::rescale_image< TImageType, TImageType >(output, 0, value);

    bip::utils::WriteImage< TImageType >(output, cmdLine.GetArgument("-o", 0));

  return 0;
}
