#include "CmdLine.h"
#include "bipUtils.h"
#include "itkExtractImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkFlipImageFilter.h"

void show_usage(char *argv[])
{
    printf("---------------------------------------------------------------------------------\n"
           "                               Slice Image                                       \n"
           "---------------------------------------------------------------------------------\n"
           " Usage: %s        -i input-image -o output-name                                  \n"
           "---------------------------------------------------------------------------------\n\n",
        argv[0]
    );
}


int main(int argc, char *argv[])
{
    CCmdLine cmdLine;
    if( cmdLine.HasSwitch("-h") || cmdLine.SplitLine(argc, argv) < 1 )
    {
        show_usage(argv);
        return EXIT_FAILURE;
    }

    // Setup types
    using ImageType = itk::Image<float, 3>;
    using ImageType2D = itk::Image<unsigned char, 2>;

    // Read input image
    ImageType::Pointer input_img = bip::utils::ReadImage< ImageType >(cmdLine.GetArgument("-i", 0));
    std::string output_name = cmdLine.GetArgument("-o", 0);

    unsigned int which_dimension = 1;
    unsigned int slice_number = input_img->GetLargestPossibleRegion().GetSize()[which_dimension] / 2;

    // Get the size and set which_dimension to zero.
    ImageType::RegionType inputRegion = input_img->GetLargestPossibleRegion();
    ImageType::SizeType size = inputRegion.GetSize();
    size[which_dimension] = 0;

    // Get the index and set which_dimension to the correct slice.
    ImageType::IndexType start = inputRegion.GetIndex();
    start[which_dimension] = slice_number;

    // Create a desired extraction region and set it into the extractor.
    ImageType::RegionType desiredRegion;
    desiredRegion.SetSize(size);
    desiredRegion.SetIndex(start);

    ImageType2D::Pointer slicePtr;
    {
      itk::FixedArray<bool, 2> flipAxes;
      flipAxes[0] = 0;
      flipAxes[1] = 1;

      typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescaleIntensityImageFilterType;
      typename RescaleIntensityImageFilterType::Pointer rescaler = RescaleIntensityImageFilterType::New();
      rescaler->SetInput(input_img);
      rescaler->SetOutputMinimum(0);
      rescaler->SetOutputMaximum(255);

      typedef itk::ExtractImageFilter<ImageType, ImageType2D> ExtractImageFilterType;
      typename ExtractImageFilterType::Pointer extractor = ExtractImageFilterType::New();
      extractor->SetInput(rescaler->GetOutput());
      extractor->SetDirectionCollapseToSubmatrix();
      extractor->SetExtractionRegion(desiredRegion);

      typedef itk::FlipImageFilter<ImageType2D> FlipImageFilterType;
      typename FlipImageFilterType::Pointer fliper = FlipImageFilterType::New();
      fliper->SetInput(extractor->GetOutput());
      fliper->SetFlipAxes(flipAxes);
      fliper->Update();

      slicePtr = fliper->GetOutput();
    }

    using WriterType = itk::ImageFileWriter<ImageType2D>;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(output_name);
    writer->SetInput(slicePtr);
    writer->Update();


    return EXIT_SUCCESS;
}
