#include <itkGrayscaleDilateImageFilter.h>
#include <itkFlatStructuringElement.h>

#include "CmdLine.h"
#include "bipUtils.h"

void show_usage(char *argv[])
{
    printf("---------------------------------------------------------------------------------\n"
           "                              Dilate mask                                        \n"
           "---------------------------------------------------------------------------------\n"
           "              Usage: %s -i image-filename  -r ratio   -o output-name             \n"
           "---------------------------------------------------------------------------------\n\n",
        argv[0]
    );
}
int main(int argc, char *argv[])
{
  CCmdLine cmdLine;
  if( cmdLine.HasSwitch("-h") || cmdLine.SplitLine(argc, argv) < 3 )
  {
      show_usage(argv);
      return EXIT_FAILURE;
  }

  // Setup types
  using ImageType = itk::Image< double, 3 >;

  // Read input image
  ImageType::Pointer input_img = bip::utils::ReadImage< ImageType >(cmdLine.GetArgument("-i", 0));

  using StructuringElementType = itk::FlatStructuringElement< 3 >;
  StructuringElementType::RadiusType radius;
  radius.Fill( atoi( cmdLine.GetArgument("-r", 0).c_str() ) );
  StructuringElementType structuringElement = StructuringElementType::Ball( radius );

  using GrayscaleDilateImageFilterType = itk::GrayscaleDilateImageFilter< ImageType, ImageType, StructuringElementType >;
  GrayscaleDilateImageFilterType::Pointer dilateFilter = GrayscaleDilateImageFilterType::New();
  dilateFilter->SetInput( input_img );
  dilateFilter->SetKernel( structuringElement );
  dilateFilter->Update();

  // Write
  bip::utils::WriteImage< ImageType >(dilateFilter->GetOutput(), cmdLine.GetArgument("-o", 0));

  return EXIT_SUCCESS;
}
