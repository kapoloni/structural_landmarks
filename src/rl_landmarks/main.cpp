#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

#include <itkImage.h>
#include <itkRescaleIntensityImageFilter.h>

#include <itkRegionOfInterestImageFilter.h>
#include <itkImageMomentsCalculator.h>

#include <itkHistogramToRunLengthFeaturesFilter.h>
#include <itkScalarImageToRunLengthMatrixFilter.h>

#include "CmdLine.h"
#include "bipMiscellaneous.h"
#include "bipUtils.h"

#define DEFAULT_BINS "64"
#define DEFAULT_PIXEL_MIN "0"
#define DEFAULT_PIXEL_MAX "255"
#define DEFAULT_RADIUS "2"
#define DEFAULT_PATCH "5"
// #define DEFAULT_SIZE "64"


void showUsage(char *argv[])
{
    printf("---------------------------------------------------------------------------------\n"
           "                    Texture Analysis RL (position [x,y,z])                       \n"
           "---------------------------------------------------------------------------------\n"
           " Usage: %s                                                                       \n"
           "        -i input image                                                           \n"
           "        -p file points                                                           \n"
           "---------------------------------------------------------------------------------\n"
           " Options [description = default values]                                          \n"
           "---------------------------------------------------------------------------------\n"
           " RUN LENGHT parameters:                                                          \n"
          //  "   -s  [img roi size, default =  %s]                                             \n"
           "   -nb [number of bins per axis, default = %s]                                   \n"
           "   -pv [pixel value (min, max), default = (%s, %s)]                              \n"
           "   -nr [neighborhood radius, default =  %s]                                      \n"
           "   -ps [patch size, default =  %s]                                               \n"
           "---------------------------------------------------------------------------------\n\n",
        argv[0],
        DEFAULT_BINS,
        DEFAULT_PIXEL_MIN,
        DEFAULT_PIXEL_MAX,
        DEFAULT_RADIUS,
        DEFAULT_PATCH
        // DEFAULT_SIZE
    );
}

// Setup types
using InputImageType = itk::Image< double, 3 >;
using OutputImageType = itk::Image< itk::Vector< double, 8 > , 3 >;
using ReaderType = itk::ImageFileReader< InputImageType >;
using OffsetVector = itk::VectorContainer< unsigned char, InputImageType::OffsetType >;

void runlenght_texture_features(
                              OffsetVector::Pointer offset,
                              InputImageType::Pointer inputImage,
                              int number_bins, int patch_size,
                              float pixel_min, float pixel_max,
                              vector<InputImageType::IndexType> positions
                            );

InputImageType::Pointer rescale_image(
                        InputImageType::Pointer inputImage,
                        float pixel_min, float pixel_max
                        );

vector<InputImageType::IndexType> get_point_positions(
                                                    std::string landmarks_parameters_filename
                                                  );

int main(int argc, char **argv)
{
  CCmdLine cmdLine;

  if(cmdLine.SplitLine(argc, argv) < 3 || (!cmdLine.HasSwitch("-i") || !cmdLine.HasSwitch("-p")))
  {
    showUsage(argv);
    return -1;
  }

  // Read input image
  InputImageType::Pointer input_img = bip::utils::ReadImage< InputImageType >(cmdLine.GetArgument("-i", 0));// Read input mask

  // Read input points
  std::string input_points = cmdLine.GetArgument("-p", 0);

  // Set default values for those not given by the user.
  int radius = std::atoi(cmdLine.GetSafeArgument("-nr", 0, DEFAULT_RADIUS).c_str());
  int number_bins = std::atoi(cmdLine.GetSafeArgument("-nb", 0, DEFAULT_BINS).c_str());
  int patch_s = std::atoi(cmdLine.GetSafeArgument("-ps", 0, DEFAULT_PATCH).c_str());
  float min_value = std::atof(cmdLine.GetSafeArgument("-pv", 0, DEFAULT_PIXEL_MIN).c_str());
  float max_value = std::atof(cmdLine.GetSafeArgument("-pv", 1, DEFAULT_PIXEL_MAX).c_str());
  // int dim = atoi(cmdLine.GetSafeArgument("-s", 0, DEFAULT_SIZE).c_str());

  // RescaleImage
  input_img = rescale_image(input_img, min_value, max_value);

  // Get the points in the image
  vector<InputImageType::IndexType> voxels = get_point_positions(input_points);

  // Define offsets
  using NeighborhoodType = itk::Neighborhood< typename InputImageType::PixelType, InputImageType::ImageDimension >;
  NeighborhoodType neighborhood;
  neighborhood.SetRadius( radius );

  unsigned int centerIndex = neighborhood.GetCenterNeighborhoodIndex();

  OffsetVector::Pointer offsets = OffsetVector::New();
  for( unsigned int d = 0; d < centerIndex; ++d )
  {
    offsets->push_back( neighborhood.GetOffset( d ) );
  }

  runlenght_texture_features(offsets, input_img, number_bins, patch_s, min_value, max_value, voxels);

  return EXIT_SUCCESS;
}

void runlenght_texture_features(
                                OffsetVector::Pointer offset,
                                InputImageType::Pointer inputImage,
                                int number_bins, int patch_size,
                                float pixel_min, float pixel_max,
                                vector<InputImageType::IndexType> positions
                                )
{
  int sx = inputImage->GetLargestPossibleRegion().GetSize(0);
  int sy = inputImage->GetLargestPossibleRegion().GetSize(1);
  int sz = inputImage->GetLargestPossibleRegion().GetSize(2);

  using Image2RunLengthType = itk::Statistics::ScalarImageToRunLengthMatrixFilter< InputImageType >;
  Image2RunLengthType::Pointer rlmGenerator = Image2RunLengthType::New();
  rlmGenerator->SetOffsets( offset );
  rlmGenerator->SetNumberOfBinsPerAxis( number_bins );
  rlmGenerator->SetPixelValueMinMax( pixel_min, pixel_max );

  using Hist2FeaturesType = itk::Statistics::HistogramToRunLengthFeaturesFilter< Image2RunLengthType::HistogramType >;
  Hist2FeaturesType::Pointer featureCalculator = Hist2FeaturesType::New();

  using ROIType = itk::RegionOfInterestImageFilter< InputImageType, InputImageType >;
  ROIType::Pointer roi = ROIType::New();
  roi->SetInput( inputImage );

  InputImageType::RegionType window;
  InputImageType::RegionType::SizeType size;

  std::cout << 100 << " " << 0 << " " << 8 << "\n";
  for (auto& voxels : positions) {
      size.Fill( patch_size );
      window.SetSize( size );
      window.SetIndex( voxels );

      roi->SetRegionOfInterest( window );
      roi->Update();
      rlmGenerator->SetInput( roi->GetOutput() );
      rlmGenerator->Update();

      featureCalculator->SetInput( rlmGenerator->GetOutput() );
      featureCalculator->Update();

      std::cout << voxels[0] << " " << voxels[1] << " " << voxels[2] << " ";

      std::cout<< std::setprecision(9) << featureCalculator->GetFeature(Hist2FeaturesType::ShortRunEmphasis) << " ";
      std::cout<< std::setprecision(9) << featureCalculator->GetFeature(Hist2FeaturesType::LongRunEmphasis) << " ";
      std::cout<< std::setprecision(9) << featureCalculator->GetFeature(Hist2FeaturesType::GreyLevelNonuniformity) << " ";
      std::cout<< std::setprecision(9) << featureCalculator->GetFeature(Hist2FeaturesType::RunLengthNonuniformity) << " ";
      std::cout<< std::setprecision(9) << featureCalculator->GetFeature(Hist2FeaturesType::LowGreyLevelRunEmphasis) << " ";
      std::cout<< std::setprecision(9) << featureCalculator->GetFeature(Hist2FeaturesType::HighGreyLevelRunEmphasis) << " ";
      std::cout<< std::setprecision(9) << featureCalculator->GetFeature(Hist2FeaturesType::ShortRunLowGreyLevelEmphasis) << " ";
      std::cout<< std::setprecision(9) << featureCalculator->GetFeature(Hist2FeaturesType::ShortRunHighGreyLevelEmphasis) << " ";
      std::cout<< std::setprecision(9) << featureCalculator->GetFeature(Hist2FeaturesType::LongRunLowGreyLevelEmphasis) << " ";
      std::cout<< std::setprecision(9) << featureCalculator->GetFeature(Hist2FeaturesType::LongRunHighGreyLevelEmphasis) << "\n";
  }
}

InputImageType::Pointer rescale_image(
                        InputImageType::Pointer inputImage,
                        float pixel_min, float pixel_max
                        )
{
  typedef itk::RescaleIntensityImageFilter< InputImageType, InputImageType > RescaleFilterType;
  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput(inputImage);
  rescaleFilter->SetOutputMinimum(pixel_min);
  rescaleFilter->SetOutputMaximum(pixel_max);
  rescaleFilter->Update();
  return rescaleFilter->GetOutput();
}

std::string remove_extension(std::string filename)
{
  size_t last_dot = filename.find_last_of(".");
  return last_dot == std::string::npos //npos is reached when find_last_of() doesnt find the symbol
                   ? filename // if filename has no extension, return itself
                   : filename.substr(0, last_dot); // return filename without extesion
}

vector<InputImageType::IndexType> get_point_positions(
                                                      std::string landmarks_parameters_filename
                                                     )
{
  // add patch diference position
  std::ifstream landmark_points(landmarks_parameters_filename.c_str());
  double parameters[6];
  unsigned int count, ind = 0;
  std::string line, value;
  vector<InputImageType::IndexType> voxels;
  InputImageType::IndexType posic;
  int strpos, strpos1;
  while(landmark_points.is_open())
  {
    getline(landmark_points, line);
    while(getline(landmark_points, line)){
      strpos = line.find(" ");
      strpos1 = line.find(" ", strpos+1);
      posic[0] =  atof(line.substr(0, strpos).c_str());
      posic[1] =  atof(line.substr(strpos+1, strpos1-(strpos+1)).c_str());
      strpos = line.find(" ", strpos1+1);
      posic[2] =  atof(line.substr(strpos1+1, strpos-(strpos1+1)).c_str());
      // std::cout << posic;
      voxels.push_back(posic);
    }
    landmark_points.close();
  }
  return voxels;
}