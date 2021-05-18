#include <stdio.h>
#include <stdlib.h>
#include <CmdLine.h>
#include <bipUtils.h>
#include <omp.h>
#include <string.h>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <itkDenseFrequencyContainer2.h>

#include "itkHistogramToRunLengthFeaturesFilter.h"
#include "itkScalarImageToRunLengthMatrixFilter.h"

#include "itkVectorContainer.h"
#include "itkAddImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"

//definitions of used types
typedef itk::Image<unsigned char, 3> InternalImageType; // colocar unsigned char?
typedef itk::Image<int, 3> MaskImageType;
typedef itk::Image<unsigned char, 3> VisualizingImageType;
typedef itk::Neighborhood<float, 3> NeighborhoodType;
typedef itk::Statistics::ScalarImageToRunLengthMatrixFilter<InternalImageType>  Image2RunLengthType;
typedef Image2RunLengthType::HistogramType HistogramType;
typedef itk::Statistics::HistogramToRunLengthFeaturesFilter<HistogramType> Hist2FeaturesType;
typedef InternalImageType::OffsetType OffsetType;
typedef itk::AddImageFilter<InternalImageType> AddImageFilterType;
typedef itk::MultiplyImageFilter<InternalImageType> MultiplyImageFilterType;

// Parse command line arguments
CCmdLine cmd;

//calculate features for one offset
void calcTextureFeatureImage (OffsetType offset,
                  InternalImageType::Pointer inputImage){

  Image2RunLengthType::Pointer rlmGenerator = Image2RunLengthType::New();
  rlmGenerator->SetOffset(offset);
  rlmGenerator->SetNumberOfBinsPerAxis(16);
  rlmGenerator->SetPixelValueMinMax(0, 63);

  Hist2FeaturesType::Pointer featureCalc = Hist2FeaturesType::New();

  typedef itk::RegionOfInterestImageFilter< TImageType, TImageType > ROIType;
  ROIType::Pointer roi = ROIType::New();
  roi->SetInput( inputImage );

  TImageType::RegionType window_patch;
  TImageType::RegionType::SizeType size;
  size.Fill( patch_size );
  window_patch.SetSize( size );
  std::cout<< "Inertia\t" << "Energy\t" << "Entropy\t" << "Correlation\t" << "InverseDifferenceMoment\t"
           << "ClusterShade\t" << "ClusterProminence\t" << "HaralickCorrelation\n";
  for (auto l : voxels) {
    window_patch.SetIndex(l);
    roi->SetRegionOfInterest( window_patch );
    roi->Update();
    rlmGenerator->SetInput(roi->GetOutput());
    rlmGenerator->Update();
    featureCalc->SetInput(rlmGenerator->GetOutput());
    featureCalc->Update();
    outShortRunEmphasisstd::cout<< featureCalc->GetFeature(Hist2FeaturesType::ShortRunEmphasis));
    outLongRunEmphasisstd::cout<< featureCalc->GetFeature(Hist2FeaturesType::LongRunEmphasis));
    outGreyLevelNonuniformitystd::cout<< featureCalc->GetFeature(Hist2FeaturesType::GreyLevelNonuniformity));
    outRunLengthNonuniformitystd::cout<< featureCalc->GetFeature(Hist2FeaturesType::RunLengthNonuniformity));
    outLowGreyLevelRunEmphasisstd::cout<< featureCalc->GetFeature(Hist2FeaturesType::LowGreyLevelRunEmphasis));
    outHighGreyLevelRunEmphasisstd::cout<< featureCalc->GetFeature(Hist2FeaturesType::HighGreyLevelRunEmphasis));
    outShortRunLowGreyLevelEmphasisstd::cout<< featureCalc->GetFeature(Hist2FeaturesType::ShortRunLowGreyLevelEmphasis));
    outShortRunHighGreyLevelEmphasisstd::cout<< featureCalc->GetFeature(Hist2FeaturesType::ShortRunHighGreyLevelEmphasis));
    std::cout<< featureCalc->GetFeature(Hist2FeaturesType::LongRunLowGreyLevelEmphasis));
    std::cout<< featureCalc->GetFeature(Hist2FeaturesType::LongRunHighGreyLevelEmphasis) << "\t";

  }


              /*
              List of possible features:
              - Hist2FeaturesType::ShortRunEmphasis
              - Hist2FeaturesType::LongRunEmphasis
              - Hist2FeaturesType::GreyLevelNonuniformity
              - Hist2FeaturesType::RunLengthNonuniformity
              - Hist2FeaturesType::LowGreyLevelRunEmphasis
              - Hist2FeaturesType::HighGreyLevelRunEmphasis
              - Hist2FeaturesType::ShortRunLowGreyLevelEmphasis
              - Hist2FeaturesType::ShortRunHighGreyLevelEmphasis
              - Hist2FeaturesType::LongRunLowGreyLevelEmphasis
              - Hist2FeaturesType::LongRunHighGreyLevelEmphasis
              */

              featureCalc->SetInput(rlmGenerator->GetOutput());
              featureCalc->Update();
              std::cout<< featureCalc->GetFeature(Hist2FeaturesType::ShortRunEmphasis) << "\t";
              std::cout<< featureCalc->GetFeature(Hist2FeaturesType::LongRunEmphasis) << "\t";
              std::cout<< featureCalc->GetFeature(Hist2FeaturesType::GreyLevelNonuniformity) << "\t";
              std::cout<< featureCalc->GetFeature(Hist2FeaturesType::RunLengthNonuniformity) << "\t";
              std::cout<< featureCalc->GetFeature(Hist2FeaturesType::LowGreyLevelRunEmphasis) << "\t";
              std::cout<< featureCalc->GetFeature(Hist2FeaturesType::HighGreyLevelRunEmphasis) << "\t";
              std::cout<< featureCalc->GetFeature(Hist2FeaturesType::ShortRunLowGreyLevelEmphasis) << "\t";
              std::cout<< featureCalc->GetFeature(Hist2FeaturesType::ShortRunHighGreyLevelEmphasis) << "\t";
              std::cout<< featureCalc->GetFeature(Hist2FeaturesType::LongRunLowGreyLevelEmphasis) << "\t";
              std::cout<< featureCalc->GetFeature(Hist2FeaturesType::LongRunHighGreyLevelEmphasis) << "\n";
          }
      }
      std::cout<< "Done for the " << x << "th voxel.\n";
  }
  std::cout << "\n";
}

void showUsage(std::string str)
{
  std::cout << "\n";
  std::cout << "----------------------\n"
        << "Usage: " << str << "\n"
        << "-i inputImage\n"
        << "[-m maskImage]\n"
        << "----------------------\n\n";
}

int main(int argc, char **argv)
{
  // Type definitions
  typedef double InputPixelType; // change it according to your needs
  const int Dimension = 3;

  typedef itk::Image<InputPixelType, Dimension> TScalarInputImage; // input image type

  // SplitLine returns the number of switches
  // If we have input (-i) and output (-o) images
  // then we must have two switches.
  if(cmd.SplitLine(argc, argv) < 1) // change this according to your needs
  {
    showUsage(argv[0]);
    return -1;
  }

  int nInputImages = cmd.GetArgumentCount("-i"); // get number of input images;

  if(nInputImages != 1)
  {
    std::cout << "Please input only ONE image.\n";
    return -1;
  }

  // Get filenames
  StringType inputFileName = cmd.GetArgument("-i", 0);
  char buffer[50] = "rlm_";

  // Get the filename without extension
  int i = 0, j = 4;
  for(i = 0; inputFileName[i] != '.'; ++i)
    buffer[j++] = inputFileName[i];
  buffer[j] = '_';

  // Read the input image
  //TScalarInputImage::Pointer inputImg = bip::utils::ReadImage<TScalarInputImage>(inputFileName);
  //InternalImageType::Pointer image = bip::utils::ReadImage<InternalImageType>(inputFileName);

    NeighborhoodType neighborhood;
    neighborhood.SetRadius(1);
    unsigned int centerIndex = neighborhood.GetCenterNeighborhoodIndex();
    OffsetType offset;
    typedef itk::ImageFileWriter<InternalImageType> WriterType;
    //WriterType::Pointer writer=WriterType::New();

    #pragma omp parallel for private(offset) num_threads(8)
    for ( unsigned int d = 0; d < centerIndex; d++ )
    {
        InternalImageType::Pointer image = bip::utils::ReadImage<InternalImageType>(inputFileName);
        WriterType::Pointer writer=WriterType::New();

        offset = neighborhood.GetOffset(d);
        InternalImageType::Pointer short_run_emphasis = InternalImageType::New();
        InternalImageType::Pointer long_run_emphasis = InternalImageType::New();
        InternalImageType::Pointer greylevel_non_uniformity = InternalImageType::New();
        InternalImageType::Pointer run_length_non_uniformity = InternalImageType::New();
        InternalImageType::Pointer low_greylevel_run_emphasis = InternalImageType::New();
        InternalImageType::Pointer high_grey_level_run_emphais = InternalImageType::New();
        InternalImageType::Pointer short_run_low_greylevel_emphasis = InternalImageType::New();
        InternalImageType::Pointer short_run_high_greylevel_emphasis = InternalImageType::New();
        InternalImageType::Pointer long_run_low_greylevel_emphasis = InternalImageType::New();
        InternalImageType::Pointer long_run_high_greylevel_emphasis = InternalImageType::New();

        calcTextureFeatureImage(offset,
                      image,
                    short_run_emphasis,
                    long_run_emphasis,
                    greylevel_non_uniformity,
                    run_length_non_uniformity,
                    low_greylevel_run_emphasis,
                    high_grey_level_run_emphais,
                    short_run_low_greylevel_emphasis,
                    short_run_high_greylevel_emphasis,
                    long_run_low_greylevel_emphasis,
                    long_run_high_greylevel_emphasis);

        writer->SetInput(short_run_emphasis);
        std::stringstream sshort_run_emphasis;
        sshort_run_emphasis << buffer << "ShortRunEmphasis_" << d << ".nii";
        writer->SetFileName(sshort_run_emphasis.str());
        writer->Update();

        writer->SetInput(long_run_emphasis);
        std::stringstream slong_run_emphasis;
        slong_run_emphasis << buffer << "LongRunEmphasis_" << d << ".nii";
        writer->SetFileName(slong_run_emphasis.str());
        writer->Update();

        writer->SetInput(greylevel_non_uniformity);
        std::stringstream sgreylevel_non_uniformity;
        sgreylevel_non_uniformity << buffer << "GreyLevelNonUniformity_" << d << ".nii";
        writer->SetFileName(sgreylevel_non_uniformity.str());
        writer->Update();

        writer->SetInput(run_length_non_uniformity);
        std::stringstream srun_length_non_uniformity;
        srun_length_non_uniformity << buffer << "RunLengthNonUniformity_" << d << ".nii";
        writer->SetFileName(srun_length_non_uniformity.str());
        writer->Update();

        writer->SetInput(low_greylevel_run_emphasis);
        std::stringstream slow_greylevel_run_emphasis;
        slow_greylevel_run_emphasis << buffer << "LowGreyLevelRunEmphasis_" << d << ".nii";
        writer->SetFileName(slow_greylevel_run_emphasis.str());
        writer->Update();

        writer->SetInput(high_grey_level_run_emphais);
        std::stringstream shigh_grey_level_run_emphais;
        shigh_grey_level_run_emphais << buffer << "HighGreyLevelRunEmphasis_" << d << ".nii";
        writer->SetFileName(shigh_grey_level_run_emphais.str());
        writer->Update();

        writer->SetInput(short_run_low_greylevel_emphasis);
        std::stringstream sshort_run_low_greylevel_emphasis;
        sshort_run_low_greylevel_emphasis << buffer << "ShortRunLowGreyLevelRunEmphasis_" << d << ".nii";
        writer->SetFileName(sshort_run_low_greylevel_emphasis.str());
        writer->Update();

        writer->SetInput(short_run_high_greylevel_emphasis);
        std::stringstream sshort_run_high_greylevel_emphasis;
        sshort_run_high_greylevel_emphasis << buffer << "ShortRunHighGreyLevelRunEmphasis_" << d << ".nii";
        writer->SetFileName(sshort_run_high_greylevel_emphasis.str());
        writer->Update();

        writer->SetInput(long_run_low_greylevel_emphasis);
        std::stringstream slong_run_low_greylevel_emphasis;
        slong_run_low_greylevel_emphasis << buffer << "LongRunLowGreyLevelEmphasis_" << d << ".nii";
        writer->SetFileName(slong_run_low_greylevel_emphasis.str());
        writer->Update();

        writer->SetInput(long_run_high_greylevel_emphasis);
        std::stringstream slong_run_high_greylevel_emphasis;
        slong_run_high_greylevel_emphasis << buffer << "LongRunHighGreyLevelEmphasis_" << d << ".nii";
        writer->SetFileName(slong_run_high_greylevel_emphasis.str());
        writer->Update();

      std::cout<<'\n';
    }

  return 0;
}
