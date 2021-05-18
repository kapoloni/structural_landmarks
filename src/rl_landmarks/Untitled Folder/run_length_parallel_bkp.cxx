#include <stdio.h>
#include <stdlib.h>
#include <CmdLine.h>
#include <bipUtils.h>
#include <omp.h>
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
typedef itk::Statistics::ScalarImageToRunLengthMatrixFilter<InternalImageType>  Image2CoOccuranceType;
typedef Image2CoOccuranceType::HistogramType HistogramType;
typedef itk::Statistics::HistogramToRunLengthFeaturesFilter<HistogramType> Hist2FeaturesType;
typedef InternalImageType::OffsetType OffsetType;
typedef itk::AddImageFilter <InternalImageType> AddImageFilterType;
typedef itk::MultiplyImageFilter<InternalImageType> MultiplyImageFilterType;

// Parse command line arguments
CCmdLine cmd;


//calculate features for one offset
void calcTextureFeatureImage (OffsetType offset,
    InternalImageType::Pointer inputImage, InternalImageType::Pointer outInertia,
    InternalImageType::Pointer outCorrelation, InternalImageType::Pointer outEnergy)
{
	#pragma omp parallel
	{

	
    //allocate output images
    outInertia->CopyInformation(inputImage);
    outInertia->SetRegions(inputImage->GetLargestPossibleRegion());
    outInertia->Allocate();
    outInertia->FillBuffer(0);
    
    outCorrelation->CopyInformation(inputImage);
    outCorrelation->SetRegions(inputImage->GetLargestPossibleRegion());
    outCorrelation->Allocate();
    outCorrelation->FillBuffer(0);
    
    outEnergy->CopyInformation(inputImage);
    outEnergy->SetRegions(inputImage->GetLargestPossibleRegion());
    outEnergy->Allocate();
    outEnergy->FillBuffer(0);

    InternalImageType::Pointer mask;
    if(cmd.HasSwitch("-m"))
    	mask = bip::utils::ReadImage<InternalImageType>(cmd.GetArgument("-m", 0));
    
    Image2CoOccuranceType::Pointer glcmGenerator=Image2CoOccuranceType::New();
    
    glcmGenerator->SetOffset(offset);
    glcmGenerator->SetNumberOfBinsPerAxis(16); //reasonable number of bins
    glcmGenerator->SetPixelValueMinMax(0, 255); //for input UCHAR pixel type
    Hist2FeaturesType::Pointer featureCalc = Hist2FeaturesType::New();
    typedef itk::RegionOfInterestImageFilter<InternalImageType,InternalImageType> roiType;
    roiType::Pointer roi=roiType::New();
    roi->SetInput(inputImage);
    InternalImageType::RegionType window;
    InternalImageType::RegionType::SizeType size;
    size.Fill(3); //window size=3x3x3
    window.SetSize(size);
    InternalImageType::IndexType pi; //pixel index
    
    //slide window over the entire image
    //#pragma omp parallel for private(pi, window, roi, glcmGenerator, featureCalc) num_threads(8)
    for (unsigned x = 1; x < inputImage->GetLargestPossibleRegion().GetSize(0)-1; ++x)
    {
        pi.SetElement(0, x);
        window.SetIndex(0, x - 1);
        for (unsigned y = 1; y < inputImage->GetLargestPossibleRegion().GetSize(1) - 1; ++y)
        {
            pi.SetElement(1, y);
            window.SetIndex(1, y - 1);
            for (unsigned z = 1; z < inputImage->GetLargestPossibleRegion().GetSize(2) - 1; ++z)
            {
                pi.SetElement(2, z);
                window.SetIndex(2, z - 1);
                roi->SetRegionOfInterest(window);
                roi->Update();

                glcmGenerator->SetInput(roi->GetOutput());
                glcmGenerator->Update();
                /*if(cmd.HasSwitch("-m"))
                	glcmGenerator->SetMaskImage(mask);*/

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

                featureCalc->SetInput( glcmGenerator->GetOutput() );
                featureCalc->Update();
                outInertia->SetPixel(pi, featureCalc->GetFeature(Hist2FeaturesType::LongRunHighGreyLevelEmphasis));
                //outCorrelation->SetPixel(pi, featureCalc->GetFeature(Hist2FeaturesType::ShortRunHighGreyLevelEmphasis));
                //outEnergy->SetPixel(pi, featureCalc->GetFeature(Hist2FeaturesType::LongRunLowGreyLevelEmphasis));
            }
        }
        std::cout<< "Done for the " << x << "th voxel.\n";
    }
    std::cout << "\n";
	}
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
     	InternalImageType::Pointer inertia = InternalImageType::New();
     	InternalImageType::Pointer correlation = InternalImageType::New();
     	InternalImageType::Pointer energy = InternalImageType::New();
     	calcTextureFeatureImage(offset, image, inertia, correlation, energy);
	     
     	writer->SetInput(inertia);
	   	std::stringstream ssInertia;
     	ssInertia << "LongRunHigh" << d << ".nii";
     	writer->SetFileName(ssInertia.str());
     	writer->Update();
     	
     	/*writer->SetInput(correlation);
     	std::stringstream ssCorrelation;
     	ssCorrelation << "ShortRunHigh" << d << ".nii";
     	writer->SetFileName(ssCorrelation.str());
     	writer->Update();
     	
     	writer->SetInput(energy);
     	std::stringstream ssEnergy;
     	ssEnergy << "LongRunLow" << d << ".nii";
     	writer->SetFileName(ssEnergy.str());
     	writer->Update();*/
     	
     	std::cout<<'\n';
  }

	return 0;
}