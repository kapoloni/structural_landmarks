
#ifndef __bipMask_h
#define __bipMask_h

/*=========================================================================
 *
 *  Copyright David Doria 2011 daviddoria@gmail.com
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

/*
 * This class is a subclass of itkImage that provides the concept of "valid" pixels
 * and "hole" pixels. Pixels that are any other value are never used in computations.
 */

#include "itkImage.h"
#include "itkImageRegionIterator.h"

namespace bip
{

namespace mask
{

class Mask : public itk::Image< unsigned char, 2>
{
public:
  /** Standard typedefs. */
  typedef Mask                                   Self;
  typedef itk::Image< unsigned char, 2>          Superclass;
  typedef itk::SmartPointer< Self >              Pointer;
  typedef itk::SmartPointer< const Self >        ConstPointer;
  typedef itk::WeakPointer< const Self >         ConstWeakPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(Mask, Image);

  /** Dimension of the image. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Types derived from the Superclass */
  typedef typename Superclass::IndexType IndexType;
  typedef typename Superclass::IOPixelType IOPixelType;

  /** Tyepdef for the functor used to access a neighborhood of pixel
  * pointers. */
  typedef itk::NeighborhoodAccessorFunctor< Self >  NeighborhoodAccessorFunctorType;

  /** Return the NeighborhoodAccessor functor. This method is called by the
  * neighborhood iterators. */
  NeighborhoodAccessorFunctorType GetNeighborhoodAccessor()
  { 
    return NeighborhoodAccessorFunctorType(); 
  }

  /** Return the NeighborhoodAccessor functor. This method is called by the
  * neighborhood iterators. */
  const NeighborhoodAccessorFunctorType GetNeighborhoodAccessor() const
  { 
    return NeighborhoodAccessorFunctorType(); 
  }

  bool IsHole(itk::Index<2> index) const
  {
    if(this->GetPixel(index) == this->HoleValue)
    {
      return true;
    }
    return false;
  }

  bool IsValid(const itk::ImageRegion<2> region) const
  {
    // If any of the pixels in the region are invalid, the region is invalid.
    itk::ImageRegionConstIterator<Mask> maskIterator(this, region);

    while(!maskIterator.IsAtEnd())
    {
      if(!this->IsValid(maskIterator.GetIndex()))
      {
      //std::cout << "Mask::IsValid - Pixel " << maskIterator.GetIndex() << " has value " << static_cast<unsigned int>(maskIterator.Get()) 
      //          << " which makes the region invalid because Mask::ValidValue = " << static_cast<unsigned int>(this->ValidValue) << std::endl;
      return false;
      }

      ++maskIterator;
    }
    return true;
  }

  bool IsValid(itk::Index<2> index) const
  {
    if(this->GetPixel(index) == this->ValidValue)
    {
      return true;
    }
    return false;
  }

  void Invert()
  {
    // Exchange HoleValue and ValidValue, but leave everything else alone.
    itk::ImageRegionIterator<Mask> maskIterator(this, this->GetLargestPossibleRegion());
    unsigned int invertedCounter = 0;
    while(!maskIterator.IsAtEnd())
    {
      if(this->IsValid(maskIterator.GetIndex()))
      {
        maskIterator.Set(this->HoleValue);
        invertedCounter++;
      }
      else if(this->IsHole(maskIterator.GetIndex()))
      {
        maskIterator.Set(this->ValidValue);
        invertedCounter++;
      }
      ++maskIterator;
    }
    std::cout << "Inverted " << invertedCounter << " in the mask." << std::endl;
  }

  void Cleanup()
  {
    // We want to interpret pixels that are "pretty much hole value" as holes, and pixels that
    // are "pretty much valid value" as valid. The "do not use" pixels must be very far away from both of these values.
    itk::ImageRegionIterator<Mask> maskIterator(this, this->GetLargestPossibleRegion());

    float tolerance = 4;
    while(!maskIterator.IsAtEnd())
    {
      if(fabs(maskIterator.Get() - this->ValidValue) < tolerance)
      {
        //std::cout << "Setting valid pixel to " << static_cast<unsigned int>(this->ValidValue) << std::endl;
        maskIterator.Set(this->ValidValue);
      }
      else if(fabs(maskIterator.Get() - this->HoleValue) < tolerance)
      {
        //std::cout << "Setting hole pixel to " << static_cast<unsigned int>(this->HoleValue) << std::endl;
        maskIterator.Set(this->HoleValue);
      }
      ++maskIterator;
    }
  }

  void SetHoleValue(unsigned char value)
  {
    this->HoleValue = value; 
  }

  void SetValidValue(unsigned char value)
  {
    this->ValidValue = value; 
  }

  unsigned char GetHoleValue() const
  {
    return this->HoleValue;
  }

  unsigned char GetValidValue() const
  {
    return this->ValidValue;
  }

  void OutputMembers() const
  {
    std::cout << "HoleValue: " << static_cast<unsigned int>(this->HoleValue) << std::endl;
    std::cout << "ValidValue: " << static_cast<unsigned int>(this->ValidValue) << std::endl;
  }

  void DeepCopyFrom(Mask::Pointer inputMask)
  {
    this->SetRegions(inputMask->GetLargestPossibleRegion());
    this->Allocate();

    itk::ImageRegionConstIterator<Mask> inputIterator(inputMask, inputMask->GetLargestPossibleRegion());
    itk::ImageRegionIterator<Mask> thisIterator(this, this->GetLargestPossibleRegion());

    while(!inputIterator.IsAtEnd())
    {
      thisIterator.Set(inputIterator.Get());
      ++inputIterator;
      ++thisIterator;
    }
    this->SetHoleValue(inputMask->GetHoleValue());
    this->SetValidValue(inputMask->GetValidValue());
  }

protected:
  Mask()
  {
    this->HoleValue = 255;
    this->ValidValue = 0;
  }

  unsigned char HoleValue; // Pixels with this value will be filled.
  unsigned char ValidValue; // Pixels with this value will not be filled - they are the source region.

private:
  Mask(const Self &);    //purposely not implemented
  void operator=(const Self &); //purposely not implemented
};


/**

*/
itk::ImageRegion<2> 
GetRegionInRadiusAroundPixel( const itk::Index<2> pixel, const unsigned int radius )
{
    // This function returns a Region with the specified 'radius' centered at 'pixel'. 
  //By the definition of the radius of a square patch, the output region is (radius*2 + 1)x(radius*2 + 1).
    // Note: This region is not necessarily entirely inside the image!
  
    // The "index" is the lower left corner, so we need to subtract the radius from the center to obtain it
    itk::Index<2> lowerLeft;
    lowerLeft[0] = pixel[0] - radius;
    lowerLeft[1] = pixel[1] - radius;

    itk::ImageRegion<2> region;
    region.SetIndex(lowerLeft);
    itk::Size<2> size;
    size[0] = radius*2 + 1;
    size[1] = radius*2 + 1;
    region.SetSize(size);

    return region;
}

/**

*/
itk::Index<2> 
GetRegionCenter( const itk::ImageRegion<2> region )
{
    itk::Index<2> center;
    center[0] = region.GetIndex()[0] + region.GetSize()[0] / 2;
    center[1] = region.GetIndex()[1] + region.GetSize()[1] / 2;

    return center;
}



} // End namespace utils

} // End namespace bip

#endif
