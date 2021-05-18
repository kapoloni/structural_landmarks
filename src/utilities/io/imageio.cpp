/**
 * @file   imageio.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup io
 * @ingroup    io
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

namespace bip
{


template <class TPixel, size_t Dimensions>
typename itk::Image<TPixel, Dimensions>::Pointer
read_image(std::string filename)
{
    bip::debug::assert2(!filename.empty());

    typedef itk::Image<TPixel, Dimensions> TImage;
    typedef itk::ImageFileReader<TImage>   TReader;

    typename TReader::Pointer reader = TReader::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    return reader->GetOutput();
}


template <class TPixel, size_t Dimensions>
void
write_image(std::string filename, typename itk::Image<TPixel, Dimensions>::Pointer image)
{
    bip::debug::assert2(!filename.empty());
    bip::debug::assert2(image.IsNotNull());

    typedef itk::Image<TPixel, Dimensions> TImage;
    typedef itk::ImageFileWriter<TImage>   TWriter;

    typename TWriter::Pointer writer = TWriter::New();
    writer->SetInput(image);
    writer->SetFileName(filename.c_str());
    writer->Update();
}


template <class TPixel, size_t Dimensions>
TPixel*
image2array(typename itk::Image<TPixel, Dimensions>::Pointer image)
{
    bip::debug::assert2(image.IsNotNull());

    typedef itk::Image<TPixel, Dimensions> TImage;

    typename TImage::SizeType sizes = image->GetBufferedRegion().GetSize();

    size_t total_size = 1;
    for (size_t d = 0; d < Dimensions; ++d)
        total_size *= sizes[d];

    TPixel *data = new TPixel[total_size]();

    itk::ImageRegionConstIterator<TImage> iterator(image, image->GetBufferedRegion());
    size_t i = 0;

    while (!iterator.IsAtEnd())
    {
        data[i] = iterator.Get();
        ++iterator;
        ++i;
    }

    return data;
}


template <class TPixel, size_t Dimensions>
typename itk::Image<TPixel, Dimensions>::Pointer
array2image(TPixel *data, bip::triple<size_t> sizes,
            typename itk::Image<TPixel, Dimensions>::Pointer reference)
{
    bip::debug::assert2(data != nullptr);

    typedef itk::Image<TPixel, Dimensions> TImage;

    typename TImage::SizeType sizes2;
    for (size_t d = 0; d < Dimensions; ++d)
        sizes2[d] = sizes[d];

    typename TImage::IndexType start;
    for (size_t d = 0; d < Dimensions; ++d)
        start[d] = 0;

    typename TImage::RegionType region;
    region.SetSize(sizes2);
    region.SetIndex(start);

    typename TImage::Pointer image = TImage::New();
    image->SetRegions(region);
    if (reference.IsNotNull()) {
        image->SetOrigin(reference->GetOrigin());
        image->SetSpacing(reference->GetSpacing());
        image->SetDirection(reference->GetDirection());
    }
    image->Allocate();
    image->FillBuffer(TPixel());

    itk::ImageRegionIterator<TImage> iterator(image, image->GetBufferedRegion());
    size_t i = 0;

    while (!iterator.IsAtEnd())
    {
        iterator.Set(data[i]);
        ++iterator;
        ++i;
    }

    return image;
}


}
