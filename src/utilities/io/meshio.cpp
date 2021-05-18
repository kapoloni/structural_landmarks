/**
 * @file   meshio.cpp
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
typename itk::Mesh<TPixel, Dimensions>::Pointer
read_mesh(std::string filename)
{
    bip::debug::assert2(!filename.empty());

    typedef itk::Mesh<TPixel, Dimensions> TMesh;
    typedef itk::MeshFileReader<TMesh>    TReader;

    typename TReader::Pointer reader = TReader::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    return reader->GetOutput();
}


template <class TPixel, size_t Dimensions>
void
write_mesh(std::string filename, typename itk::Mesh<TPixel, Dimensions>::Pointer mesh)
{
    bip::debug::assert2(!filename.empty());
    bip::debug::assert2(mesh.IsNotNull());

    typedef itk::Mesh<TPixel, Dimensions> TMesh;
    typedef itk::MeshFileWriter<TMesh>    TWriter;

    typename TWriter::Pointer writer = TWriter::New();
    writer->SetInput(mesh);
    writer->SetFileName(filename.c_str());
    writer->Update();
}


}
