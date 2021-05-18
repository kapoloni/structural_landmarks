/**
 * @file   meshio.hpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup io
 * @ingroup    io
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#ifndef MESHIO_HPP
#define MESHIO_HPP

#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <itkMesh.h>
#include <itkMeshFileReader.h>
#include <itkMeshFileWriter.h>
#include "assert.hpp"


namespace bip
{


/**
 * @fn {function-name}
 *
 * @brief {brief-description}
 * {detailed-description}
 * 
 * @param[in]  {param-name} {param-description}
 * @param[out] {param-name} {param-description}
 *
 * @returns {return-value-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
template <class TPixel, size_t Dimensions>
typename itk::Mesh<TPixel, Dimensions>::Pointer
read_mesh(std::string filename);

/**
 * @fn {function-name}
 *
 * @brief {brief-description}
 * {detailed-description}
 * 
 * @param[in]  {param-name} {param-description}
 * @param[out] {param-name} {param-description}
 *
 * @returns {return-value-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
template <class TPixel, size_t Dimensions>
void
write_mesh(std::string filename, typename itk::Mesh<TPixel, Dimensions>::Pointer mesh);


}

#include "meshio.cpp"

#endif
