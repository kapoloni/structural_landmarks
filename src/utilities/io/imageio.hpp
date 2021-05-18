/**
 * @file   imageio.hpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup io
 * @ingroup    io
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#ifndef IMAGEIO_HPP
#define IMAGEIO_HPP

#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "assert.hpp"
#include "triple.hpp"


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
typename itk::Image<TPixel, Dimensions>::Pointer
read_image(std::string filename);

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
write_image(std::string filename, typename itk::Image<TPixel, Dimensions>::Pointer image);

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
TPixel*
image2array(typename itk::Image<TPixel, Dimensions>::Pointer image);

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
typename itk::Image<TPixel, Dimensions>::Pointer
array2image(TPixel *array, bip::triple<size_t> sizes,
            typename itk::Image<TPixel, Dimensions>::Pointer reference = nullptr);


}

#include "imageio.cpp"

#endif
