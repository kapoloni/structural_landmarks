/**
 * @file   landmarkio.hpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup landmark-io
 * @ingroup    landmark-io
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#ifndef LANDMARKIO_HPP
#define LANDMARKIO_HPP

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "assert.hpp"
#include "triple.hpp"
#include "landmark.hpp"

// #define BIP_DEBUG_MODE
// #define BIP_VERBOSE_MODE


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
std::vector<bip::landmark> read_landmarks(std::string filename);


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
void write_landmarks(std::string filename, const std::vector<bip::landmark> &landmarks);


}

#endif
