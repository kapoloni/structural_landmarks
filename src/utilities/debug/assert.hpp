/**
 * @file   assert.hpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup debug
 * @ingroup    debug
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#ifndef ASSERT_HPP
#define ASSERT_HPP

#include <cstddef>
#include <cassert>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>


namespace bip
{
namespace debug
{


/**
 * @fn {definition-name}
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
void _assert(bool passed, const char *assertion, const char *file, long line);

/**
 * @def {definition-name}
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
// #ifdef NDEBUG
    // #define assert2(expr) _assert(true, "", "", 0)
// #else
    #define assert2(expr) _assert(expr, #expr, __FILE__, __LINE__)
// #endif


}
}

#endif
