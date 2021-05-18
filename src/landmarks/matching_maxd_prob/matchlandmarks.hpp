/**
 * @file   matchlandmarks.hpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup landmark-matching
 * @ingroup    landmark-matching
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#ifndef MATCHLANDMARKS_HPP
#define MATCHLANDMARKS_HPP

#define BIP_VERBOSE_MODE
// #define BIP_DEBUG_MODE

#define _USE_MATH_DEFINES
#include <cstddef>
#include <cfloat>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <algorithm>
#include "assert.hpp"
#include "mathfunctions.hpp"
#include "triple.hpp"
#include "landmark.hpp"


namespace bip
{


typedef std::vector<std::pair<size_t, size_t> > matches_vector;


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
matches_vector
match_landmarks(const std::vector<landmark> &landmarks1,
                const std::vector<landmark> &landmarks2,
                float                       descriptors_tradeoff = 0.5,
                float                       max_descriptor_dist  = FLT_MAX,
                float                       max_location_dist    = FLT_MAX,
	            float                       min_location_dist    = 0.0);
}

#endif
