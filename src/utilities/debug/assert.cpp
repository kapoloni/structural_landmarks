/**
 * @file   assert.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup debug
 * @ingroup    debug
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#include "assert.hpp"


namespace bip
{
namespace debug
{


void
_assert(bool passed, const char *assertion, const char *file, long line)
{
    if (!passed) {
        std::stringstream ss;
        ss << "Failed assertion " << assertion << " in file " << file << " at line " << line;

        std::string message(ss.str());
        std::cerr << message << std::endl;

        throw std::logic_error(message);
    }
}


}
}
