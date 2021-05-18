/**
 * @file   __test.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup debug
 * @ingroup    debug
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include "assert.hpp"


int main(int argc, char *argv[])
{
    bip::debug::assert2(10 < 20); // Nothing happens.
    bip::debug::assert2(10 > 20); // Throws std::logic_error.

    return EXIT_SUCCESS;
}
