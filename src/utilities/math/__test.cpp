/**
 * @file   __test.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup math
 * @ingroup    math
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#define _USE_MATH_DEFINES
#include <cstddef>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "triple.hpp"
#include "mathfunctions.hpp"


int main(int argc, char *argv[])
{
    std::cout << bip::rad2deg(0)      << std::endl; // 0
    std::cout << bip::rad2deg(M_PI_4) << std::endl; // 45
    std::cout << bip::rad2deg(M_PI_2) << std::endl; // 90
    std::cout << bip::rad2deg(M_PI)   << std::endl; // 180

    std::cout << bip::deg2rad(0)   << std::endl; // 0
    std::cout << bip::deg2rad(45)  << std::endl; // 0.785398
    std::cout << bip::deg2rad(90)  << std::endl; // 1.5708
    std::cout << bip::deg2rad(180) << std::endl; // 3.14159

    // --------------------------------------------------------------------------------------------

    std::cout << bip::sph2cart(bip::triple<float>(1, 0, 0))           << std::endl;
    std::cout << bip::sph2cart(bip::triple<float>(1, 0, M_PI_2))      << std::endl;
    std::cout << bip::sph2cart(bip::triple<float>(1, M_PI_2, M_PI_4)) << std::endl;
    std::cout << bip::sph2cart(bip::triple<float>(1, 0, -M_PI_2))     << std::endl;

    std::cout << bip::cart2sph(bip::triple<float>(0, 0, 0)) << std::endl;
    std::cout << bip::cart2sph(bip::triple<float>(1, 0, 0)) << std::endl;
    std::cout << bip::cart2sph(bip::triple<float>(0, 1, 0)) << std::endl;
    std::cout << bip::cart2sph(bip::triple<float>(0, 0, 1)) << std::endl;

    // --------------------------------------------------------------------------------------------

    // ...

    return EXIT_SUCCESS;
}
