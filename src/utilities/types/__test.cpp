/**
 * @file   __test.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup types
 * @ingroup    types
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include "triple.hpp"
#include "landmark.hpp"


int main(int argc, char *argv[])
{
    // Constructor.
    bip::triple<size_t> t0;
    bip::triple<size_t> t1(10);
    bip::triple<size_t> t2(10, 20);
    bip::triple<size_t> t3(10, 20, 30);

    // Assignment.
    bip::triple<size_t> t4 = t2;

    // Access of elements.
    t4[0] = t4[1] *= t2[0];
    t4[2] = t3[2] * 10;

    // Stream output.
    std::cout << "t0 = " << t0 << std::endl; // (0, 0, 0)
    std::cout << "t1 = " << t1 << std::endl; // (10, 0, 0)
    std::cout << "t2 = " << t2 << std::endl; // (10, 20, 0)
    std::cout << "t3 = " << t3 << std::endl; // (10, 20, 30)
    std::cout << "t4 = " << t4 << std::endl; // (200, 200, 300)

    // --------------------------------------------------------------------------------------------

    bip::triple<size_t> location(70, 80, 90);
    bip::triple<float>  features(0.75, 0.5236, 1.0472);
    std::vector<float>  local_descriptor;
    std::vector<float>  global_descriptor;

    for (size_t i = 1; i < 10; ++i) {
        local_descriptor.push_back(i * 0.5);
        global_descriptor.push_back(i * 0.25);
    }

    // Constructor.
    bip::landmark l0;
    bip::landmark l1(location);
    bip::landmark l2(location, features);
    bip::landmark l3(location, features, local_descriptor);
    bip::landmark l4(location, features, local_descriptor, global_descriptor);

    // Stream output.
    std::cout << "l0 = " << l0 << std::endl;
    std::cout << "l1 = " << l1 << std::endl;
    std::cout << "l2 = " << l2 << std::endl;
    std::cout << "l3 = " << l3 << std::endl;
    std::cout << "l4 = " << l4 << std::endl;

    return EXIT_SUCCESS;
}
