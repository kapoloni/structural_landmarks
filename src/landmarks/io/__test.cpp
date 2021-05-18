/**
 * @file   __test.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup landmark-io
 * @ingroup    landmark-io
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

// #define BIP_DEBUG_MODE
// #define BIP_VERBOSE_MODE

#include <cstdlib>
#include <iostream>
#include <vector>
#include "landmark.hpp"
#include "landmarkio.hpp"


int main(int argc, char *argv[])
{
    std::vector<bip::landmark> landmarks;
    try
    {
        landmarks = bip::read_landmarks("data/landmarks.txt");        
    }
    catch (const char *exception)
    {
        std::cerr << exception << std::endl;
    }

    std::cout << "Num landmarks = " << landmarks.size() << "\n";

    for (size_t i = 0; i < landmarks.size(); ++i)
        std::cout << landmarks[i] << "\n";
    std::cout << std::endl;


    return EXIT_SUCCESS;
}
