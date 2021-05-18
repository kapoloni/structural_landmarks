/**
 * @file   __test.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup filterbank
 * @ingroup    filterbank
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

// #define BIP_DEBUG_MODE
// #define BIP_VERBOSE_MODE

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <string>
#include "triple.hpp"
#include "loggaborfilterbank.hpp"


int main(int argc, char const *argv[])
{
    bip::triple<size_t> sizes2(32, 32, 1);

    try
    {

    bip::loggabor_filter_bank bof2a("loggabor2d", sizes2, 2, 4, 1, 0.33, 2.1, 0.55, 1.2, 15, 0.45, false);
    bip::loggabor_filter_bank::write_parameters(bof2a);

    std::cout << bof2a << std::endl;
    bof2a.compute();
    
    std::string filename2("loggabor2d_#_032_032_001_02_04_01_033_209_055_120_0.bof");
    bip::loggabor_filter_bank *bof2b = bip::loggabor_filter_bank::read_parameters(filename2);

    std::cout << *bof2b << std::endl;

    delete bof2b;

    // --------------------------------------------------------------------------------------------

    bip::triple<size_t> sizes3(32, 32, 32);
    
    bip::loggabor_filter_bank bof3a("loggabor3d", sizes3, 2, 4, 3, 0.33, 2.1, 0.55, 1.2, 15, 0.45, false);
    bip::loggabor_filter_bank::write_parameters(bof3a);

    std::cout << bof3a << std::endl;
    bof3a.compute();
    
    std::string filename3("loggabor3d_#_032_032_032_02_04_03_033_209_055_120_0.bof");
    bip::loggabor_filter_bank *bof3b = bip::loggabor_filter_bank::read_parameters(filename3);

    std::cout << *bof3b << std::endl;

    delete bof3b;

    }
    catch(const char *e)
    {
        std::cerr << e << std::endl;
    }


    return EXIT_SUCCESS;
}
