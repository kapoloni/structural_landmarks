CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

PROJECT("mscproject")

SET(INCLUDE_OUTPUT_PATH
    ${CMAKE_BINARY_DIR}/include
    CACHE PATH "Single output directory for all includes."
)

SET(EXECUTABLE_OUTPUT_PATH
    ${CMAKE_BINARY_DIR}/bin
    CACHE PATH "Single output directory for building all executables."
)

SET(LIBRARY_OUTPUT_PATH
    ${CMAKE_BINARY_DIR}/lib
    CACHE PATH "Single output directory for building all libraries."
)

INCLUDE_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty
    ${CMAKE_CURRENT_SOURCE_DIR}/utilities/debug
    ${CMAKE_CURRENT_SOURCE_DIR}/utilities/types
    ${CMAKE_CURRENT_SOURCE_DIR}/utilities/math
    ${CMAKE_CURRENT_SOURCE_DIR}/utilities/io
    ${CMAKE_CURRENT_SOURCE_DIR}/filterbank
    ${CMAKE_CURRENT_SOURCE_DIR}/phasecongruency
    ${CMAKE_CURRENT_SOURCE_DIR}/landmarks/io
    ${CMAKE_CURRENT_SOURCE_DIR}/landmarks/detection
    ${CMAKE_CURRENT_SOURCE_DIR}/landmarks/matching
    ${CMAKE_CURRENT_SOURCE_DIR}/landmarks/registration
    ${CMAKE_CURRENT_SOURCE_DIR}/landmarks/atlas
    ${CMAKE_CURRENT_SOURCE_DIR}/utils
)

ADD_SUBDIRECTORY(utilities)
ADD_SUBDIRECTORY(landmarks)
ADD_SUBDIRECTORY(filterbank)
ADD_SUBDIRECTORY(phasecongruency)
