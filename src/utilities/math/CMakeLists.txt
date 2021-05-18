CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT("math")

# -------------------------------------------------------------------------------------------------

SET(LIBRARY_NAME
    "bip-math"
)
SET(LIBRARY_SOURCE_FILES
    "mathfunctions.cpp"
)
ADD_LIBRARY(
    ${LIBRARY_NAME}
    ${LIBRARY_SOURCE_FILES}
)
TARGET_LINK_LIBRARIES(
    ${LIBRARY_NAME}
    m
    fftw3f
    fftw3f_threads
    gsl
    gslcblas
    bip-debug
    bip-types
)
TARGET_INCLUDE_DIRECTORIES(
    ${LIBRARY_NAME}
    PUBLIC
    ${CMAKE_SOURCE_DIR}/${PROJECT_NAME}
)

# -------------------------------------------------------------------------------------------------

SET(PROGRAM_NAME
    "test-math"
)
SET(PROGRAM_SOURCE_FILES
    "__test.cpp"
)
ADD_EXECUTABLE(
    ${PROGRAM_NAME}
    ${PROGRAM_SOURCE_FILES}
)
TARGET_LINK_LIBRARIES(
    ${PROGRAM_NAME}
    m
    fftw3f
    fftw3f_threads
    gsl
    gslcblas
    bip-debug
    bip-types
    ${LIBRARY_NAME}
)
