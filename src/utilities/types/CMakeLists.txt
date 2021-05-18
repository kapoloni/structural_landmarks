CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT("types")

# -------------------------------------------------------------------------------------------------

SET(LIBRARY_NAME
    "bip-types"
)
SET(LIBRARY_SOURCE_FILES
    "triple.hpp"
    "landmark.cpp"
)
ADD_LIBRARY(
    ${LIBRARY_NAME}
    ${LIBRARY_SOURCE_FILES}
)
TARGET_LINK_LIBRARIES(
    ${LIBRARY_NAME}
    bip-debug
)
TARGET_INCLUDE_DIRECTORIES(
    ${LIBRARY_NAME}
    PUBLIC
    ${CMAKE_SOURCE_DIR}/${PROJECT_NAME}
)

# -------------------------------------------------------------------------------------------------

SET(PROGRAM_NAME
    "test-types"
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
    bip-debug
    ${LIBRARY_NAME}
)
