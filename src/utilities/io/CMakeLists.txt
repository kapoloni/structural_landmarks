CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT("io")

# -------------------------------------------------------------------------------------------------

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

# -------------------------------------------------------------------------------------------------

SET(LIBRARY_NAME
    "bip-io"
)
SET(LIBRARY_SOURCE_FILES
    "imageio.hpp"
    "meshio.hpp"
)
ADD_LIBRARY(
    ${LIBRARY_NAME}
    ${LIBRARY_SOURCE_FILES}
)
TARGET_LINK_LIBRARIES(
    ${LIBRARY_NAME}
    ${ITK_LIBRARIES}
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
    "test-io"
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
    ${ITK_LIBRARIES}
    bip-debug
    bip-types
    ${LIBRARY_NAME}
)
