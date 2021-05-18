CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT("landmark-detection")

# -------------------------------------------------------------------------------------------------

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

# -------------------------------------------------------------------------------------------------

SET(LIBRARY_NAME
    "bip-landmark-detection"
)
SET(LIBRARY_SOURCE_FILES
    "landmarkdetector.cpp"
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
    bip-math
    bip-io
    bip-landmark-io
)
TARGET_INCLUDE_DIRECTORIES(
    ${LIBRARY_NAME}
    PUBLIC
    ${CMAKE_SOURCE_DIR}/${PROJECT_NAME}
)


# -------------------------------------------------------------------------------------------------

SET(PROGRAM_NAME
    "landmark-detector"
)
SET(PROGRAM_SOURCE_FILES
  "__main.cpp"
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
    bip-math
    bip-io
    bip-landmark-io
    ${LIBRARY_NAME}
)
