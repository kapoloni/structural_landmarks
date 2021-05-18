CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

PROJECT("filterbank")

# -------------------------------------------------------------------------------------------------

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

# -------------------------------------------------------------------------------------------------

SET(LIBRARY_NAME
    "bip-filterbank"
)
SET(LIBRARY_SOURCE_FILES
    "loggaborfilterbank.cpp"
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
)
TARGET_INCLUDE_DIRECTORIES(
    ${LIBRARY_NAME}
    PUBLIC
    ${CMAKE_SOURCE_DIR}/${PROJECT_NAME}
)

# -------------------------------------------------------------------------------------------------

SET(PROGRAM_NAME
    "test-filterbank"
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
    bip-math
    bip-io
    ${LIBRARY_NAME}
)


# -------------------------------------------------------------------------------------------------

SET(PROGRAM_NAME
    "log-gabor-filter-bank"
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
    ${LIBRARY_NAME}
)
