CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT("debug")

# -------------------------------------------------------------------------------------------------

SET(LIBRARY_NAME
    "bip-debug"
)
SET(LIBRARY_SOURCE_FILES
    "assert.cpp"
)
ADD_LIBRARY(
    ${LIBRARY_NAME}
    ${LIBRARY_SOURCE_FILES}
)
TARGET_LINK_LIBRARIES(
    ${LIBRARY_NAME}
)
TARGET_INCLUDE_DIRECTORIES(
    ${LIBRARY_NAME}
    PUBLIC
    ${CMAKE_SOURCE_DIR}/${PROJECT_NAME}
)

# -------------------------------------------------------------------------------------------------

SET(PROGRAM_NAME
    "test-debug"
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
    ${LIBRARY_NAME}
)
