CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(osfi-1)

ADD_DEFINITIONS(-std=c++11 -Wall)

FILE(GLOB_RECURSE LIBRARY_SOURCES_CPP src/lib/*.cpp)
FILE(GLOB_RECURSE LIBRARY_SOURCES_C src/lib/*.c)

ADD_LIBRARY(lib_cpp STATIC ${LIBRARY_SOURCES_CPP})

ADD_EXECUTABLE(test src/test.cpp)
TARGET_LINK_LIBRARIES(test lib_cpp)

ADD_EXECUTABLE(gen-hdr src/gen-hdr.cpp)
TARGET_LINK_LIBRARIES(gen-hdr lib_cpp)

ADD_EXECUTABLE(convert-to-ldr src/convert-to-ldr.cpp)
TARGET_LINK_LIBRARIES(convert-to-ldr lib_cpp)
