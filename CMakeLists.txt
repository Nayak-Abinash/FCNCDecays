cmake_minimum_required(VERSION 3.5)

# set project name
project(FCNCDecays VERSION 1.0)

# external libraries
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}/build)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS}/etc/cmake)
set(CMAKE_MACOSX_RPATH TRUE) # use @rpath keyword for the macOS rpath  

# ROOT:
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)
include(${ROOT_USE_FILE})
include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})




include_directories(include)

file(GLOB source_files CONFIGURE_DEPENDS "main/sample_q2plot.cc" "src/*.cc")

# add executable
add_executable(output ${source_files})

target_link_libraries(output ${ROOT_LIBRARIES} )
