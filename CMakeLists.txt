cmake_minimum_required(VERSION 3.5)
# set project name
project(FCNCDecays VERSION 1.0)

# add external libraries
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}/build)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS}/etc/cmake)
set(CMAKE_MACOSX_RPATH TRUE) # use @rpath keyword for the macOS rpath  

# add ROOT:
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)
include(${ROOT_USE_FILE})
add_definitions(${ROOT_CXX_FLAGS})

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

# make file system
file(GLOB source_files CONFIGURE_DEPENDS "src/*.cc")

# create library 
add_library(FCNC_shared ${source_files})

# create executables
add_executable(sample_q2plot main/sample_q2plot.cc)
add_executable(sample_mcerror main/sample_mcerror.cc)
add_executable(sample_observable main/sample_observable.cc)
add_executable(sample_covariance main/sample_covariance.cc)
add_executable(sample_flavioobs main/sample_flavioobs.cc)
add_executable(sample_flavioerr main/sample_flavioerr.cc)
add_executable(sample_flavioplot main/sample_flavioplot.cc)

# add custom libraries
target_link_libraries(sample_q2plot FCNC_shared ${ROOT_LIBRARIES})
target_link_libraries(sample_mcerror FCNC_shared ${ROOT_LIBRARIES})
target_link_libraries(sample_observable FCNC_shared ${ROOT_LIBRARIES})
target_link_libraries(sample_covariance FCNC_shared ${ROOT_LIBRARIES})
target_link_libraries(sample_flavioobs FCNC_shared ${ROOT_LIBRARIES})
target_link_libraries(sample_flavioerr FCNC_shared ${ROOT_LIBRARIES})
target_link_libraries(sample_flavioplot FCNC_shared ${ROOT_LIBRARIES})
