cmake_minimum_required (VERSION 2.8.11)
project (FCNCDecays)
enable_testing()

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}/build)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS}/etc/cmake)
set(CMAKE_MACOSX_RPATH TRUE) # use @rpath keyword for the macOS rpath

###########################
#
#    Compile settings
#
###########################
# Add -fPIC for 64 bit systems
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif()
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()



###########################
#
#    Packages search
#
###########################
# looking for yaml-cpp
find_package(Yaml-cpp)
# looking for boost
find_package(Boost REQUIRED)
# GSL
find_package(GSL REQUIRED)
# ROOT:
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)

include(${ROOT_USE_FILE})
include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})


###########################
#
#    Config   Variables
#
###########################
if(NOT DEFINED EXECUTABLE)
  set(EXECUTABLE FALSE)
endif()

if(NOT DEFINED USE_ROOT)
  set(USE_ROOT TRUE)
endif()


###########################
#
#    Libraries
#
###########################
if(USE_ROOT)

  file(GLOB FD_SOURCES ${CMAKE_SOURCE_DIR}/src/*.cc)
else()
  file(GLOB FD_SOURCES ${CMAKE_SOURCE_DIR}/src/*.cc)

endif()  

add_library(FCNCDecays_static ${FD_SOURCES})
add_library(FCNCDecays_shared SHARED ${FD_SOURCES})
set_target_properties(FCNCDecays_static PROPERTIES OUTPUT_NAME FCNCDecays)
set_target_properties(FCNCDecays_shared PROPERTIES OUTPUT_NAME FCNCDecays)
target_link_libraries(FCNCDecays_shared PRIVATE ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${GSL_LIBRARIES} )


###########################
#
#    Executables
#
###########################

if(EXECUTABLE AND USE_ROOT)
  
  add_executable(output main/main.cc)


  target_link_libraries(output FCNCDecays_shared yaml-cpp ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )


  
endif()
