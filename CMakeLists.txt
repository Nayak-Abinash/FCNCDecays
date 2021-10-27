cmake_minimum_required(VERSION 3.7)

# set project name
project(FCNCDecays VERSION 1.0)

include_directories(include)

file(GLOB source_files CONFIGURE_DEPENDS "main/main.cc" "src/*.cc")

# add executable
add_executable(output ${source_files})
