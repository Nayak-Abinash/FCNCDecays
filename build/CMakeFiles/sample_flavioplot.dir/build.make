# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/akn/Documents/GitHub/FCNCDecays

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/akn/Documents/GitHub/FCNCDecays/build

# Include any dependencies generated for this target.
include CMakeFiles/sample_flavioplot.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/sample_flavioplot.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/sample_flavioplot.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/sample_flavioplot.dir/flags.make

CMakeFiles/sample_flavioplot.dir/main/sample_flavioplot.cc.o: CMakeFiles/sample_flavioplot.dir/flags.make
CMakeFiles/sample_flavioplot.dir/main/sample_flavioplot.cc.o: ../main/sample_flavioplot.cc
CMakeFiles/sample_flavioplot.dir/main/sample_flavioplot.cc.o: CMakeFiles/sample_flavioplot.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/akn/Documents/GitHub/FCNCDecays/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/sample_flavioplot.dir/main/sample_flavioplot.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/sample_flavioplot.dir/main/sample_flavioplot.cc.o -MF CMakeFiles/sample_flavioplot.dir/main/sample_flavioplot.cc.o.d -o CMakeFiles/sample_flavioplot.dir/main/sample_flavioplot.cc.o -c /home/akn/Documents/GitHub/FCNCDecays/main/sample_flavioplot.cc

CMakeFiles/sample_flavioplot.dir/main/sample_flavioplot.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/sample_flavioplot.dir/main/sample_flavioplot.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/akn/Documents/GitHub/FCNCDecays/main/sample_flavioplot.cc > CMakeFiles/sample_flavioplot.dir/main/sample_flavioplot.cc.i

CMakeFiles/sample_flavioplot.dir/main/sample_flavioplot.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/sample_flavioplot.dir/main/sample_flavioplot.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/akn/Documents/GitHub/FCNCDecays/main/sample_flavioplot.cc -o CMakeFiles/sample_flavioplot.dir/main/sample_flavioplot.cc.s

# Object files for target sample_flavioplot
sample_flavioplot_OBJECTS = \
"CMakeFiles/sample_flavioplot.dir/main/sample_flavioplot.cc.o"

# External object files for target sample_flavioplot
sample_flavioplot_EXTERNAL_OBJECTS =

sample_flavioplot: CMakeFiles/sample_flavioplot.dir/main/sample_flavioplot.cc.o
sample_flavioplot: CMakeFiles/sample_flavioplot.dir/build.make
sample_flavioplot: libFCNC_shared.a
sample_flavioplot: /home/akn/root/lib/libCore.so
sample_flavioplot: /home/akn/root/lib/libImt.so
sample_flavioplot: /home/akn/root/lib/libRIO.so
sample_flavioplot: /home/akn/root/lib/libNet.so
sample_flavioplot: /home/akn/root/lib/libHist.so
sample_flavioplot: /home/akn/root/lib/libGraf.so
sample_flavioplot: /home/akn/root/lib/libGraf3d.so
sample_flavioplot: /home/akn/root/lib/libGpad.so
sample_flavioplot: /home/akn/root/lib/libROOTDataFrame.so
sample_flavioplot: /home/akn/root/lib/libTree.so
sample_flavioplot: /home/akn/root/lib/libTreePlayer.so
sample_flavioplot: /home/akn/root/lib/libRint.so
sample_flavioplot: /home/akn/root/lib/libPostscript.so
sample_flavioplot: /home/akn/root/lib/libMatrix.so
sample_flavioplot: /home/akn/root/lib/libPhysics.so
sample_flavioplot: /home/akn/root/lib/libMathCore.so
sample_flavioplot: /home/akn/root/lib/libThread.so
sample_flavioplot: /home/akn/root/lib/libMultiProc.so
sample_flavioplot: /home/akn/root/lib/libROOTVecOps.so
sample_flavioplot: CMakeFiles/sample_flavioplot.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/akn/Documents/GitHub/FCNCDecays/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable sample_flavioplot"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sample_flavioplot.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/sample_flavioplot.dir/build: sample_flavioplot
.PHONY : CMakeFiles/sample_flavioplot.dir/build

CMakeFiles/sample_flavioplot.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/sample_flavioplot.dir/cmake_clean.cmake
.PHONY : CMakeFiles/sample_flavioplot.dir/clean

CMakeFiles/sample_flavioplot.dir/depend:
	cd /home/akn/Documents/GitHub/FCNCDecays/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/akn/Documents/GitHub/FCNCDecays /home/akn/Documents/GitHub/FCNCDecays /home/akn/Documents/GitHub/FCNCDecays/build /home/akn/Documents/GitHub/FCNCDecays/build /home/akn/Documents/GitHub/FCNCDecays/build/CMakeFiles/sample_flavioplot.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/sample_flavioplot.dir/depend

