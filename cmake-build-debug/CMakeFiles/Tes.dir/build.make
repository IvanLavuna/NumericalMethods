# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

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

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /snap/clion/129/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /snap/clion/129/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ivan/CLionProjects/NumericalMethodsLabs

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ivan/CLionProjects/NumericalMethodsLabs/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/Tes.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Tes.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Tes.dir/flags.make

CMakeFiles/Tes.dir/testing/TestGenerator.cpp.o: CMakeFiles/Tes.dir/flags.make
CMakeFiles/Tes.dir/testing/TestGenerator.cpp.o: ../testing/TestGenerator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ivan/CLionProjects/NumericalMethodsLabs/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Tes.dir/testing/TestGenerator.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Tes.dir/testing/TestGenerator.cpp.o -c /home/ivan/CLionProjects/NumericalMethodsLabs/testing/TestGenerator.cpp

CMakeFiles/Tes.dir/testing/TestGenerator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Tes.dir/testing/TestGenerator.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ivan/CLionProjects/NumericalMethodsLabs/testing/TestGenerator.cpp > CMakeFiles/Tes.dir/testing/TestGenerator.cpp.i

CMakeFiles/Tes.dir/testing/TestGenerator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Tes.dir/testing/TestGenerator.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ivan/CLionProjects/NumericalMethodsLabs/testing/TestGenerator.cpp -o CMakeFiles/Tes.dir/testing/TestGenerator.cpp.s

CMakeFiles/Tes.dir/testing/Gauss_Jordan.cpp.o: CMakeFiles/Tes.dir/flags.make
CMakeFiles/Tes.dir/testing/Gauss_Jordan.cpp.o: ../testing/Gauss_Jordan.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ivan/CLionProjects/NumericalMethodsLabs/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/Tes.dir/testing/Gauss_Jordan.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Tes.dir/testing/Gauss_Jordan.cpp.o -c /home/ivan/CLionProjects/NumericalMethodsLabs/testing/Gauss_Jordan.cpp

CMakeFiles/Tes.dir/testing/Gauss_Jordan.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Tes.dir/testing/Gauss_Jordan.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ivan/CLionProjects/NumericalMethodsLabs/testing/Gauss_Jordan.cpp > CMakeFiles/Tes.dir/testing/Gauss_Jordan.cpp.i

CMakeFiles/Tes.dir/testing/Gauss_Jordan.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Tes.dir/testing/Gauss_Jordan.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ivan/CLionProjects/NumericalMethodsLabs/testing/Gauss_Jordan.cpp -o CMakeFiles/Tes.dir/testing/Gauss_Jordan.cpp.s

# Object files for target Tes
Tes_OBJECTS = \
"CMakeFiles/Tes.dir/testing/TestGenerator.cpp.o" \
"CMakeFiles/Tes.dir/testing/Gauss_Jordan.cpp.o"

# External object files for target Tes
Tes_EXTERNAL_OBJECTS =

Tes: CMakeFiles/Tes.dir/testing/TestGenerator.cpp.o
Tes: CMakeFiles/Tes.dir/testing/Gauss_Jordan.cpp.o
Tes: CMakeFiles/Tes.dir/build.make
Tes: CMakeFiles/Tes.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ivan/CLionProjects/NumericalMethodsLabs/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable Tes"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Tes.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Tes.dir/build: Tes

.PHONY : CMakeFiles/Tes.dir/build

CMakeFiles/Tes.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Tes.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Tes.dir/clean

CMakeFiles/Tes.dir/depend:
	cd /home/ivan/CLionProjects/NumericalMethodsLabs/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ivan/CLionProjects/NumericalMethodsLabs /home/ivan/CLionProjects/NumericalMethodsLabs /home/ivan/CLionProjects/NumericalMethodsLabs/cmake-build-debug /home/ivan/CLionProjects/NumericalMethodsLabs/cmake-build-debug /home/ivan/CLionProjects/NumericalMethodsLabs/cmake-build-debug/CMakeFiles/Tes.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Tes.dir/depend
