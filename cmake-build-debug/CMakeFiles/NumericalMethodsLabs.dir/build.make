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
CMAKE_COMMAND = /snap/clion/126/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /snap/clion/126/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ivan/CLionProjects/NumericalMethodsLabs

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ivan/CLionProjects/NumericalMethodsLabs/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/NumericalMethodsLabs.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/NumericalMethodsLabs.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/NumericalMethodsLabs.dir/flags.make

CMakeFiles/NumericalMethodsLabs.dir/main.cpp.o: CMakeFiles/NumericalMethodsLabs.dir/flags.make
CMakeFiles/NumericalMethodsLabs.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ivan/CLionProjects/NumericalMethodsLabs/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/NumericalMethodsLabs.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/NumericalMethodsLabs.dir/main.cpp.o -c /home/ivan/CLionProjects/NumericalMethodsLabs/main.cpp

CMakeFiles/NumericalMethodsLabs.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NumericalMethodsLabs.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ivan/CLionProjects/NumericalMethodsLabs/main.cpp > CMakeFiles/NumericalMethodsLabs.dir/main.cpp.i

CMakeFiles/NumericalMethodsLabs.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NumericalMethodsLabs.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ivan/CLionProjects/NumericalMethodsLabs/main.cpp -o CMakeFiles/NumericalMethodsLabs.dir/main.cpp.s

CMakeFiles/NumericalMethodsLabs.dir/SLAE.cpp.o: CMakeFiles/NumericalMethodsLabs.dir/flags.make
CMakeFiles/NumericalMethodsLabs.dir/SLAE.cpp.o: ../SLAE.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ivan/CLionProjects/NumericalMethodsLabs/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/NumericalMethodsLabs.dir/SLAE.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/NumericalMethodsLabs.dir/SLAE.cpp.o -c /home/ivan/CLionProjects/NumericalMethodsLabs/SLAE.cpp

CMakeFiles/NumericalMethodsLabs.dir/SLAE.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NumericalMethodsLabs.dir/SLAE.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ivan/CLionProjects/NumericalMethodsLabs/SLAE.cpp > CMakeFiles/NumericalMethodsLabs.dir/SLAE.cpp.i

CMakeFiles/NumericalMethodsLabs.dir/SLAE.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NumericalMethodsLabs.dir/SLAE.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ivan/CLionProjects/NumericalMethodsLabs/SLAE.cpp -o CMakeFiles/NumericalMethodsLabs.dir/SLAE.cpp.s

# Object files for target NumericalMethodsLabs
NumericalMethodsLabs_OBJECTS = \
"CMakeFiles/NumericalMethodsLabs.dir/main.cpp.o" \
"CMakeFiles/NumericalMethodsLabs.dir/SLAE.cpp.o"

# External object files for target NumericalMethodsLabs
NumericalMethodsLabs_EXTERNAL_OBJECTS =

NumericalMethodsLabs: CMakeFiles/NumericalMethodsLabs.dir/main.cpp.o
NumericalMethodsLabs: CMakeFiles/NumericalMethodsLabs.dir/SLAE.cpp.o
NumericalMethodsLabs: CMakeFiles/NumericalMethodsLabs.dir/build.make
NumericalMethodsLabs: CMakeFiles/NumericalMethodsLabs.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ivan/CLionProjects/NumericalMethodsLabs/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable NumericalMethodsLabs"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/NumericalMethodsLabs.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/NumericalMethodsLabs.dir/build: NumericalMethodsLabs

.PHONY : CMakeFiles/NumericalMethodsLabs.dir/build

CMakeFiles/NumericalMethodsLabs.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/NumericalMethodsLabs.dir/cmake_clean.cmake
.PHONY : CMakeFiles/NumericalMethodsLabs.dir/clean

CMakeFiles/NumericalMethodsLabs.dir/depend:
	cd /home/ivan/CLionProjects/NumericalMethodsLabs/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ivan/CLionProjects/NumericalMethodsLabs /home/ivan/CLionProjects/NumericalMethodsLabs /home/ivan/CLionProjects/NumericalMethodsLabs/cmake-build-debug /home/ivan/CLionProjects/NumericalMethodsLabs/cmake-build-debug /home/ivan/CLionProjects/NumericalMethodsLabs/cmake-build-debug/CMakeFiles/NumericalMethodsLabs.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/NumericalMethodsLabs.dir/depend

