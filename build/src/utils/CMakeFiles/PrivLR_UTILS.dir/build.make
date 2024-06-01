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
CMAKE_SOURCE_DIR = /root/PrivLR

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /root/PrivLR/build

# Include any dependencies generated for this target.
include src/utils/CMakeFiles/PrivLR_UTILS.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/utils/CMakeFiles/PrivLR_UTILS.dir/compiler_depend.make

# Include the progress variables for this target.
include src/utils/CMakeFiles/PrivLR_UTILS.dir/progress.make

# Include the compile flags for this target's objects.
include src/utils/CMakeFiles/PrivLR_UTILS.dir/flags.make

src/utils/CMakeFiles/PrivLR_UTILS.dir/data.cpp.o: src/utils/CMakeFiles/PrivLR_UTILS.dir/flags.make
src/utils/CMakeFiles/PrivLR_UTILS.dir/data.cpp.o: ../src/utils/data.cpp
src/utils/CMakeFiles/PrivLR_UTILS.dir/data.cpp.o: src/utils/CMakeFiles/PrivLR_UTILS.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/root/PrivLR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/utils/CMakeFiles/PrivLR_UTILS.dir/data.cpp.o"
	cd /root/PrivLR/build/src/utils && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/utils/CMakeFiles/PrivLR_UTILS.dir/data.cpp.o -MF CMakeFiles/PrivLR_UTILS.dir/data.cpp.o.d -o CMakeFiles/PrivLR_UTILS.dir/data.cpp.o -c /root/PrivLR/src/utils/data.cpp

src/utils/CMakeFiles/PrivLR_UTILS.dir/data.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PrivLR_UTILS.dir/data.cpp.i"
	cd /root/PrivLR/build/src/utils && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /root/PrivLR/src/utils/data.cpp > CMakeFiles/PrivLR_UTILS.dir/data.cpp.i

src/utils/CMakeFiles/PrivLR_UTILS.dir/data.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PrivLR_UTILS.dir/data.cpp.s"
	cd /root/PrivLR/build/src/utils && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /root/PrivLR/src/utils/data.cpp -o CMakeFiles/PrivLR_UTILS.dir/data.cpp.s

src/utils/CMakeFiles/PrivLR_UTILS.dir/io.cpp.o: src/utils/CMakeFiles/PrivLR_UTILS.dir/flags.make
src/utils/CMakeFiles/PrivLR_UTILS.dir/io.cpp.o: ../src/utils/io.cpp
src/utils/CMakeFiles/PrivLR_UTILS.dir/io.cpp.o: src/utils/CMakeFiles/PrivLR_UTILS.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/root/PrivLR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/utils/CMakeFiles/PrivLR_UTILS.dir/io.cpp.o"
	cd /root/PrivLR/build/src/utils && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/utils/CMakeFiles/PrivLR_UTILS.dir/io.cpp.o -MF CMakeFiles/PrivLR_UTILS.dir/io.cpp.o.d -o CMakeFiles/PrivLR_UTILS.dir/io.cpp.o -c /root/PrivLR/src/utils/io.cpp

src/utils/CMakeFiles/PrivLR_UTILS.dir/io.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PrivLR_UTILS.dir/io.cpp.i"
	cd /root/PrivLR/build/src/utils && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /root/PrivLR/src/utils/io.cpp > CMakeFiles/PrivLR_UTILS.dir/io.cpp.i

src/utils/CMakeFiles/PrivLR_UTILS.dir/io.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PrivLR_UTILS.dir/io.cpp.s"
	cd /root/PrivLR/build/src/utils && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /root/PrivLR/src/utils/io.cpp -o CMakeFiles/PrivLR_UTILS.dir/io.cpp.s

src/utils/CMakeFiles/PrivLR_UTILS.dir/paillier.cpp.o: src/utils/CMakeFiles/PrivLR_UTILS.dir/flags.make
src/utils/CMakeFiles/PrivLR_UTILS.dir/paillier.cpp.o: ../src/utils/paillier.cpp
src/utils/CMakeFiles/PrivLR_UTILS.dir/paillier.cpp.o: src/utils/CMakeFiles/PrivLR_UTILS.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/root/PrivLR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/utils/CMakeFiles/PrivLR_UTILS.dir/paillier.cpp.o"
	cd /root/PrivLR/build/src/utils && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/utils/CMakeFiles/PrivLR_UTILS.dir/paillier.cpp.o -MF CMakeFiles/PrivLR_UTILS.dir/paillier.cpp.o.d -o CMakeFiles/PrivLR_UTILS.dir/paillier.cpp.o -c /root/PrivLR/src/utils/paillier.cpp

src/utils/CMakeFiles/PrivLR_UTILS.dir/paillier.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PrivLR_UTILS.dir/paillier.cpp.i"
	cd /root/PrivLR/build/src/utils && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /root/PrivLR/src/utils/paillier.cpp > CMakeFiles/PrivLR_UTILS.dir/paillier.cpp.i

src/utils/CMakeFiles/PrivLR_UTILS.dir/paillier.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PrivLR_UTILS.dir/paillier.cpp.s"
	cd /root/PrivLR/build/src/utils && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /root/PrivLR/src/utils/paillier.cpp -o CMakeFiles/PrivLR_UTILS.dir/paillier.cpp.s

# Object files for target PrivLR_UTILS
PrivLR_UTILS_OBJECTS = \
"CMakeFiles/PrivLR_UTILS.dir/data.cpp.o" \
"CMakeFiles/PrivLR_UTILS.dir/io.cpp.o" \
"CMakeFiles/PrivLR_UTILS.dir/paillier.cpp.o"

# External object files for target PrivLR_UTILS
PrivLR_UTILS_EXTERNAL_OBJECTS =

src/utils/libPrivLR_UTILS.a: src/utils/CMakeFiles/PrivLR_UTILS.dir/data.cpp.o
src/utils/libPrivLR_UTILS.a: src/utils/CMakeFiles/PrivLR_UTILS.dir/io.cpp.o
src/utils/libPrivLR_UTILS.a: src/utils/CMakeFiles/PrivLR_UTILS.dir/paillier.cpp.o
src/utils/libPrivLR_UTILS.a: src/utils/CMakeFiles/PrivLR_UTILS.dir/build.make
src/utils/libPrivLR_UTILS.a: src/utils/CMakeFiles/PrivLR_UTILS.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/root/PrivLR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX static library libPrivLR_UTILS.a"
	cd /root/PrivLR/build/src/utils && $(CMAKE_COMMAND) -P CMakeFiles/PrivLR_UTILS.dir/cmake_clean_target.cmake
	cd /root/PrivLR/build/src/utils && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/PrivLR_UTILS.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/utils/CMakeFiles/PrivLR_UTILS.dir/build: src/utils/libPrivLR_UTILS.a
.PHONY : src/utils/CMakeFiles/PrivLR_UTILS.dir/build

src/utils/CMakeFiles/PrivLR_UTILS.dir/clean:
	cd /root/PrivLR/build/src/utils && $(CMAKE_COMMAND) -P CMakeFiles/PrivLR_UTILS.dir/cmake_clean.cmake
.PHONY : src/utils/CMakeFiles/PrivLR_UTILS.dir/clean

src/utils/CMakeFiles/PrivLR_UTILS.dir/depend:
	cd /root/PrivLR/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /root/PrivLR /root/PrivLR/src/utils /root/PrivLR/build /root/PrivLR/build/src/utils /root/PrivLR/build/src/utils/CMakeFiles/PrivLR_UTILS.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/utils/CMakeFiles/PrivLR_UTILS.dir/depend
