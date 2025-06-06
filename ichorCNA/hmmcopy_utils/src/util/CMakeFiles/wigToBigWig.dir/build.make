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
CMAKE_SOURCE_DIR = /home/zjp/projects/202312_cfDNA_integrated_tools/ichorCNA/hmmcopy_utils

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/zjp/projects/202312_cfDNA_integrated_tools/ichorCNA/hmmcopy_utils

# Include any dependencies generated for this target.
include src/util/CMakeFiles/wigToBigWig.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/util/CMakeFiles/wigToBigWig.dir/compiler_depend.make

# Include the progress variables for this target.
include src/util/CMakeFiles/wigToBigWig.dir/progress.make

# Include the compile flags for this target's objects.
include src/util/CMakeFiles/wigToBigWig.dir/flags.make

src/util/CMakeFiles/wigToBigWig.dir/bigwig/wigToBigWig.c.o: src/util/CMakeFiles/wigToBigWig.dir/flags.make
src/util/CMakeFiles/wigToBigWig.dir/bigwig/wigToBigWig.c.o: src/util/bigwig/wigToBigWig.c
src/util/CMakeFiles/wigToBigWig.dir/bigwig/wigToBigWig.c.o: src/util/CMakeFiles/wigToBigWig.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zjp/projects/202312_cfDNA_integrated_tools/ichorCNA/hmmcopy_utils/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object src/util/CMakeFiles/wigToBigWig.dir/bigwig/wigToBigWig.c.o"
	cd /home/zjp/projects/202312_cfDNA_integrated_tools/ichorCNA/hmmcopy_utils/src/util && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT src/util/CMakeFiles/wigToBigWig.dir/bigwig/wigToBigWig.c.o -MF CMakeFiles/wigToBigWig.dir/bigwig/wigToBigWig.c.o.d -o CMakeFiles/wigToBigWig.dir/bigwig/wigToBigWig.c.o -c /home/zjp/projects/202312_cfDNA_integrated_tools/ichorCNA/hmmcopy_utils/src/util/bigwig/wigToBigWig.c

src/util/CMakeFiles/wigToBigWig.dir/bigwig/wigToBigWig.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/wigToBigWig.dir/bigwig/wigToBigWig.c.i"
	cd /home/zjp/projects/202312_cfDNA_integrated_tools/ichorCNA/hmmcopy_utils/src/util && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/zjp/projects/202312_cfDNA_integrated_tools/ichorCNA/hmmcopy_utils/src/util/bigwig/wigToBigWig.c > CMakeFiles/wigToBigWig.dir/bigwig/wigToBigWig.c.i

src/util/CMakeFiles/wigToBigWig.dir/bigwig/wigToBigWig.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/wigToBigWig.dir/bigwig/wigToBigWig.c.s"
	cd /home/zjp/projects/202312_cfDNA_integrated_tools/ichorCNA/hmmcopy_utils/src/util && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/zjp/projects/202312_cfDNA_integrated_tools/ichorCNA/hmmcopy_utils/src/util/bigwig/wigToBigWig.c -o CMakeFiles/wigToBigWig.dir/bigwig/wigToBigWig.c.s

# Object files for target wigToBigWig
wigToBigWig_OBJECTS = \
"CMakeFiles/wigToBigWig.dir/bigwig/wigToBigWig.c.o"

# External object files for target wigToBigWig
wigToBigWig_EXTERNAL_OBJECTS =

util/bigwig/wigToBigWig: src/util/CMakeFiles/wigToBigWig.dir/bigwig/wigToBigWig.c.o
util/bigwig/wigToBigWig: src/util/CMakeFiles/wigToBigWig.dir/build.make
util/bigwig/wigToBigWig: lib/libkent.a
util/bigwig/wigToBigWig: src/util/CMakeFiles/wigToBigWig.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/zjp/projects/202312_cfDNA_integrated_tools/ichorCNA/hmmcopy_utils/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable ../../util/bigwig/wigToBigWig"
	cd /home/zjp/projects/202312_cfDNA_integrated_tools/ichorCNA/hmmcopy_utils/src/util && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/wigToBigWig.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/util/CMakeFiles/wigToBigWig.dir/build: util/bigwig/wigToBigWig
.PHONY : src/util/CMakeFiles/wigToBigWig.dir/build

src/util/CMakeFiles/wigToBigWig.dir/clean:
	cd /home/zjp/projects/202312_cfDNA_integrated_tools/ichorCNA/hmmcopy_utils/src/util && $(CMAKE_COMMAND) -P CMakeFiles/wigToBigWig.dir/cmake_clean.cmake
.PHONY : src/util/CMakeFiles/wigToBigWig.dir/clean

src/util/CMakeFiles/wigToBigWig.dir/depend:
	cd /home/zjp/projects/202312_cfDNA_integrated_tools/ichorCNA/hmmcopy_utils && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/zjp/projects/202312_cfDNA_integrated_tools/ichorCNA/hmmcopy_utils /home/zjp/projects/202312_cfDNA_integrated_tools/ichorCNA/hmmcopy_utils/src/util /home/zjp/projects/202312_cfDNA_integrated_tools/ichorCNA/hmmcopy_utils /home/zjp/projects/202312_cfDNA_integrated_tools/ichorCNA/hmmcopy_utils/src/util /home/zjp/projects/202312_cfDNA_integrated_tools/ichorCNA/hmmcopy_utils/src/util/CMakeFiles/wigToBigWig.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/util/CMakeFiles/wigToBigWig.dir/depend

