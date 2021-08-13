# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/enrg/Software/WMUtils/ratpac

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/enrg/Software/WMUtils/ratpac/build

# Include any dependencies generated for this target.
include tools/root/CMakeFiles/root.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include tools/root/CMakeFiles/root.dir/compiler_depend.make

# Include the progress variables for this target.
include tools/root/CMakeFiles/root.dir/progress.make

# Include the compile flags for this target's objects.
include tools/root/CMakeFiles/root.dir/flags.make

tools/root/CMakeFiles/root.dir/root.cc.o: tools/root/CMakeFiles/root.dir/flags.make
tools/root/CMakeFiles/root.dir/root.cc.o: ../tools/root/root.cc
tools/root/CMakeFiles/root.dir/root.cc.o: tools/root/CMakeFiles/root.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enrg/Software/WMUtils/ratpac/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tools/root/CMakeFiles/root.dir/root.cc.o"
	cd /home/enrg/Software/WMUtils/ratpac/build/tools/root && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT tools/root/CMakeFiles/root.dir/root.cc.o -MF CMakeFiles/root.dir/root.cc.o.d -o CMakeFiles/root.dir/root.cc.o -c /home/enrg/Software/WMUtils/ratpac/tools/root/root.cc

tools/root/CMakeFiles/root.dir/root.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/root.dir/root.cc.i"
	cd /home/enrg/Software/WMUtils/ratpac/build/tools/root && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/enrg/Software/WMUtils/ratpac/tools/root/root.cc > CMakeFiles/root.dir/root.cc.i

tools/root/CMakeFiles/root.dir/root.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/root.dir/root.cc.s"
	cd /home/enrg/Software/WMUtils/ratpac/build/tools/root && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/enrg/Software/WMUtils/ratpac/tools/root/root.cc -o CMakeFiles/root.dir/root.cc.s

# Object files for target root
root_OBJECTS = \
"CMakeFiles/root.dir/root.cc.o"

# External object files for target root
root_EXTERNAL_OBJECTS =

bin/root: tools/root/CMakeFiles/root.dir/root.cc.o
bin/root: tools/root/CMakeFiles/root.dir/build.make
bin/root: /home/enrg/Software/WMUtils/local/lib/libCore.so
bin/root: /home/enrg/Software/WMUtils/local/lib/libImt.so
bin/root: /home/enrg/Software/WMUtils/local/lib/libRIO.so
bin/root: /home/enrg/Software/WMUtils/local/lib/libNet.so
bin/root: /home/enrg/Software/WMUtils/local/lib/libHist.so
bin/root: /home/enrg/Software/WMUtils/local/lib/libGraf.so
bin/root: /home/enrg/Software/WMUtils/local/lib/libGraf3d.so
bin/root: /home/enrg/Software/WMUtils/local/lib/libGpad.so
bin/root: /home/enrg/Software/WMUtils/local/lib/libROOTDataFrame.so
bin/root: /home/enrg/Software/WMUtils/local/lib/libTree.so
bin/root: /home/enrg/Software/WMUtils/local/lib/libTreePlayer.so
bin/root: /home/enrg/Software/WMUtils/local/lib/libRint.so
bin/root: /home/enrg/Software/WMUtils/local/lib/libPostscript.so
bin/root: /home/enrg/Software/WMUtils/local/lib/libMatrix.so
bin/root: /home/enrg/Software/WMUtils/local/lib/libPhysics.so
bin/root: /home/enrg/Software/WMUtils/local/lib/libMathCore.so
bin/root: /home/enrg/Software/WMUtils/local/lib/libThread.so
bin/root: /home/enrg/Software/WMUtils/local/lib/libMultiProc.so
bin/root: /home/enrg/Software/WMUtils/local/lib/libMinuit2.so
bin/root: /home/enrg/Software/WMUtils/local/lib/libPyROOT.so
bin/root: tools/root/CMakeFiles/root.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/enrg/Software/WMUtils/ratpac/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../bin/root"
	cd /home/enrg/Software/WMUtils/ratpac/build/tools/root && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/root.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tools/root/CMakeFiles/root.dir/build: bin/root
.PHONY : tools/root/CMakeFiles/root.dir/build

tools/root/CMakeFiles/root.dir/clean:
	cd /home/enrg/Software/WMUtils/ratpac/build/tools/root && $(CMAKE_COMMAND) -P CMakeFiles/root.dir/cmake_clean.cmake
.PHONY : tools/root/CMakeFiles/root.dir/clean

tools/root/CMakeFiles/root.dir/depend:
	cd /home/enrg/Software/WMUtils/ratpac/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/enrg/Software/WMUtils/ratpac /home/enrg/Software/WMUtils/ratpac/tools/root /home/enrg/Software/WMUtils/ratpac/build /home/enrg/Software/WMUtils/ratpac/build/tools/root /home/enrg/Software/WMUtils/ratpac/build/tools/root/CMakeFiles/root.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tools/root/CMakeFiles/root.dir/depend

