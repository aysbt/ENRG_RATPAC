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
include src/core/CMakeFiles/core.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/core/CMakeFiles/core.dir/compiler_depend.make

# Include the progress variables for this target.
include src/core/CMakeFiles/core.dir/progress.make

# Include the compile flags for this target's objects.
include src/core/CMakeFiles/core.dir/flags.make

src/core/CMakeFiles/core.dir/src/GLG4HitPhoton.cc.o: src/core/CMakeFiles/core.dir/flags.make
src/core/CMakeFiles/core.dir/src/GLG4HitPhoton.cc.o: ../src/core/src/GLG4HitPhoton.cc
src/core/CMakeFiles/core.dir/src/GLG4HitPhoton.cc.o: src/core/CMakeFiles/core.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enrg/Software/WMUtils/ratpac/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/core/CMakeFiles/core.dir/src/GLG4HitPhoton.cc.o"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/core/CMakeFiles/core.dir/src/GLG4HitPhoton.cc.o -MF CMakeFiles/core.dir/src/GLG4HitPhoton.cc.o.d -o CMakeFiles/core.dir/src/GLG4HitPhoton.cc.o -c /home/enrg/Software/WMUtils/ratpac/src/core/src/GLG4HitPhoton.cc

src/core/CMakeFiles/core.dir/src/GLG4HitPhoton.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/core.dir/src/GLG4HitPhoton.cc.i"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/enrg/Software/WMUtils/ratpac/src/core/src/GLG4HitPhoton.cc > CMakeFiles/core.dir/src/GLG4HitPhoton.cc.i

src/core/CMakeFiles/core.dir/src/GLG4HitPhoton.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/core.dir/src/GLG4HitPhoton.cc.s"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/enrg/Software/WMUtils/ratpac/src/core/src/GLG4HitPhoton.cc -o CMakeFiles/core.dir/src/GLG4HitPhoton.cc.s

src/core/CMakeFiles/core.dir/src/GLG4HitPMT.cc.o: src/core/CMakeFiles/core.dir/flags.make
src/core/CMakeFiles/core.dir/src/GLG4HitPMT.cc.o: ../src/core/src/GLG4HitPMT.cc
src/core/CMakeFiles/core.dir/src/GLG4HitPMT.cc.o: src/core/CMakeFiles/core.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enrg/Software/WMUtils/ratpac/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/core/CMakeFiles/core.dir/src/GLG4HitPMT.cc.o"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/core/CMakeFiles/core.dir/src/GLG4HitPMT.cc.o -MF CMakeFiles/core.dir/src/GLG4HitPMT.cc.o.d -o CMakeFiles/core.dir/src/GLG4HitPMT.cc.o -c /home/enrg/Software/WMUtils/ratpac/src/core/src/GLG4HitPMT.cc

src/core/CMakeFiles/core.dir/src/GLG4HitPMT.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/core.dir/src/GLG4HitPMT.cc.i"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/enrg/Software/WMUtils/ratpac/src/core/src/GLG4HitPMT.cc > CMakeFiles/core.dir/src/GLG4HitPMT.cc.i

src/core/CMakeFiles/core.dir/src/GLG4HitPMT.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/core.dir/src/GLG4HitPMT.cc.s"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/enrg/Software/WMUtils/ratpac/src/core/src/GLG4HitPMT.cc -o CMakeFiles/core.dir/src/GLG4HitPMT.cc.s

src/core/CMakeFiles/core.dir/src/GLG4SteppingAction.cc.o: src/core/CMakeFiles/core.dir/flags.make
src/core/CMakeFiles/core.dir/src/GLG4SteppingAction.cc.o: ../src/core/src/GLG4SteppingAction.cc
src/core/CMakeFiles/core.dir/src/GLG4SteppingAction.cc.o: src/core/CMakeFiles/core.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enrg/Software/WMUtils/ratpac/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/core/CMakeFiles/core.dir/src/GLG4SteppingAction.cc.o"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/core/CMakeFiles/core.dir/src/GLG4SteppingAction.cc.o -MF CMakeFiles/core.dir/src/GLG4SteppingAction.cc.o.d -o CMakeFiles/core.dir/src/GLG4SteppingAction.cc.o -c /home/enrg/Software/WMUtils/ratpac/src/core/src/GLG4SteppingAction.cc

src/core/CMakeFiles/core.dir/src/GLG4SteppingAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/core.dir/src/GLG4SteppingAction.cc.i"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/enrg/Software/WMUtils/ratpac/src/core/src/GLG4SteppingAction.cc > CMakeFiles/core.dir/src/GLG4SteppingAction.cc.i

src/core/CMakeFiles/core.dir/src/GLG4SteppingAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/core.dir/src/GLG4SteppingAction.cc.s"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/enrg/Software/WMUtils/ratpac/src/core/src/GLG4SteppingAction.cc -o CMakeFiles/core.dir/src/GLG4SteppingAction.cc.s

src/core/CMakeFiles/core.dir/src/PythonProc.cc.o: src/core/CMakeFiles/core.dir/flags.make
src/core/CMakeFiles/core.dir/src/PythonProc.cc.o: ../src/core/src/PythonProc.cc
src/core/CMakeFiles/core.dir/src/PythonProc.cc.o: src/core/CMakeFiles/core.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enrg/Software/WMUtils/ratpac/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object src/core/CMakeFiles/core.dir/src/PythonProc.cc.o"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/core/CMakeFiles/core.dir/src/PythonProc.cc.o -MF CMakeFiles/core.dir/src/PythonProc.cc.o.d -o CMakeFiles/core.dir/src/PythonProc.cc.o -c /home/enrg/Software/WMUtils/ratpac/src/core/src/PythonProc.cc

src/core/CMakeFiles/core.dir/src/PythonProc.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/core.dir/src/PythonProc.cc.i"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/enrg/Software/WMUtils/ratpac/src/core/src/PythonProc.cc > CMakeFiles/core.dir/src/PythonProc.cc.i

src/core/CMakeFiles/core.dir/src/PythonProc.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/core.dir/src/PythonProc.cc.s"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/enrg/Software/WMUtils/ratpac/src/core/src/PythonProc.cc -o CMakeFiles/core.dir/src/PythonProc.cc.s

src/core/CMakeFiles/core.dir/src/ProcBlock.cc.o: src/core/CMakeFiles/core.dir/flags.make
src/core/CMakeFiles/core.dir/src/ProcBlock.cc.o: ../src/core/src/ProcBlock.cc
src/core/CMakeFiles/core.dir/src/ProcBlock.cc.o: src/core/CMakeFiles/core.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enrg/Software/WMUtils/ratpac/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object src/core/CMakeFiles/core.dir/src/ProcBlock.cc.o"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/core/CMakeFiles/core.dir/src/ProcBlock.cc.o -MF CMakeFiles/core.dir/src/ProcBlock.cc.o.d -o CMakeFiles/core.dir/src/ProcBlock.cc.o -c /home/enrg/Software/WMUtils/ratpac/src/core/src/ProcBlock.cc

src/core/CMakeFiles/core.dir/src/ProcBlock.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/core.dir/src/ProcBlock.cc.i"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/enrg/Software/WMUtils/ratpac/src/core/src/ProcBlock.cc > CMakeFiles/core.dir/src/ProcBlock.cc.i

src/core/CMakeFiles/core.dir/src/ProcBlock.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/core.dir/src/ProcBlock.cc.s"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/enrg/Software/WMUtils/ratpac/src/core/src/ProcBlock.cc -o CMakeFiles/core.dir/src/ProcBlock.cc.s

src/core/CMakeFiles/core.dir/src/RunManager.cc.o: src/core/CMakeFiles/core.dir/flags.make
src/core/CMakeFiles/core.dir/src/RunManager.cc.o: ../src/core/src/RunManager.cc
src/core/CMakeFiles/core.dir/src/RunManager.cc.o: src/core/CMakeFiles/core.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enrg/Software/WMUtils/ratpac/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object src/core/CMakeFiles/core.dir/src/RunManager.cc.o"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/core/CMakeFiles/core.dir/src/RunManager.cc.o -MF CMakeFiles/core.dir/src/RunManager.cc.o.d -o CMakeFiles/core.dir/src/RunManager.cc.o -c /home/enrg/Software/WMUtils/ratpac/src/core/src/RunManager.cc

src/core/CMakeFiles/core.dir/src/RunManager.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/core.dir/src/RunManager.cc.i"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/enrg/Software/WMUtils/ratpac/src/core/src/RunManager.cc > CMakeFiles/core.dir/src/RunManager.cc.i

src/core/CMakeFiles/core.dir/src/RunManager.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/core.dir/src/RunManager.cc.s"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/enrg/Software/WMUtils/ratpac/src/core/src/RunManager.cc -o CMakeFiles/core.dir/src/RunManager.cc.s

src/core/CMakeFiles/core.dir/src/CountProc.cc.o: src/core/CMakeFiles/core.dir/flags.make
src/core/CMakeFiles/core.dir/src/CountProc.cc.o: ../src/core/src/CountProc.cc
src/core/CMakeFiles/core.dir/src/CountProc.cc.o: src/core/CMakeFiles/core.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enrg/Software/WMUtils/ratpac/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object src/core/CMakeFiles/core.dir/src/CountProc.cc.o"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/core/CMakeFiles/core.dir/src/CountProc.cc.o -MF CMakeFiles/core.dir/src/CountProc.cc.o.d -o CMakeFiles/core.dir/src/CountProc.cc.o -c /home/enrg/Software/WMUtils/ratpac/src/core/src/CountProc.cc

src/core/CMakeFiles/core.dir/src/CountProc.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/core.dir/src/CountProc.cc.i"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/enrg/Software/WMUtils/ratpac/src/core/src/CountProc.cc > CMakeFiles/core.dir/src/CountProc.cc.i

src/core/CMakeFiles/core.dir/src/CountProc.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/core.dir/src/CountProc.cc.s"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/enrg/Software/WMUtils/ratpac/src/core/src/CountProc.cc -o CMakeFiles/core.dir/src/CountProc.cc.s

src/core/CMakeFiles/core.dir/src/PruneProc.cc.o: src/core/CMakeFiles/core.dir/flags.make
src/core/CMakeFiles/core.dir/src/PruneProc.cc.o: ../src/core/src/PruneProc.cc
src/core/CMakeFiles/core.dir/src/PruneProc.cc.o: src/core/CMakeFiles/core.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enrg/Software/WMUtils/ratpac/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object src/core/CMakeFiles/core.dir/src/PruneProc.cc.o"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/core/CMakeFiles/core.dir/src/PruneProc.cc.o -MF CMakeFiles/core.dir/src/PruneProc.cc.o.d -o CMakeFiles/core.dir/src/PruneProc.cc.o -c /home/enrg/Software/WMUtils/ratpac/src/core/src/PruneProc.cc

src/core/CMakeFiles/core.dir/src/PruneProc.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/core.dir/src/PruneProc.cc.i"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/enrg/Software/WMUtils/ratpac/src/core/src/PruneProc.cc > CMakeFiles/core.dir/src/PruneProc.cc.i

src/core/CMakeFiles/core.dir/src/PruneProc.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/core.dir/src/PruneProc.cc.s"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/enrg/Software/WMUtils/ratpac/src/core/src/PruneProc.cc -o CMakeFiles/core.dir/src/PruneProc.cc.s

src/core/CMakeFiles/core.dir/src/Gsim.cc.o: src/core/CMakeFiles/core.dir/flags.make
src/core/CMakeFiles/core.dir/src/Gsim.cc.o: ../src/core/src/Gsim.cc
src/core/CMakeFiles/core.dir/src/Gsim.cc.o: src/core/CMakeFiles/core.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enrg/Software/WMUtils/ratpac/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object src/core/CMakeFiles/core.dir/src/Gsim.cc.o"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/core/CMakeFiles/core.dir/src/Gsim.cc.o -MF CMakeFiles/core.dir/src/Gsim.cc.o.d -o CMakeFiles/core.dir/src/Gsim.cc.o -c /home/enrg/Software/WMUtils/ratpac/src/core/src/Gsim.cc

src/core/CMakeFiles/core.dir/src/Gsim.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/core.dir/src/Gsim.cc.i"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/enrg/Software/WMUtils/ratpac/src/core/src/Gsim.cc > CMakeFiles/core.dir/src/Gsim.cc.i

src/core/CMakeFiles/core.dir/src/Gsim.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/core.dir/src/Gsim.cc.s"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/enrg/Software/WMUtils/ratpac/src/core/src/Gsim.cc -o CMakeFiles/core.dir/src/Gsim.cc.s

src/core/CMakeFiles/core.dir/src/ConstructUserProc.cc.o: src/core/CMakeFiles/core.dir/flags.make
src/core/CMakeFiles/core.dir/src/ConstructUserProc.cc.o: ../src/core/src/ConstructUserProc.cc
src/core/CMakeFiles/core.dir/src/ConstructUserProc.cc.o: src/core/CMakeFiles/core.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enrg/Software/WMUtils/ratpac/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object src/core/CMakeFiles/core.dir/src/ConstructUserProc.cc.o"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/core/CMakeFiles/core.dir/src/ConstructUserProc.cc.o -MF CMakeFiles/core.dir/src/ConstructUserProc.cc.o.d -o CMakeFiles/core.dir/src/ConstructUserProc.cc.o -c /home/enrg/Software/WMUtils/ratpac/src/core/src/ConstructUserProc.cc

src/core/CMakeFiles/core.dir/src/ConstructUserProc.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/core.dir/src/ConstructUserProc.cc.i"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/enrg/Software/WMUtils/ratpac/src/core/src/ConstructUserProc.cc > CMakeFiles/core.dir/src/ConstructUserProc.cc.i

src/core/CMakeFiles/core.dir/src/ConstructUserProc.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/core.dir/src/ConstructUserProc.cc.s"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/enrg/Software/WMUtils/ratpac/src/core/src/ConstructUserProc.cc -o CMakeFiles/core.dir/src/ConstructUserProc.cc.s

src/core/CMakeFiles/core.dir/src/GLG4VisManager.cc.o: src/core/CMakeFiles/core.dir/flags.make
src/core/CMakeFiles/core.dir/src/GLG4VisManager.cc.o: ../src/core/src/GLG4VisManager.cc
src/core/CMakeFiles/core.dir/src/GLG4VisManager.cc.o: src/core/CMakeFiles/core.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enrg/Software/WMUtils/ratpac/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object src/core/CMakeFiles/core.dir/src/GLG4VisManager.cc.o"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/core/CMakeFiles/core.dir/src/GLG4VisManager.cc.o -MF CMakeFiles/core.dir/src/GLG4VisManager.cc.o.d -o CMakeFiles/core.dir/src/GLG4VisManager.cc.o -c /home/enrg/Software/WMUtils/ratpac/src/core/src/GLG4VisManager.cc

src/core/CMakeFiles/core.dir/src/GLG4VisManager.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/core.dir/src/GLG4VisManager.cc.i"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/enrg/Software/WMUtils/ratpac/src/core/src/GLG4VisManager.cc > CMakeFiles/core.dir/src/GLG4VisManager.cc.i

src/core/CMakeFiles/core.dir/src/GLG4VisManager.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/core.dir/src/GLG4VisManager.cc.s"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/enrg/Software/WMUtils/ratpac/src/core/src/GLG4VisManager.cc -o CMakeFiles/core.dir/src/GLG4VisManager.cc.s

src/core/CMakeFiles/core.dir/src/GLG4VEventAction.cc.o: src/core/CMakeFiles/core.dir/flags.make
src/core/CMakeFiles/core.dir/src/GLG4VEventAction.cc.o: ../src/core/src/GLG4VEventAction.cc
src/core/CMakeFiles/core.dir/src/GLG4VEventAction.cc.o: src/core/CMakeFiles/core.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enrg/Software/WMUtils/ratpac/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object src/core/CMakeFiles/core.dir/src/GLG4VEventAction.cc.o"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/core/CMakeFiles/core.dir/src/GLG4VEventAction.cc.o -MF CMakeFiles/core.dir/src/GLG4VEventAction.cc.o.d -o CMakeFiles/core.dir/src/GLG4VEventAction.cc.o -c /home/enrg/Software/WMUtils/ratpac/src/core/src/GLG4VEventAction.cc

src/core/CMakeFiles/core.dir/src/GLG4VEventAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/core.dir/src/GLG4VEventAction.cc.i"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/enrg/Software/WMUtils/ratpac/src/core/src/GLG4VEventAction.cc > CMakeFiles/core.dir/src/GLG4VEventAction.cc.i

src/core/CMakeFiles/core.dir/src/GLG4VEventAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/core.dir/src/GLG4VEventAction.cc.s"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/enrg/Software/WMUtils/ratpac/src/core/src/GLG4VEventAction.cc -o CMakeFiles/core.dir/src/GLG4VEventAction.cc.s

src/core/CMakeFiles/core.dir/src/StackingAction.cc.o: src/core/CMakeFiles/core.dir/flags.make
src/core/CMakeFiles/core.dir/src/StackingAction.cc.o: ../src/core/src/StackingAction.cc
src/core/CMakeFiles/core.dir/src/StackingAction.cc.o: src/core/CMakeFiles/core.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enrg/Software/WMUtils/ratpac/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object src/core/CMakeFiles/core.dir/src/StackingAction.cc.o"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/core/CMakeFiles/core.dir/src/StackingAction.cc.o -MF CMakeFiles/core.dir/src/StackingAction.cc.o.d -o CMakeFiles/core.dir/src/StackingAction.cc.o -c /home/enrg/Software/WMUtils/ratpac/src/core/src/StackingAction.cc

src/core/CMakeFiles/core.dir/src/StackingAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/core.dir/src/StackingAction.cc.i"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/enrg/Software/WMUtils/ratpac/src/core/src/StackingAction.cc > CMakeFiles/core.dir/src/StackingAction.cc.i

src/core/CMakeFiles/core.dir/src/StackingAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/core.dir/src/StackingAction.cc.s"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/enrg/Software/WMUtils/ratpac/src/core/src/StackingAction.cc -o CMakeFiles/core.dir/src/StackingAction.cc.s

src/core/CMakeFiles/core.dir/src/Producer.cc.o: src/core/CMakeFiles/core.dir/flags.make
src/core/CMakeFiles/core.dir/src/Producer.cc.o: ../src/core/src/Producer.cc
src/core/CMakeFiles/core.dir/src/Producer.cc.o: src/core/CMakeFiles/core.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enrg/Software/WMUtils/ratpac/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building CXX object src/core/CMakeFiles/core.dir/src/Producer.cc.o"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/core/CMakeFiles/core.dir/src/Producer.cc.o -MF CMakeFiles/core.dir/src/Producer.cc.o.d -o CMakeFiles/core.dir/src/Producer.cc.o -c /home/enrg/Software/WMUtils/ratpac/src/core/src/Producer.cc

src/core/CMakeFiles/core.dir/src/Producer.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/core.dir/src/Producer.cc.i"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/enrg/Software/WMUtils/ratpac/src/core/src/Producer.cc > CMakeFiles/core.dir/src/Producer.cc.i

src/core/CMakeFiles/core.dir/src/Producer.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/core.dir/src/Producer.cc.s"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/enrg/Software/WMUtils/ratpac/src/core/src/Producer.cc -o CMakeFiles/core.dir/src/Producer.cc.s

src/core/CMakeFiles/core.dir/src/Trajectory.cc.o: src/core/CMakeFiles/core.dir/flags.make
src/core/CMakeFiles/core.dir/src/Trajectory.cc.o: ../src/core/src/Trajectory.cc
src/core/CMakeFiles/core.dir/src/Trajectory.cc.o: src/core/CMakeFiles/core.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enrg/Software/WMUtils/ratpac/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Building CXX object src/core/CMakeFiles/core.dir/src/Trajectory.cc.o"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/core/CMakeFiles/core.dir/src/Trajectory.cc.o -MF CMakeFiles/core.dir/src/Trajectory.cc.o.d -o CMakeFiles/core.dir/src/Trajectory.cc.o -c /home/enrg/Software/WMUtils/ratpac/src/core/src/Trajectory.cc

src/core/CMakeFiles/core.dir/src/Trajectory.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/core.dir/src/Trajectory.cc.i"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/enrg/Software/WMUtils/ratpac/src/core/src/Trajectory.cc > CMakeFiles/core.dir/src/Trajectory.cc.i

src/core/CMakeFiles/core.dir/src/Trajectory.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/core.dir/src/Trajectory.cc.s"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/enrg/Software/WMUtils/ratpac/src/core/src/Trajectory.cc -o CMakeFiles/core.dir/src/Trajectory.cc.s

src/core/CMakeFiles/core.dir/src/Processor.cc.o: src/core/CMakeFiles/core.dir/flags.make
src/core/CMakeFiles/core.dir/src/Processor.cc.o: ../src/core/src/Processor.cc
src/core/CMakeFiles/core.dir/src/Processor.cc.o: src/core/CMakeFiles/core.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enrg/Software/WMUtils/ratpac/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_16) "Building CXX object src/core/CMakeFiles/core.dir/src/Processor.cc.o"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/core/CMakeFiles/core.dir/src/Processor.cc.o -MF CMakeFiles/core.dir/src/Processor.cc.o.d -o CMakeFiles/core.dir/src/Processor.cc.o -c /home/enrg/Software/WMUtils/ratpac/src/core/src/Processor.cc

src/core/CMakeFiles/core.dir/src/Processor.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/core.dir/src/Processor.cc.i"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/enrg/Software/WMUtils/ratpac/src/core/src/Processor.cc > CMakeFiles/core.dir/src/Processor.cc.i

src/core/CMakeFiles/core.dir/src/Processor.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/core.dir/src/Processor.cc.s"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/enrg/Software/WMUtils/ratpac/src/core/src/Processor.cc -o CMakeFiles/core.dir/src/Processor.cc.s

src/core/CMakeFiles/core.dir/src/SignalHandler.cc.o: src/core/CMakeFiles/core.dir/flags.make
src/core/CMakeFiles/core.dir/src/SignalHandler.cc.o: ../src/core/src/SignalHandler.cc
src/core/CMakeFiles/core.dir/src/SignalHandler.cc.o: src/core/CMakeFiles/core.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enrg/Software/WMUtils/ratpac/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_17) "Building CXX object src/core/CMakeFiles/core.dir/src/SignalHandler.cc.o"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/core/CMakeFiles/core.dir/src/SignalHandler.cc.o -MF CMakeFiles/core.dir/src/SignalHandler.cc.o.d -o CMakeFiles/core.dir/src/SignalHandler.cc.o -c /home/enrg/Software/WMUtils/ratpac/src/core/src/SignalHandler.cc

src/core/CMakeFiles/core.dir/src/SignalHandler.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/core.dir/src/SignalHandler.cc.i"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/enrg/Software/WMUtils/ratpac/src/core/src/SignalHandler.cc > CMakeFiles/core.dir/src/SignalHandler.cc.i

src/core/CMakeFiles/core.dir/src/SignalHandler.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/core.dir/src/SignalHandler.cc.s"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/enrg/Software/WMUtils/ratpac/src/core/src/SignalHandler.cc -o CMakeFiles/core.dir/src/SignalHandler.cc.s

src/core/CMakeFiles/core.dir/src/TrackInfo.cc.o: src/core/CMakeFiles/core.dir/flags.make
src/core/CMakeFiles/core.dir/src/TrackInfo.cc.o: ../src/core/src/TrackInfo.cc
src/core/CMakeFiles/core.dir/src/TrackInfo.cc.o: src/core/CMakeFiles/core.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enrg/Software/WMUtils/ratpac/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_18) "Building CXX object src/core/CMakeFiles/core.dir/src/TrackInfo.cc.o"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/core/CMakeFiles/core.dir/src/TrackInfo.cc.o -MF CMakeFiles/core.dir/src/TrackInfo.cc.o.d -o CMakeFiles/core.dir/src/TrackInfo.cc.o -c /home/enrg/Software/WMUtils/ratpac/src/core/src/TrackInfo.cc

src/core/CMakeFiles/core.dir/src/TrackInfo.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/core.dir/src/TrackInfo.cc.i"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/enrg/Software/WMUtils/ratpac/src/core/src/TrackInfo.cc > CMakeFiles/core.dir/src/TrackInfo.cc.i

src/core/CMakeFiles/core.dir/src/TrackInfo.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/core.dir/src/TrackInfo.cc.s"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/enrg/Software/WMUtils/ratpac/src/core/src/TrackInfo.cc -o CMakeFiles/core.dir/src/TrackInfo.cc.s

src/core/CMakeFiles/core.dir/src/Log.cc.o: src/core/CMakeFiles/core.dir/flags.make
src/core/CMakeFiles/core.dir/src/Log.cc.o: ../src/core/src/Log.cc
src/core/CMakeFiles/core.dir/src/Log.cc.o: src/core/CMakeFiles/core.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enrg/Software/WMUtils/ratpac/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_19) "Building CXX object src/core/CMakeFiles/core.dir/src/Log.cc.o"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/core/CMakeFiles/core.dir/src/Log.cc.o -MF CMakeFiles/core.dir/src/Log.cc.o.d -o CMakeFiles/core.dir/src/Log.cc.o -c /home/enrg/Software/WMUtils/ratpac/src/core/src/Log.cc

src/core/CMakeFiles/core.dir/src/Log.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/core.dir/src/Log.cc.i"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/enrg/Software/WMUtils/ratpac/src/core/src/Log.cc > CMakeFiles/core.dir/src/Log.cc.i

src/core/CMakeFiles/core.dir/src/Log.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/core.dir/src/Log.cc.s"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/enrg/Software/WMUtils/ratpac/src/core/src/Log.cc -o CMakeFiles/core.dir/src/Log.cc.s

src/core/CMakeFiles/core.dir/src/GLG4HitPMTCollection.cc.o: src/core/CMakeFiles/core.dir/flags.make
src/core/CMakeFiles/core.dir/src/GLG4HitPMTCollection.cc.o: ../src/core/src/GLG4HitPMTCollection.cc
src/core/CMakeFiles/core.dir/src/GLG4HitPMTCollection.cc.o: src/core/CMakeFiles/core.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enrg/Software/WMUtils/ratpac/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_20) "Building CXX object src/core/CMakeFiles/core.dir/src/GLG4HitPMTCollection.cc.o"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/core/CMakeFiles/core.dir/src/GLG4HitPMTCollection.cc.o -MF CMakeFiles/core.dir/src/GLG4HitPMTCollection.cc.o.d -o CMakeFiles/core.dir/src/GLG4HitPMTCollection.cc.o -c /home/enrg/Software/WMUtils/ratpac/src/core/src/GLG4HitPMTCollection.cc

src/core/CMakeFiles/core.dir/src/GLG4HitPMTCollection.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/core.dir/src/GLG4HitPMTCollection.cc.i"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/enrg/Software/WMUtils/ratpac/src/core/src/GLG4HitPMTCollection.cc > CMakeFiles/core.dir/src/GLG4HitPMTCollection.cc.i

src/core/CMakeFiles/core.dir/src/GLG4HitPMTCollection.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/core.dir/src/GLG4HitPMTCollection.cc.s"
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/enrg/Software/WMUtils/ratpac/src/core/src/GLG4HitPMTCollection.cc -o CMakeFiles/core.dir/src/GLG4HitPMTCollection.cc.s

core: src/core/CMakeFiles/core.dir/src/GLG4HitPhoton.cc.o
core: src/core/CMakeFiles/core.dir/src/GLG4HitPMT.cc.o
core: src/core/CMakeFiles/core.dir/src/GLG4SteppingAction.cc.o
core: src/core/CMakeFiles/core.dir/src/PythonProc.cc.o
core: src/core/CMakeFiles/core.dir/src/ProcBlock.cc.o
core: src/core/CMakeFiles/core.dir/src/RunManager.cc.o
core: src/core/CMakeFiles/core.dir/src/CountProc.cc.o
core: src/core/CMakeFiles/core.dir/src/PruneProc.cc.o
core: src/core/CMakeFiles/core.dir/src/Gsim.cc.o
core: src/core/CMakeFiles/core.dir/src/ConstructUserProc.cc.o
core: src/core/CMakeFiles/core.dir/src/GLG4VisManager.cc.o
core: src/core/CMakeFiles/core.dir/src/GLG4VEventAction.cc.o
core: src/core/CMakeFiles/core.dir/src/StackingAction.cc.o
core: src/core/CMakeFiles/core.dir/src/Producer.cc.o
core: src/core/CMakeFiles/core.dir/src/Trajectory.cc.o
core: src/core/CMakeFiles/core.dir/src/Processor.cc.o
core: src/core/CMakeFiles/core.dir/src/SignalHandler.cc.o
core: src/core/CMakeFiles/core.dir/src/TrackInfo.cc.o
core: src/core/CMakeFiles/core.dir/src/Log.cc.o
core: src/core/CMakeFiles/core.dir/src/GLG4HitPMTCollection.cc.o
core: src/core/CMakeFiles/core.dir/build.make
.PHONY : core

# Rule to build all files generated by this target.
src/core/CMakeFiles/core.dir/build: core
.PHONY : src/core/CMakeFiles/core.dir/build

src/core/CMakeFiles/core.dir/clean:
	cd /home/enrg/Software/WMUtils/ratpac/build/src/core && $(CMAKE_COMMAND) -P CMakeFiles/core.dir/cmake_clean.cmake
.PHONY : src/core/CMakeFiles/core.dir/clean

src/core/CMakeFiles/core.dir/depend:
	cd /home/enrg/Software/WMUtils/ratpac/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/enrg/Software/WMUtils/ratpac /home/enrg/Software/WMUtils/ratpac/src/core /home/enrg/Software/WMUtils/ratpac/build /home/enrg/Software/WMUtils/ratpac/build/src/core /home/enrg/Software/WMUtils/ratpac/build/src/core/CMakeFiles/core.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/core/CMakeFiles/core.dir/depend

