# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.25.1/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.25.1/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3.

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./build

# Include any dependencies generated for this target.
include CMakeFiles/driven_QME.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/driven_QME.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/driven_QME.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/driven_QME.dir/flags.make

CMakeFiles/driven_QME.dir/main.o: CMakeFiles/driven_QME.dir/flags.make
CMakeFiles/driven_QME.dir/main.o: /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./main.cpp
CMakeFiles/driven_QME.dir/main.o: CMakeFiles/driven_QME.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/driven_QME.dir/main.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/driven_QME.dir/main.o -MF CMakeFiles/driven_QME.dir/main.o.d -o CMakeFiles/driven_QME.dir/main.o -c /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./main.cpp

CMakeFiles/driven_QME.dir/main.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/driven_QME.dir/main.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./main.cpp > CMakeFiles/driven_QME.dir/main.i

CMakeFiles/driven_QME.dir/main.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/driven_QME.dir/main.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./main.cpp -o CMakeFiles/driven_QME.dir/main.s

CMakeFiles/driven_QME.dir/constant.o: CMakeFiles/driven_QME.dir/flags.make
CMakeFiles/driven_QME.dir/constant.o: /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./constant.cpp
CMakeFiles/driven_QME.dir/constant.o: CMakeFiles/driven_QME.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/driven_QME.dir/constant.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/driven_QME.dir/constant.o -MF CMakeFiles/driven_QME.dir/constant.o.d -o CMakeFiles/driven_QME.dir/constant.o -c /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./constant.cpp

CMakeFiles/driven_QME.dir/constant.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/driven_QME.dir/constant.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./constant.cpp > CMakeFiles/driven_QME.dir/constant.i

CMakeFiles/driven_QME.dir/constant.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/driven_QME.dir/constant.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./constant.cpp -o CMakeFiles/driven_QME.dir/constant.s

CMakeFiles/driven_QME.dir/ioput.o: CMakeFiles/driven_QME.dir/flags.make
CMakeFiles/driven_QME.dir/ioput.o: /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./ioput.cpp
CMakeFiles/driven_QME.dir/ioput.o: CMakeFiles/driven_QME.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/driven_QME.dir/ioput.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/driven_QME.dir/ioput.o -MF CMakeFiles/driven_QME.dir/ioput.o.d -o CMakeFiles/driven_QME.dir/ioput.o -c /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./ioput.cpp

CMakeFiles/driven_QME.dir/ioput.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/driven_QME.dir/ioput.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./ioput.cpp > CMakeFiles/driven_QME.dir/ioput.i

CMakeFiles/driven_QME.dir/ioput.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/driven_QME.dir/ioput.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./ioput.cpp -o CMakeFiles/driven_QME.dir/ioput.s

CMakeFiles/driven_QME.dir/initialize.o: CMakeFiles/driven_QME.dir/flags.make
CMakeFiles/driven_QME.dir/initialize.o: /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./initialize.cpp
CMakeFiles/driven_QME.dir/initialize.o: CMakeFiles/driven_QME.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/driven_QME.dir/initialize.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/driven_QME.dir/initialize.o -MF CMakeFiles/driven_QME.dir/initialize.o.d -o CMakeFiles/driven_QME.dir/initialize.o -c /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./initialize.cpp

CMakeFiles/driven_QME.dir/initialize.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/driven_QME.dir/initialize.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./initialize.cpp > CMakeFiles/driven_QME.dir/initialize.i

CMakeFiles/driven_QME.dir/initialize.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/driven_QME.dir/initialize.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./initialize.cpp -o CMakeFiles/driven_QME.dir/initialize.s

CMakeFiles/driven_QME.dir/evolution.o: CMakeFiles/driven_QME.dir/flags.make
CMakeFiles/driven_QME.dir/evolution.o: /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./evolution.cpp
CMakeFiles/driven_QME.dir/evolution.o: CMakeFiles/driven_QME.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/driven_QME.dir/evolution.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/driven_QME.dir/evolution.o -MF CMakeFiles/driven_QME.dir/evolution.o.d -o CMakeFiles/driven_QME.dir/evolution.o -c /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./evolution.cpp

CMakeFiles/driven_QME.dir/evolution.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/driven_QME.dir/evolution.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./evolution.cpp > CMakeFiles/driven_QME.dir/evolution.i

CMakeFiles/driven_QME.dir/evolution.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/driven_QME.dir/evolution.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./evolution.cpp -o CMakeFiles/driven_QME.dir/evolution.s

# Object files for target driven_QME
driven_QME_OBJECTS = \
"CMakeFiles/driven_QME.dir/main.o" \
"CMakeFiles/driven_QME.dir/constant.o" \
"CMakeFiles/driven_QME.dir/ioput.o" \
"CMakeFiles/driven_QME.dir/initialize.o" \
"CMakeFiles/driven_QME.dir/evolution.o"

# External object files for target driven_QME
driven_QME_EXTERNAL_OBJECTS =

/Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./driven_QME: CMakeFiles/driven_QME.dir/main.o
/Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./driven_QME: CMakeFiles/driven_QME.dir/constant.o
/Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./driven_QME: CMakeFiles/driven_QME.dir/ioput.o
/Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./driven_QME: CMakeFiles/driven_QME.dir/initialize.o
/Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./driven_QME: CMakeFiles/driven_QME.dir/evolution.o
/Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./driven_QME: CMakeFiles/driven_QME.dir/build.make
/Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./driven_QME: /opt/homebrew/Cellar/armadillo/11.4.2/lib/libarmadillo.11.4.2.dylib
/Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./driven_QME: CMakeFiles/driven_QME.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./driven_QME"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/driven_QME.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/driven_QME.dir/build: /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./driven_QME
.PHONY : CMakeFiles/driven_QME.dir/build

CMakeFiles/driven_QME.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/driven_QME.dir/cmake_clean.cmake
.PHONY : CMakeFiles/driven_QME.dir/clean

CMakeFiles/driven_QME.dir/depend:
	cd /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3. /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3. /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./build /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./build /Users/huangzongwei/Documents/Academic/Umich/2022/Geva/driven_QME/cpp/3./build/CMakeFiles/driven_QME.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/driven_QME.dir/depend

