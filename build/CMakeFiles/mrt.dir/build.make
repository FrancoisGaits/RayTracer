# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/francois/Documents/C/IGTAI-RayTracer

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/francois/Documents/C/IGTAI-RayTracer/build

# Include any dependencies generated for this target.
include CMakeFiles/mrt.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mrt.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mrt.dir/flags.make

CMakeFiles/mrt.dir/image.cpp.o: CMakeFiles/mrt.dir/flags.make
CMakeFiles/mrt.dir/image.cpp.o: ../image.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/francois/Documents/C/IGTAI-RayTracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/mrt.dir/image.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mrt.dir/image.cpp.o -c /home/francois/Documents/C/IGTAI-RayTracer/image.cpp

CMakeFiles/mrt.dir/image.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mrt.dir/image.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/francois/Documents/C/IGTAI-RayTracer/image.cpp > CMakeFiles/mrt.dir/image.cpp.i

CMakeFiles/mrt.dir/image.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mrt.dir/image.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/francois/Documents/C/IGTAI-RayTracer/image.cpp -o CMakeFiles/mrt.dir/image.cpp.s

CMakeFiles/mrt.dir/image.cpp.o.requires:

.PHONY : CMakeFiles/mrt.dir/image.cpp.o.requires

CMakeFiles/mrt.dir/image.cpp.o.provides: CMakeFiles/mrt.dir/image.cpp.o.requires
	$(MAKE) -f CMakeFiles/mrt.dir/build.make CMakeFiles/mrt.dir/image.cpp.o.provides.build
.PHONY : CMakeFiles/mrt.dir/image.cpp.o.provides

CMakeFiles/mrt.dir/image.cpp.o.provides.build: CMakeFiles/mrt.dir/image.cpp.o


CMakeFiles/mrt.dir/kdtree.cpp.o: CMakeFiles/mrt.dir/flags.make
CMakeFiles/mrt.dir/kdtree.cpp.o: ../kdtree.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/francois/Documents/C/IGTAI-RayTracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/mrt.dir/kdtree.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mrt.dir/kdtree.cpp.o -c /home/francois/Documents/C/IGTAI-RayTracer/kdtree.cpp

CMakeFiles/mrt.dir/kdtree.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mrt.dir/kdtree.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/francois/Documents/C/IGTAI-RayTracer/kdtree.cpp > CMakeFiles/mrt.dir/kdtree.cpp.i

CMakeFiles/mrt.dir/kdtree.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mrt.dir/kdtree.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/francois/Documents/C/IGTAI-RayTracer/kdtree.cpp -o CMakeFiles/mrt.dir/kdtree.cpp.s

CMakeFiles/mrt.dir/kdtree.cpp.o.requires:

.PHONY : CMakeFiles/mrt.dir/kdtree.cpp.o.requires

CMakeFiles/mrt.dir/kdtree.cpp.o.provides: CMakeFiles/mrt.dir/kdtree.cpp.o.requires
	$(MAKE) -f CMakeFiles/mrt.dir/build.make CMakeFiles/mrt.dir/kdtree.cpp.o.provides.build
.PHONY : CMakeFiles/mrt.dir/kdtree.cpp.o.provides

CMakeFiles/mrt.dir/kdtree.cpp.o.provides.build: CMakeFiles/mrt.dir/kdtree.cpp.o


CMakeFiles/mrt.dir/lodepng-master/lodepng.cpp.o: CMakeFiles/mrt.dir/flags.make
CMakeFiles/mrt.dir/lodepng-master/lodepng.cpp.o: ../lodepng-master/lodepng.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/francois/Documents/C/IGTAI-RayTracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/mrt.dir/lodepng-master/lodepng.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mrt.dir/lodepng-master/lodepng.cpp.o -c /home/francois/Documents/C/IGTAI-RayTracer/lodepng-master/lodepng.cpp

CMakeFiles/mrt.dir/lodepng-master/lodepng.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mrt.dir/lodepng-master/lodepng.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/francois/Documents/C/IGTAI-RayTracer/lodepng-master/lodepng.cpp > CMakeFiles/mrt.dir/lodepng-master/lodepng.cpp.i

CMakeFiles/mrt.dir/lodepng-master/lodepng.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mrt.dir/lodepng-master/lodepng.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/francois/Documents/C/IGTAI-RayTracer/lodepng-master/lodepng.cpp -o CMakeFiles/mrt.dir/lodepng-master/lodepng.cpp.s

CMakeFiles/mrt.dir/lodepng-master/lodepng.cpp.o.requires:

.PHONY : CMakeFiles/mrt.dir/lodepng-master/lodepng.cpp.o.requires

CMakeFiles/mrt.dir/lodepng-master/lodepng.cpp.o.provides: CMakeFiles/mrt.dir/lodepng-master/lodepng.cpp.o.requires
	$(MAKE) -f CMakeFiles/mrt.dir/build.make CMakeFiles/mrt.dir/lodepng-master/lodepng.cpp.o.provides.build
.PHONY : CMakeFiles/mrt.dir/lodepng-master/lodepng.cpp.o.provides

CMakeFiles/mrt.dir/lodepng-master/lodepng.cpp.o.provides.build: CMakeFiles/mrt.dir/lodepng-master/lodepng.cpp.o


CMakeFiles/mrt.dir/main.cpp.o: CMakeFiles/mrt.dir/flags.make
CMakeFiles/mrt.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/francois/Documents/C/IGTAI-RayTracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/mrt.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mrt.dir/main.cpp.o -c /home/francois/Documents/C/IGTAI-RayTracer/main.cpp

CMakeFiles/mrt.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mrt.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/francois/Documents/C/IGTAI-RayTracer/main.cpp > CMakeFiles/mrt.dir/main.cpp.i

CMakeFiles/mrt.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mrt.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/francois/Documents/C/IGTAI-RayTracer/main.cpp -o CMakeFiles/mrt.dir/main.cpp.s

CMakeFiles/mrt.dir/main.cpp.o.requires:

.PHONY : CMakeFiles/mrt.dir/main.cpp.o.requires

CMakeFiles/mrt.dir/main.cpp.o.provides: CMakeFiles/mrt.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/mrt.dir/build.make CMakeFiles/mrt.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/mrt.dir/main.cpp.o.provides

CMakeFiles/mrt.dir/main.cpp.o.provides.build: CMakeFiles/mrt.dir/main.cpp.o


CMakeFiles/mrt.dir/raytracer.cpp.o: CMakeFiles/mrt.dir/flags.make
CMakeFiles/mrt.dir/raytracer.cpp.o: ../raytracer.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/francois/Documents/C/IGTAI-RayTracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/mrt.dir/raytracer.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mrt.dir/raytracer.cpp.o -c /home/francois/Documents/C/IGTAI-RayTracer/raytracer.cpp

CMakeFiles/mrt.dir/raytracer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mrt.dir/raytracer.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/francois/Documents/C/IGTAI-RayTracer/raytracer.cpp > CMakeFiles/mrt.dir/raytracer.cpp.i

CMakeFiles/mrt.dir/raytracer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mrt.dir/raytracer.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/francois/Documents/C/IGTAI-RayTracer/raytracer.cpp -o CMakeFiles/mrt.dir/raytracer.cpp.s

CMakeFiles/mrt.dir/raytracer.cpp.o.requires:

.PHONY : CMakeFiles/mrt.dir/raytracer.cpp.o.requires

CMakeFiles/mrt.dir/raytracer.cpp.o.provides: CMakeFiles/mrt.dir/raytracer.cpp.o.requires
	$(MAKE) -f CMakeFiles/mrt.dir/build.make CMakeFiles/mrt.dir/raytracer.cpp.o.provides.build
.PHONY : CMakeFiles/mrt.dir/raytracer.cpp.o.provides

CMakeFiles/mrt.dir/raytracer.cpp.o.provides.build: CMakeFiles/mrt.dir/raytracer.cpp.o


CMakeFiles/mrt.dir/scene.cpp.o: CMakeFiles/mrt.dir/flags.make
CMakeFiles/mrt.dir/scene.cpp.o: ../scene.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/francois/Documents/C/IGTAI-RayTracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/mrt.dir/scene.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mrt.dir/scene.cpp.o -c /home/francois/Documents/C/IGTAI-RayTracer/scene.cpp

CMakeFiles/mrt.dir/scene.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mrt.dir/scene.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/francois/Documents/C/IGTAI-RayTracer/scene.cpp > CMakeFiles/mrt.dir/scene.cpp.i

CMakeFiles/mrt.dir/scene.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mrt.dir/scene.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/francois/Documents/C/IGTAI-RayTracer/scene.cpp -o CMakeFiles/mrt.dir/scene.cpp.s

CMakeFiles/mrt.dir/scene.cpp.o.requires:

.PHONY : CMakeFiles/mrt.dir/scene.cpp.o.requires

CMakeFiles/mrt.dir/scene.cpp.o.provides: CMakeFiles/mrt.dir/scene.cpp.o.requires
	$(MAKE) -f CMakeFiles/mrt.dir/build.make CMakeFiles/mrt.dir/scene.cpp.o.provides.build
.PHONY : CMakeFiles/mrt.dir/scene.cpp.o.provides

CMakeFiles/mrt.dir/scene.cpp.o.provides.build: CMakeFiles/mrt.dir/scene.cpp.o


# Object files for target mrt
mrt_OBJECTS = \
"CMakeFiles/mrt.dir/image.cpp.o" \
"CMakeFiles/mrt.dir/kdtree.cpp.o" \
"CMakeFiles/mrt.dir/lodepng-master/lodepng.cpp.o" \
"CMakeFiles/mrt.dir/main.cpp.o" \
"CMakeFiles/mrt.dir/raytracer.cpp.o" \
"CMakeFiles/mrt.dir/scene.cpp.o"

# External object files for target mrt
mrt_EXTERNAL_OBJECTS =

mrt: CMakeFiles/mrt.dir/image.cpp.o
mrt: CMakeFiles/mrt.dir/kdtree.cpp.o
mrt: CMakeFiles/mrt.dir/lodepng-master/lodepng.cpp.o
mrt: CMakeFiles/mrt.dir/main.cpp.o
mrt: CMakeFiles/mrt.dir/raytracer.cpp.o
mrt: CMakeFiles/mrt.dir/scene.cpp.o
mrt: CMakeFiles/mrt.dir/build.make
mrt: CMakeFiles/mrt.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/francois/Documents/C/IGTAI-RayTracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking CXX executable mrt"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mrt.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mrt.dir/build: mrt

.PHONY : CMakeFiles/mrt.dir/build

CMakeFiles/mrt.dir/requires: CMakeFiles/mrt.dir/image.cpp.o.requires
CMakeFiles/mrt.dir/requires: CMakeFiles/mrt.dir/kdtree.cpp.o.requires
CMakeFiles/mrt.dir/requires: CMakeFiles/mrt.dir/lodepng-master/lodepng.cpp.o.requires
CMakeFiles/mrt.dir/requires: CMakeFiles/mrt.dir/main.cpp.o.requires
CMakeFiles/mrt.dir/requires: CMakeFiles/mrt.dir/raytracer.cpp.o.requires
CMakeFiles/mrt.dir/requires: CMakeFiles/mrt.dir/scene.cpp.o.requires

.PHONY : CMakeFiles/mrt.dir/requires

CMakeFiles/mrt.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mrt.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mrt.dir/clean

CMakeFiles/mrt.dir/depend:
	cd /home/francois/Documents/C/IGTAI-RayTracer/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/francois/Documents/C/IGTAI-RayTracer /home/francois/Documents/C/IGTAI-RayTracer /home/francois/Documents/C/IGTAI-RayTracer/build /home/francois/Documents/C/IGTAI-RayTracer/build /home/francois/Documents/C/IGTAI-RayTracer/build/CMakeFiles/mrt.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mrt.dir/depend

