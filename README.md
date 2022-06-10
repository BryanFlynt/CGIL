# CGIL
Computational Geometry Interpolation Library

## Prerequisite
- C++ Compiler supporting C++17 or newer
- [CMake](https://cmake.org/) >= 3.19.0
- [Boost C++](https://www.boost.org/) headers only portion
- [CGAL](https://www.cgal.org/) headers only portion
  - Automaticaly downloaded during Build phase

## Configuration
1. Create a build directory (for instance inside the source tree): 
```console
cd CGIL
mkdir build
cd build
```

2. Configure your system to have the proper libraries visable to the CMake build system:  
**Note:** These rules and variables used by CMake and not part of the CGIL application
  - Boost 
    - CMake will check BOOST_ROOT then system paths for Boost libraries
```console
> export BOOST_ROOT=/my/dir/to/boost/1.71.0     # Boost Installation
```

3. Run cmake to generate the Makefile:
```console
cmake ..
```

## Build
Building the application has one extra complexity the way it's currently configured.
If the CGAL library is not already found it will download and configure CGAL within 
your build directory.  Therefore, you will need <ins>internet access</ins> and ports that allow 
access to the GitHub repository of CGAL.

To build the command is the same.
```console
> make
```

## Installation
Standard process to install into the bin directory
```console
> make install
```

## How to Use
At this stage several programs have been installed into the bin directory
```console
> cd bin
> ls
cgil2  cgil3  make_random_sources  make_random_targets
```
We can use these to see the program behavior and expected file formats

### Create random source values
To create random source values we can use the "make_random_sources" application.
The program accepts 3 command line arguments or will defualt them if not provided.
```console
make_random_sources <num_dim> <num_points> <num_variables>
```

Here we will create 1000 3D points with 2 made up variables and save to a file name "source.txt"
```console
> ./make_random_sources 3 1000 2 > source.txt
```

### Create random target locations
To create random target locations we can use the "make_random_targets" application.
The program accepts 2 command line arguments or will defualt them if not provided.
```console
make_random_targets <num_dim> <num_points>
```

Here we will create 100 3D points and save to a file name "target.txt"
```console
> ./make_random_targets 3 1000 > target.txt
```

### Run Interpolations
To interpolate source values to target locations we use the "cgil2" or "cgil3" programs depending 
on 2D or 3D data.  Again we can provide 2 command line arguments which in this case are the source 
and target file names.  If none are provided it will expect the default names which we just used in 
the previous steps.

```console
cgil2 <source_file> <target_file>
cgil3 <source_file> <target_file>
```

In this case we used the default names "source.txt" and "target.txt" so we can run without arguments 
and save the answers to a file named "interps.txt"
```console
> ./cgil3 > interps.txt
```

## Potential Issues
1. Interpolated values are 0's

If a target location is not found inside the convex hull of the source locations this will occur.  
In fact, special steps were taken inside the "make_random_" applications to ensure all targets would 
be inside the convex hull of the sources.  

2. Others?



