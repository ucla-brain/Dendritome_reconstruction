## Mouse Connectome Project 3D

-----------------------------------------
Code base for large 3D image data processing. Its layout is as following.   
* **mcp3d_external_repos**
  * gtest
  * nlohman_json
  * eigen

* **mcp3d_external_packages**
  * hdf5-1.10.1.tar.bz2
  * boost_1_65_1
  
* **mcp3d**
```
  * src
        * 3rd_party: external dependencies. see `utils/built_ext_libraries.py` for details
            * gtest
            * nlohman_json
            * boost
            * hdf5 (serial version)
            * eigen 
            * opencv #TODO currently built out of src_dir
        * benchmark
        * algorithm
        * common
        * external_tools
        * image: image IO module, supports cluster parallelism.
    
            MImage class: data container
            MImageIO class: dispatches io operations requested by MImage class to 
                            derived classes of MImageFormats
            MImageFormats: abstract class with pure virtual functions
               MTiffFormat: for sequence of tiff images i/o
               #TODO MOMETiffFormat
               MHdf5: reading Imaris. #TODO: virtual hdf5 io
               #TODO anything we may need
            each image is consisted of raw voxels and available resolution levels. 
            within a resolution level, all files must have same format.
            across resolution levels, file formats can differ. one can have scans 
            saved as .raw by a microscope, and build .hdf5 pyramids for the same data. 
            MImage specifies the global spatial location and resolution level of data 
            to retrieve within the global image volume. Image file formats are opaque 
            to MImage. Format specific IO handling is done by MImageIO and MImageFormats.
            This allows concise high level API and avoids data conversions.
    
        * log: for me (mzhu) is the destination for symlink /ifshome/mzhu/log, 
          where SGE output files $job_name_$job_id.o/e are written to. only relevant 
          for LONI cluster jobs
        * parallel
        * python: python wrapper generation
        * test: unit tests, parallel environment tests
        * util: build scripts, centos system related scripts
        * some scattered files that will be cleaned up in due course
  * python: python scripts
  * test_data: test data files  
```
-------------------------------------------
#### building external dependencies
The build sections are only relevant for developers in linux builds. The 
directories `mcp3d_external_packages` and `mcp3d_external_repos` should contain
downloaded external package source code or external package repo respectively.
The build script will attempt to clone repos. If error message "Permission 
denied (publickey)", refer to resources on ssh authentication 
(e.g. https://www.git-tower.com/learn/git/ebook/en/command-line/advanced-topics/ssh-public-keys).
The external packages need to be downloaded initially. Detailed package version
information can be found in `mcp3d/src/util/build_ext_libraries.py`  
Third party dependencies build command is as following:  

```python util/build_ext_libraries.py all / [components]```

-------------------------------------------------------
#### build c++ source on local machine:
Assuming the c++ source directory is `src_dir` (`src_dir = /path/to/mcp3d/src`), and the build directory is 
`build_dir` from now on. passing `-DMPI=TRUE` to cmake will build components
that require MPI. This assumes MPI library is available on target
system. if `-DMPI=FALSE` or is omitted, no MPI components will be built.   

Passing `-DPYTHON=TRUE` will build python binding. By default, python2 binding is 
made as module `pymcp3d`. Passing `-DPYTHON_VERSION=3` will create python3 module 
`pymcp3d3`. The only supported `python2` version is `python2.7`. For `python3`,
the build searches for system availability of `python3.7`, `python3.6`, 
`python3.5`, in that order, and builds wrapper for the first interpreter version 
found. For `python3` wrappers, string literals should be prefixed 
by `'b'` to pass bytestring to c++ functions (e.g. `b'foo'`).   

Other user prameters 
to cmake are `-DVERBOSE=TRUE`, which prompts verbose output during run time; 
`-DPROFILE=TRUE`, which compiles targets with -pg flag for profiling.   

Example commands to build and install c++ source:
```
# assuming 3rd party dependencies have been built
cd build_dir
cmake src_dir -DMPI=[TRUE/FALSE] -DCMAKE_BUILD_TYPE=[Debug/Release]  -DPYTHON=[TRUE/FALSE]
make
make install
```
-------------------------------------------------------
The executables are installed to `src_dir/../bin`, while libraries are installed
to `src_dir/../lib`

#### build on loni cluster:
Done by calling `util/grid_build.sh`, which sets `MPI` to TRUE and 
builds all MPI components. Python binding is currently not made, but will be 
supported in the future. 
```
# you need proper LD_LIBRARY_PATH for gcc-7.2 and mpich-3.1.4. 
# neither the cluster default gcc nor the cluster default MPI implementation 
# compiles this project
# assuming 3rd party dependencies have been built
bash util/grid_build.sh
```
-----------------------------------------------------------
####todo: packaging for linux and windows

-----------------------------------------------------------
#### unit tests.
```
# either local or cluster execution will do
cd build_dir/test
./mcp3d_unit_tests
```
The unit tests use gtest framework and report success / failure across all tests
in test cases, like this:
```
[ RUN      ] TiffFormat.ReadPartialTiff
executing TiffFormat_ReadPartialTiffImage_Impl 5 times
time to cv::imread full image with dimensions [8192, 8192]: 0.88975 seconds
loading selected region in image. will need 0.000745319 GB of memory.
loading selected region in image. will need 0.000745319 GB of memory.
...(some std output...)...
...(some std output...)...
[       OK ] TiffFormat.ReadPartialTiff (21571 ms)
[----------] 3 tests from TiffFormat (21574 ms total)

[----------] 1 test from MImage
[ RUN      ] MImage.ReadInfoTiffFormat
[       OK ] MImage.ReadInfoTiffFormat (1 ms)
[----------] 1 test from MImage (1 ms total)

[----------] Global test environment tear-down
[==========] 21 tests from 5 test cases ran. (21587 ms total)
[  PASSED  ] 21 tests.

```
--------------------------------------------------------
#### mpi tests
if buildt with `-DMPI=TRUE`, `mcp_mpi_tests` will be created. 
`proc_number` instances of tests will be run. 
they should all pass for the entire test to be successful.
```
cd build_dir/test
mpirun/mpiexec -n proc_number mcp3d_mpi_tests
```
---------------------------------------------------
#####todo gcov lcov for coverage stats

-----------------------------------------------------

#### parallel environment and parallel io tests. 
```
cd src_dir/test/infrastructure
# simple mpi environment test
bash qsub_mpi_simple_test.sh n_proc
# mpio test
bash qsub_mpi_io_test.sh n_proc
```
the command given above as is should be run on loni cluster, though for each
`qsub_[test_name].sh` file, you should be able to invoke the test on your local
machine with following if you have a local MPI install
```
mpirun/mpiexec -n n_procs build_dir/test/infrastructure/[test_name] *args
```
--------------------------------------------------------------------------------
