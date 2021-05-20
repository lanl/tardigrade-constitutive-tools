# constitutive\_tools

A collection of tools useful for constitutive modeling. These tools are
intended to reduce the burden in creating a new constitutive model from
scratch enabling a faster turn-around for model development. These tools
should be as general as possible to avoid cluttering the database with
extraneous things.

Note: In order to use the Intel compiler one must run the following command in a
bash prompt:

    source /apps/intel2016/bin/ifortvars.sh -arch intel64 -platform linux

This is the same command that the abaqus command issues. It may be that
this command will change on different platforms.

---

## Dependencies:

### Compilers

* c++11 compiler (listed version number has been tested at some point)

  * g++ >= GNU 4.8.5

### Executables

* [CMake](https://cmake.org/cmake/help/v3.14/) >= 3.14
* [Doxygen](https://www.doxygen.nl/manual/docblocks.html) >= 1.8.5
* [LaTeX](https://www.latex-project.org/help/documentation/) >= 2017

### Python Modules (for documentation)

For convenience, the minimal Python environment requirements for the
documentation build are included in ``configuration_files/environment.yaml``.
This file was created from the [pipreqs](https://github.com/bndr/pipreqs)
command line tool and Sphinx configuration inspection, e.g. the extension
packages.

    $ pwd
    path/to/cpp_stub/
    $ pipreqs --use-local --print --no-pin .

A minimal anaconda environment for building the documentation can be created
from an existing anaconda installation with the following commands.

    $ conda env create --file environment.yaml

You can learn more about Anaconda Python environment creation and management in
the [Anaconda
Documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)

### C++ Libraries

> **NOTE: Non-admin installations for Eigen and Boost are no longer required.** This project is built and deployed
> against C++ libraries managed in Conda. See the Conda environment file and README discussion for non-admin environment
> management.

* [Eigen](https://eigen.tuxfamily.org/dox/) >= 3.3.7
* [BOOST](https://www.boost.org/doc/libs/1_53_0/) >= 1.53.0

### "Internal" project libraries

* error\_tools: https://xcp-stash.lanl.gov/projects/MM/repos/error_tools
* vector\_tools: https://xcp-stash.lanl.gov/projects/MM/repos/vector_tools

All of the ``{error,vector}_tools`` libraries are pulled from their git repos by
branch name and built with their respective cmake files as part of the cmake
build for this project.

---

## Build and Test

This repository is now built completely with cmake.

> **API Health Note**: The sphinx API docs are a work-in-progress. The doxygen
> API is much more useful

A build script has been created for convenience, ``new_build.sh``. It will build
everything including the library binary, the test binary, and the documentation.
This is the same build script used by ``jenkins_build.sh`` for CI builds and
testing.

### sstelmo

1) Activate a [W-13 Python Environment](https://xcp-confluence.lanl.gov/display/PYT/The+W-13+Python+3+environment)

       $ module load python/2019.10-python-3.7
       $ sv3r

2) Build everything

       $ pwd
       /path/to/constitutive_tools/
       
       # Just perform the build
       ./new_build.sh
       
       # Build and perform tests
       ./jenkins_build.sh

3) View test results

       cat build/src/cpp/tests/results.tex

4) Display docs

       # Sphinx
       firefox build/docs/sphinx/index.html &
       
       # Doxygen
       firefox build/docs/doxygen/html/index.html &

### Local development

In some cases it is not convenient to pull down every repository required but it may be desired that local
versions of the repository are used. An example of when this may be needed is if development is across
multiple libraries and is proceeding faster than collaborators can check in results. In this case, and
outside of developers no-one should need to do this, a version of the code using local repositories can be
built.

1) Activate a [W-13 Python Environment](https://xcp-confluence.lanl.gov/display/PYT/The+W-13+Python+3+environment)

       $ module load python/2019.10-python-3.7
       $ sv3r

2) Define convenience environment variables

       $ my_error_tools=/path/to/my/error_tools
       $ my_vector_tools=/path/to/my/vector_tools

3) Perform the initial configuration

       $ pwd
       /path/to/constitutive_tools
       $ mkdir build
       $ cd build
       $ cmake .. -DFETCH_SOURCE=LOCAL -DERROR_TOOLS_PATH=${my_error_tools} -DVECTOR_TOOLS_PATH=${my_vector_tools}

4) Building the library

       $ pwd
       /path/to/constitutive_tools/build
       $ make

### Building the documentation

To build just the documentation pick up the steps here:

2) Create the build directory and move there

       $ pwd
       /path/to/constitutive_tools/
       $ mkdir build/
       $ cd build/

3) Run cmake3 configuration

       $ pwd
       /path/to/constitutive_tools/build/
       $ cmake3 ..

4) Build the docs

       $ cmake3 --build docs

5) Documentation builds to:

       constitutive_tools/build/docs/sphinx/index.html

6) Display docs

       $ pwd
       /path/to/constitutive_tools/build/
       $ firefox docs/sphinx/index.html &

7) While the Sphinx API is still a WIP, try the doxygen API

       $ pwd
       /path/to/constitutive_tools/build/
       $ firefox docs/doxygen/html/index.html &
