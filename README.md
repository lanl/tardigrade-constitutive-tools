# constitutive\_tools

A collection of tools useful for constitutive modeling. These tools are
intended to reduce the burden in creating a new constitutive model from
scratch enabling a faster turn-around for model development. These tools
should be as general as possible to avoid cluttering the database with
extraneous things.

Note: In order to use the Intel compiler one must run the following command in a
bash prompt:

```
source /apps/intel2016/bin/ifortvars.sh -arch intel64 -platform linux
```

This is the same command that the abaqus command issues. It may be that
this command will change on different platforms.

---

---

## Dependencies:

### Executables

* CMake >= 3.14
* Doxygen >= 1.8.5

### Python Modules (for documentation)

* Sphinx >= 3.0.4
* Breathe >= 4.18.1
* sphinx\_rtd\_theme >= 0.4.3

For convenience, the minimal Python environment requirements for the
documentation build are included in ``environment.yaml`` and
``requirements.txt``. A minimal anaconda environment for building the
documentation can be created from an existing anaconda installation with the
following commands.

```
$ conda env create --file environment.yaml
```

### Libraries

* eigen >= 3.3.7
* BOOST >= 1.53.0
* error\_tools: https://xcp-stash.lanl.gov/projects/MM/repos/error_tools
* vector\_tools: https://xcp-stash.lanl.gov/projects/MM/repos/vector_tools

#### "Internal" project libraries

All of the ``{error,vector}_tools`` libraries are pulled from their git repos by
branch name and built with their respective cmake files as part of the cmake
build for this project.

#### Eigen

https://gitlab.com/libeigen/eigen

Eigen must be "installed" following the ``eigen/INSTALL`` instructions. The
Eigen dependence is easiest to resolve if eigen is installed in the default
install directory.  However, if you don't have admin privileges, you can also
insall Eigen to your home directory in ``~/include`` (or possibly in
``~/.local/include``, but this is untested by this project).

#### Non-admin Eigen install for constitutive_tools
[Reference](https://unix.stackexchange.com/questions/36871/where-should-a-local-executable-be-placed)

```
# sstelmo
ssh -X sstelmo.lanl.gov
# (OPTIONAL) source Intel compilers. c++/g++ compilers from GNU 4.8.5 also work.
source /apps/intel2016/bin/ifortvars.sh -arch intel64 -platform linux
# Create personal include file directory
$ pwd
/home/$USER
$ mkdir .local/include
# Move to repository directory
$ cd /preferred/path/to/repos
# Example
$ pwd
/projects/$USER/e13repos
# Clone eigen
$ git clone https://gitlab.com/libeigen/eigen.git
$ cd eigen
$ git checkout 3.3.7
# Build eigen
$ mkdir build
$ cd build
$ export CXX=$(command -v icpc) # OPTIONAL
$ cmake3 .. -DCMAKE_INSTALL_PREFIX=$HOME/.local
$ make install
```

---

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

```
$ module load python/2019.10-python-3.7
$ sv3r
```

2) Build everything

```
$ pwd
/path/to/constitutive_tools/

# Just perform the build
./new_build.sh

# Build and perform tests
./jenkins_build.sh
```

3) View test results

```
cat build/src/cpp/tests/results.tex
```

4) Display docs

```
# Sphinx
firefox build/docs/sphinx/index.html &

# Doxygen
firefox build/docs/doxygen/html/index.html &
```

### Building the documentation

To build just the documentation pick up the steps here:

2) Create the build directory and move there

```
$ pwd
/path/to/constitutive_tools/
$ mkdir build/
$ cd build/
```

3) Run cmake3 configuration

```
$ pwd
/path/to/constitutive_tools/build/
$ cmake3 ..
```

4) Build the docs

```
$ cmake3 --build docs
```

5) Documentation builds to:

```
constitutive_tools/build/docs/sphinx/index.html
```

6) Display docs

```
$ pwd
/path/to/constitutive_tools/build/
$ firefox docs/sphinx/index.html &
```

7) While the Sphinx API is still a WIP, try the doxygen API

```
$ pwd
/path/to/constitutive_tools/build/
$ firefox docs/doxygen/html/index.html &
```
