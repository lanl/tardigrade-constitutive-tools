# USAGE:
#
# ./new_build.sh CXX_compiler

# Make bash script more like high-level languages.
set -Eeuxo pipefail

# Get this scripts file name
script=`basename "$0"`

# Parse arguments
if [ "$#" -lt 1 ]; then
    echo "${script} USAGE:"
    echo "./${script} CXX_compiler"
    echo "    CXX_compiler: desired c++ compiler"
    echo "    CMAKE_BUILD_TYPE: cmake build type string"
    exit 1
fi
cxx_path=$1  # Path to CXX compiler
cmake_build_type='None'
if [ "$#" -ge 2 ]; then
    cmake_build_type=$2
fi
if [ "$#" -gt 2 ]; then
    echo "I only know how to use the first two arguments"
fi

# Verify compiler command and set CXX env variable
if [ -x "$(command -v ${cxx_path})" ]; then
    export CXX="$(command -v ${cxx_path})"
else
    exit 2
fi

# Debugging
whoami
groups
ls -l $HOME/include || true
ls -l $HOME/.local/include || true

# Find cmake3 executable
if [ -x "$(command -v cmake3)" ]; then
    cmake_exec=$(command -v cmake3)
elif [ -x "$(command -v cmake)" ]; then
    cmake_exec=$(command -v cmake)
else
    echo "Could not find cmake executable"
    exit 3
fi

# Clean and build repo tests
rm -rf build/
mkdir build
cd build
${cmake_exec} .. -DCMAKE_BUILD_TYPE=${cmake_build_type}
${cmake_exec} --build . --verbose
