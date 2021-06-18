# Make bash script more like high-level languages.
set -Eeuxo pipefail
# Find cmake3 executable
if [ -x "$(command -v cmake3)" ]; then
    cmake_exec=$(command -v cmake3)
elif [ -x "$(command -v cmake)" ]; then
    cmake_exec=$(command -v cmake)
else
    echo "Could not find cmake executable"
    exit 3
fi
# Change to build directory and run cmake install
cd "build"
${cmake_exec} --install . --prefix ${CONDA_DEFAULT_ENV}
