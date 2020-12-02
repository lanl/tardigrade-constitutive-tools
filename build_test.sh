#!/usr/bin/env bash

# report conda environment
conda info

# Set some common shell variables
source set_vars.sh

# Clean and build repo tests
case $OSTYPE in
    darwin*)
        compiler='c++'
        ;;
    linux-gnu*)
        compiler='g++'
        ;;
esac
./new_build.sh ${compiler}

# Perform repo tests
cd "build"
ctest --verbose --output-log results.tex

#=============================================================== RUN PYTESTS ===
cd ${workdir}
./pytest_${repo}.sh
cd ${workdir}
