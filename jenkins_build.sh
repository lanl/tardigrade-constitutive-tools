#!/usr/bin/env bash

# Identify common paths automatically
case $OSTYPE in
    darwin*)
        apps="/Users/apps"
        projects="/Users/projects"
        ;;
    linux-gnu*)
        apps="/apps"
        projects="/projects"
        ;;
    *)
        echo "Detected OS $OSTYPE is not supported. Exiting."
        exit 3
        ;;
esac

# Activate conda environment
rel='develop'
if [ -f "${apps}/anaconda/2019.10-python-3.7/etc/profile.d/conda.sh" ]; then
    . "${apps}/anaconda/2019.10-python-3.7/etc/profile.d/conda.sh"
    conda activate
    conda activate "${projects}/python/${rel}"
else
    export PATH="${apps}/anaconda/5.0.1-python-3.6/bin:$PATH"
    source activate
    source activate "${projects}/python/${rel}"
fi

# Make bash script more like high-level languages.
# https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
# Have to do this after sourcing ifortvars.sh becuase the shell script has unbound variables
set -Eeuxo pipefail

# Clone and update dependencies
source set_vars.sh

# Clean and build repo tests
./new_build.sh

# Perform repo tests
cd "build/${tests}"
./test_${repo}

# Check for failed tests
if grep -i false results.tex; then
    exit 1
fi
