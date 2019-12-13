#!/usr/bin/env bash

# Set common shell script variables
repo='constitutive_tools'
workdir=${PWD}
declare -A deprepo 
deprepo['eigen']='https://gitlab.com/libeigen/eigen.git'
deprepo['error_tools']='ssh://git@xcp-stash.lanl.gov:7999/mm/error_tools.git'
deprepo['vector_tools']='ssh://git@xcp-stash.lanl.gov:7999/mm/vector_tools.git'
proxyout='proxyout.lanl.gov:8080'
