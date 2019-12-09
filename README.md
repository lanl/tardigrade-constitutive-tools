A collection of tools useful for constitutive modeling. These tools are 
intended to reduce the burden in creating a new constitutive model from 
scratch enabling a faster turn-around for model development. These tools 
should be as general as possible to avoid cluttering the database with 
extraneous things.

Note: In order to use the Intel compiler one must run the following command 
in a bash prompt:
source /apps/intel2016/bin/ifortvars.sh -arch intel64 -platform linux

This is the same command that the abaqus command issues. It may be that 
this command will change on different platforms.

---

---

Dependencies: 

These tools have several dependencies that must be available in the same parent
directory as this repo. 

* eigen: https://gitlab.com/libeigen/eigen
* error\_tools: https://xcp-stash.lanl.gov/projects/MM/repos/error_tools
* vector\_tools: https://xcp-stash.lanl.gov/projects/MM/repos/vector_tools
