###################
constitutive\_tools
###################

*******************
Project Description
*******************

A collection of tools useful for constitutive modeling. These tools are
intended to reduce the burden in creating a new constitutive model from
scratch enabling a faster turn-around for model development. These tools
should be as general as possible to avoid cluttering the database with
extraneous things.

Information
===========

* Documentation: https://aea.re-pages.lanl.gov/material-models/constitutive_tools
* Wiki: https://re-git.lanl.gov/aea/material-models/constitutive_tools/-/wikis/home

Developers
==========

* Nathan Miller nathanm@lanl.gov
* Kyle Brindley kbrindley@lanl.gov

************
Dependencies
************

Compilers
=========

* c++11 compiler (listed version number has been tested at some point)

  * g++ >= GNU 4.8.5

Executables
===========

* [CMake](https://cmake.org/cmake/help/v3.14/) >= 3.14
* [Doxygen](https://www.doxygen.nl/manual/docblocks.html) >= 1.8.5
* [LaTeX](https://www.latex-project.org/help/documentation/) >= 2017

Python Modules (for documentation)
==================================

For convenience, the minimal Python environment requirements for the
documentation build are included in ``configuration_files/environment.yaml``.
This file was created from the [pipreqs](https://github.com/bndr/pipreqs)
command line tool and Sphinx configuration inspection, e.g. the extension
packages.

.. code-block:: bash

   $ pwd
   path/to/constitutive_tools/
   $ pipreqs --use-local --print --no-pin .

A minimal anaconda environment for building the documentation can be created
from an existing anaconda installation with the following commands.

.. code-block:: bash

   $ conda env create --file configuration_files/environment.yaml

You can learn more about Anaconda Python environment creation and management in
the [Anaconda
Documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)

C++ Libraries
=============

.. note::

   **NOTE: Non-admin installations for Eigen and Boost are no longer required.** This project is built and deployed
   against C++ libraries managed in Conda. See the Conda environment file and README discussion for non-admin environment
   management.

* [Eigen](https://eigen.tuxfamily.org/dox/) >= 3.3.7
* [BOOST](https://www.boost.org/doc/libs/1_59_0/) >= 1.59.0
* error\_tools: https://re-git.lanl.gov/aea/material-models/error_tools
* vector\_tools: https://re-git.lanl.gov/aea/material-models/vector_tools

If not found on the current system or active Conda environment, all of the
``*_tools`` libraries are pulled from their git repos by branch name and built
with their respective cmake files as part of the cmake build for this project.

**************
Build and Test
**************

This project is built with [CMake](https://cmake.org/cmake/help/v3.14/) and uses
[Sphinx](https://www.sphinx-doc.org/en/master/) to build the documentation with
[Doxygen](https://www.doxygen.nl/manual/docblocks.html) +
[Breathe](https://breathe.readthedocs.io/en/latest/) for the c++ API.

.. warning::

   **API Health Note**: The sphinx API docs are a work-in-progress. The doxygen
   API is much more useful

A build script has been created for convenience, ``new_build.sh``. It will build
everything including the library binary, the test binary, and the documentation.
This is the same build script used by ``jenkins_build.sh`` for CI builds and
testing.

sstelmo
=======

1) Activate a [W-13 Python Environment](https://xcp-confluence.lanl.gov/display/PYT/The+W-13+Python+3+environment)

   .. code-block:: bash

      $ module load python/2020.07-python-3.8
      $ sv3r

2) Build everything

   .. code-block:: bash

      $ pwd
      /path/to/constitutive_tools/

      # Just perform the build. Usage arguments are "cmake_build_type"
      ./new_build.sh None

      # Build and perform tests
      ./jenkins_build.sh

3) View test results

   .. code-block:: bash

      cat build/src/cpp/tests/results.tex

4) Display docs

   .. code-block:: bash

      # Sphinx
      firefox build/docs/sphinx/html/index.html &

      # Doxygen
      firefox build/docs/doxygen/html/index.html &

Local development
=================

In some cases it is not convenient to pull down every repository required but it may be desired that local
versions of the repository are used. An example of when this may be needed is if development is across
multiple libraries and is proceeding faster than collaborators can check in results. In this case, and
outside of developers no-one should need to do this, a version of the code using local repositories can be
built.

To perform in-source builds of upstream libraries, the active Conda environment can NOT include installed versions of
the upstream libraries to be built in-source with the current project. It is possible to mix sources with some upstream
libraries coming from the active Conda environment and others built in-source from a Git repository. Developers may
build minimal working Conda environments from the Python Modules discussion.

1) Build and activate a minimal Conda development environment

   .. code-block:: bash

       $ conda env create --file configuration_files/environment.yaml
       $ conda activate environment

2) Define convenience environment variables

   .. code-block:: bash

      $ error_tools=/path/to/my/error_tools
      $ error_tools_version=origin/dev
      $ vector_tools=/path/to/my/vector_tools
      $ vector_tools_version=origin/dev

3) Perform the initial configuration. Note that the environment variables are mutually independent. Each variable can be
   used alone or in arbitrary combinations. The default values are found in the root ``CMakeLists.txt`` file. The ``PATH``
   variables can accept anything that the [``CMake``
   ``FetchContent``](https://cmake.org/cmake/help/latest/module/FetchContent.html) ``GIT_REPOSITORY`` option can accept.
   The ``GITTAG`` variables will accept anything that the [``CMake``
   ``FetchContent``](https://cmake.org/cmake/help/latest/module/FetchContent.html) ``GIT_TAG`` option can accept.

   .. code-block:: bash

      # View the defaults
      $ grep _TOOLS_ CMakeLists.txt
      set(ERROR_TOOLS_PATH "" CACHE PATH "The path to the local version of error_tools")
      set(ERROR_TOOLS_GITTAG "" CACHE PATH "The path to the local version of error_tools")
      set(VECTOR_TOOLS_PATH "" CACHE PATH "The path to the local version of vector_tools")
      set(VECTOR_TOOLS_GITTAG "" CACHE PATH "The path to the local version of vector_tools")

      $ Build against local directory paths and possible custom branch
      $ pwd
      /path/to/constitutive_tools
      $ mkdir build
      $ cd build
      $ cmake .. -DERROR_TOOLS_PATH=${my_error_tools} -DERROR_TOOLS_GITTAG=${error_tools_version} -DVECTOR_TOOLS_PATH=${my_vector_tools} -DVECTOR_TOOLS_GITTAG=${vector_tools_version}

4) Building the library

   .. code-block:: bash

      $ pwd
      /path/to/constitutive_tools/build
      $ make

Building the documentation
==========================

To build just the documentation pick up the steps here:

2) Create the build directory and move there

   .. code-block:: bash

      $ pwd
      /path/to/constitutive_tools/
      $ mkdir build/
      $ cd build/

3) Run cmake3 configuration

   .. code-block:: bash

      $ pwd
      /path/to/constitutive_tools/build/
      $ cmake3 ..

4) Build the docs

   .. code-block:: bash

      $ cmake3 --build docs

5) Documentation builds to:

   .. code-block:: bash

      constitutive_tools/build/docs/sphinx/index.html

6) Display docs

   .. code-block:: bash

      $ pwd
      /path/to/constitutive_tools/build/
      $ firefox docs/sphinx/index.html &

7) While the Sphinx API is still a WIP, try the doxygen API

   .. code-block:: bash

      $ pwd
      /path/to/constitutive_tools/build/
      $ firefox docs/doxygen/html/index.html &

*******************
Install the library
*******************

Build the entire before performing the installation.

4) Build the entire project

   .. code-block:: bash

      $ pwd
      /path/to/constitutive_tools/build
      $ cmake3 --build .

5) Install the library

   .. code-block:: bash

      $ pwd
      /path/to/constitutive_tools/build
      $ cmake --install . --prefix path/to/root/install

      # Example local user (non-admin) Linux install
      $ cmake --install . --prefix /home/$USER/.local

      # Example install to conda environment
      $ conda active my_env
      $ cmake --install . --prefix ${CONDA_DEFAULT_ENV}

      # Example install to W-13 CI/CD conda environment performed by CI/CD institutional account
      $ cmake --install . --prefix /projects/aea_compute/release

***********************
Contribution Guidelines
***********************

Git Commit Message
==================

Begin Git commit messages with one of the following headings:

* BUG: bug fix
* DOC: documentation
* FEAT: feature
* MAINT: maintenance
* TST: tests
* REL: release
* WIP: work-in-progress

For example:

.. code-block:: bash

   git commit -m "DOC: adds documentation for feature"

Git Branch Names
================

When creating branches use one of the following naming conventions. When in
doubt use ``feature/<description>``.

* ``bugfix/\<description>``
* ``feature/\<description>``
* ``release/\<description>``

reStructured Text
=================

[Sphinx](https://www.sphinx-doc.org/en/master/) reads in docstrings and other special portions of the code as
reStructured text. Developers should follow styles in this [Sphinx style
guide](https://documentation-style-guide-sphinx.readthedocs.io/en/latest/style-guide.html#).

Style Guide
===========

This project does not yet have a full style guide. Generally, wherever a style can't be
inferred from surrounding code this project falls back to
[PEP-8](https://www.python.org/dev/peps/pep-0008/)-like styles. There are two
notable exceptions to the notional PEP-8 fall back:

1. [Doxygen](https://www.doxygen.nl/manual/docblocks.html) style docstrings are
   required for automated, API from source documentation.
2. This project prefers expansive whitespace surrounding parentheses, braces, and
   brackets.
   * No leading space between a function and the argument list.
   * One space following an open paranthesis ``(``, brace ``{``, or bracket
     ``[``
   * One space leading a close paranthesis ``)``, brace ``}``, or bracket ``]``

An example of the whitespace style:

.. code-block:: bash

   my_function( arg1, { arg2, arg3 }, arg4 );

The following ``sed`` commands may be useful for updating white space, but must
be used with care. The developer is recommended to use a unique git commit
between each command with a corresponding review of the changes and a unit test
run.

* Trailing space for open paren/brace/bracket

  .. code-block:: bash

     sed -i 's/\([({[]\)\([^ ]\)/\1 \2/g' <list of files to update>

* Leading space for close paren/brace/bracket

  .. code-block:: bash

     sed -i 's/\([^ ]\)\([)}\]]\)/\1 \2/g' <list of files to update>

* White space between adjacent paren/brace/bracket

  .. code-block:: bash

     sed -i 's/\([)}\]]\)\([)}\]]\)/\1 \2/g' <list of files to update>
