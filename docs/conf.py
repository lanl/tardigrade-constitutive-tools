# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import sys
# sys.path.insert(0, os.path.abspath('.'))
import os
import re


# -- Project information -----------------------------------------------------

project = 'constitutive_tools'
copyright = '2020, Nathan A. Miller and Kyle A. Brindley'
author = 'Nathan A. Miller and Kyle A. Brindley'

# The full version, including alpha/beta/rc tags
git_describe = os.popen('git describe --always --dirty --tags').read().strip()
with open("../CMakeLists.txt") as config:
    contents = config.read()
pattern = f"{project} VERSION ([0-9]+\.[0-9]+\.[0-9])"
version_search = re.search(pattern, contents)
if version_search:
    release = version_search.group(1)
else:
    release = git_describe
if release != git_describe:
    release = release + f"+{git_describe}"
version = release


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ["breathe"]

# Breathe Configuration
breathe_projects = {project: "../build/docs/doxygen/xml"}
breathe_default_project = project

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']