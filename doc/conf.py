# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import pathlib
import shutil
import os
import sys
sys.path.insert(0, os.path.abspath(".."))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'DeBruijn'
copyright = '2022, Institut Pasteur'
author = 'Amine Ghozlane'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
extensions = ['myst_parser', 'sphinx.ext.autodoc','sphinx.ext.napoleon']
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# The readme that already exists
readme_path = pathlib.Path(__file__).parent.resolve().parent / "README.md"
# We copy a modified version here
readme_target = pathlib.Path(__file__).parent / "README.md"
shutil.copyfile(readme_path, readme_target)