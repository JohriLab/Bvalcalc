# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Bvalcalc'
copyright = '2025, Jacob Marsh, Parul Johri'
author = 'Jacob Marsh, Parul Johri'
release = '0.1.2'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))  # so it can find bvalcalc/

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",   # to support Google/NumPy style docstrings
    "sphinx.ext.viewcode",   # adds [source] links
]
