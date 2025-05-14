# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Bvalcalc'
copyright = '2025, Jacob Marsh, Parul Johri'
author = 'Jacob Marsh, Parul Johri'
release = '0.2.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']


html_css_files = [
    'custom.css',
]

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))  # so it can find bvalcalc/

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinxarg.ext",  # sphinx-argparse extension
    "myst_parser",
]

myst_enable_extensions = [
  "deflist",
  "html_image",
  # etc.
]