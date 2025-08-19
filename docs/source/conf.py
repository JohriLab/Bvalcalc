# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Bvalcalc'
copyright = '2025, Jacob Marsh, Parul Johri'
author = 'Jacob Marsh, Parul Johri'
release = '0.6.5'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'collapse_navigation': False,
    'sticky_navigation': False,
    'navigation_depth': 1,
}
html_js_files = [
    "js/onthispage-scrollspy.js",
]
html_static_path = ['_static']

templates_path = ['_templates']
html_show_sphinx = False
html_show_sourcelink = False

# html_show_copyright = False

html_css_files = [
    'custom.css',
]

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))  # so it can find Bvalcalc/
html_output = os.path.join(os.path.dirname(__file__), "..", "docs")

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "myst_parser",
]

myst_enable_extensions = [
  "deflist",
  "html_image",
  # etc.
]