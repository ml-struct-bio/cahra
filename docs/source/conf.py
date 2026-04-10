# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# to add complexity to displayed api docs can extend standard templates
# via similar to https://stackoverflow.com/questions/2701998/
# automatically-document-all-modules-recursively-with-sphinx-autodoc
import os
import sys

sys.path.insert(
    0, os.path.abspath("../../src")
)  # Source code dir relative to this file


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "cahra-internal"
copyright = "Joel Greer, Tom Burnley, David Silva Sanchez, Pilar Cossio, \
    Sonya Hanson, Geoffrey Woollard, Robert Heeter, Ryan Feathers, \
    Ellen Zhong, E.Z. Lab @ Princeton"
author = "Joel Greer, Tom Burnley, David Silva Sanchez, Pilar Cossio, \
    Sonya Hanson, Geoffrey Woollard, Robert Heeter, Ryan Feathers, \
    Ellen Zhong, E.Z. Lab @ Princeton"

release = "0.0.1"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
]
autosummary_generate = True  # Turn on sphinx.ext.autosummary
autoclass_content = "both"  # Add __init__ doc (ie. params) to class summaries
html_show_sourcelink = (
    False  # Remove 'view source code' from top of page (for html, not python)
)
autodoc_inherit_docstrings = True  # If no docstring, inherit from base class
set_type_checking_flag = (
    True  # Enable 'expensive' imports for sphinx_autodoc_typehints
)

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master/", None),
}
intersphinx_disabled_domains = ["std"]

templates_path = ["_templates"]
exclude_patterns: list = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "alabaster"
html_static_path = ["_static"]
