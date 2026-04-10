# Configuration file for the Sphinx documentation builder.

# -- Project information

project = "Roodmus"
copyright = "2026, Joel Greer(UKRI), Tom Burnley (UKRI),\
    David Silva Sanchez (Flatiron Institute), Pilar Cossio (Flatiron Institute)"
author = "Joel Greer, David Silva Sanchez, Tom Burnley, Pilar Cossio"

release = "0.1"
version = "0.0.1"

# -- General configuration

extensions = [
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master/", None),
}
intersphinx_disabled_domains = ["std"]

templates_path = ["_templates"]

# -- Options for HTML output

html_theme = "sphinx_rtd_theme"

# -- Options for EPUB output
epub_show_urls = "footnote"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
