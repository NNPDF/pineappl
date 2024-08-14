# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import pineappl

project = 'pineappl'
copyright = '2020–2024, the PineAPPL team'
author = 'the PineAPPL team'
release = pineappl.version
version = release


extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.autosummary',
    'sphinx.ext.extlinks',
    'sphinx.ext.inheritance_diagram',
    'sphinx.ext.intersphinx',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx_rtd_theme',
    'sphinxcontrib.bibtex',
]

autosummary_generate = True
autosummary_imported_members = True

extlinks = {
    "yadism": ("https://nnpdf.github.io/yadism/%s", "yadism - %s"),
    "rustdoc": ("https://docs.rs/pineappl/latest/pineappl/%s", "PineAPPL - %s"),
    "pineko": ("https://github.com/NNPDF/pineko/%s", "pineko - %s"),
}

bibtex_bibfiles = ["refs.bib"]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
