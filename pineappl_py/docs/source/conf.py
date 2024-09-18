# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import pineappl
import sys
import pathlib

project = 'pineappl'
copyright = '2020–2024, the PineAPPL team'
author = 'the PineAPPL team'
release = pineappl.version
version = release


extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.extlinks',
    'sphinx.ext.inheritance_diagram',
    'sphinx.ext.intersphinx',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx_rtd_theme',
    'nbsphinx',
]


extlinks = {
    "rustdoc": ("https://docs.rs/pineappl/latest/pineappl/%s", "PineAPPL - %s"),
}

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']


# TODO: find a way to reactivate apidoc, which doesn't seem to work for the moment.

# here = pathlib.Path(__file__).absolute().parent
# # https://github.com/readthedocs/readthedocs.org/issues/1139#issuecomment-312626491
# def run_apidoc(_):
#     from sphinx.ext.apidoc import main  # pylint: disable=import-outside-toplevel

#     sys.path.append(str(here.parent))
#     # analyse 'pineappl'
#     docs_dest = here / "modules"
#     import pineappl   # pylint: disable=import-outside-toplevel

#     # note that we can NOT point to the local directory (`here.parents[1] / "pineappl"`)
#     # but we need the package built by `maturin` and installed by `pip`
#     package = pathlib.Path(pineappl.__file__).parent / "pineappl"
#     main(["--module-first", "--no-toc", "-o", str(docs_dest), str(package)])


# def setup(app):
#     app.connect("builder-inited", run_apidoc)
