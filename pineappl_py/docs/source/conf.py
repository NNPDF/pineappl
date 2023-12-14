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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

import pathlib
import os
import sys

here = pathlib.Path(__file__).absolute().parent

# -- Project information -----------------------------------------------------

project = "pineappl"
copyright = "2020-2021, the PineAPPL team"
author = "the PineAPPL team"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.ifconfig",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.napoleon",
    "sphinxcontrib.bibtex",
    "sphinx.ext.graphviz",
    "sphinx.ext.extlinks",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
source_suffix = {
    ".rst": "restructuredtext",
    ".txt": "restructuredtext",
}

autosectionlabel_prefix_document = True
# autosectionlabel_maxdepth = 10
# Allow to embed rst syntax in  markdown files.
enable_eval_rst = True

# The master toctree document.
master_doc = "index"
bibtex_bibfiles = ["refs.bib"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["shared/*"]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = None

# A string to be included at the beginning of all files
shared = here / "shared"
rst_prolog = "\n".join([open(x).read() for x in os.scandir(shared)])

extlinks = {
    "yadism": ("https://nnpdf.github.io/yadism/%s", "yadism - %s"),
    "rustdoc": ("https://docs.rs/pineappl/latest/pineappl/%s", "PineAPPL - %s"),
    "pineko": ("https://github.com/NNPDF/pineko/%s", "pineko - %s"),
}

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# -- Extension configuration -------------------------------------------------

# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
# Thanks https://github.com/bskinn/sphobjinv
intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "numpy": ("https://numpy.org/doc/stable", None),
}
# -- Options for todo extension ----------------------------------------------

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True

mathjax3_config = {
    "tex": {
        "macros": {
            # fncs
            # "atan": [r"\text{atan}", 0],
            # "span": [r"\text{span}", 0],
        }
    }
}

# https://stackoverflow.com/questions/1871549/determine-if-python-is-running-inside-virtualenv
def get_base_prefix_compat():
    """Get base/real prefix, or sys.prefix if there is none."""
    return (
        getattr(sys, "base_prefix", None)
        or getattr(sys, "real_prefix", None)
        or sys.prefix
    )


def in_virtualenv():
    return get_base_prefix_compat() != sys.prefix


# https://github.com/readthedocs/readthedocs.org/issues/1139#issuecomment-312626491
def run_apidoc(_):
    import subprocess  # pylint: disable=import-outside-toplevel

    from sphinx.ext.apidoc import main  # pylint: disable=import-outside-toplevel

    sys.path.append(str(here.parent))
    # run maturin to have the latest stuff
    pkg_root = here.parents[1]
    # if in_virtualenv():  # in local repos we're always in a virtualenv
    #     subprocess.run(["maturin", "develop"], cwd=pkg_root)
    # else:  # on RTD we can't (for some reason we're not inside the virtualenv - or maybe only the subshell isn't)
    #     subprocess.run(["maturin", "build"], cwd=pkg_root)
    #     # On RTD we were already installing before, but of course this was fake
    #     # as it only had the raw Python stuff, so let's do it again
    #     subprocess.run(["pip", "uninstall", "pineappl", "-y"], cwd=pkg_root)
    #     wheels = list((pkg_root / "target" / "wheels").glob("pineappl*.whl"))
    #     # In case there are several wheels (as on RTD) find the one matching (and let the others happily fail)
    #     for wheel in wheels:
    #         subprocess.run(["pip", "install", str(wheel.absolute())], cwd=pkg_root)

    # analyse 'pineappl'
    docs_dest = here / "modules" / "pineappl"
    import pineappl

    # note that we can NOT point to the local directory (`here.parents[1] / "pineappl"`)
    # but we need the package built by `maturin` and installed by `pip`
    package = pathlib.Path(pineappl.__file__).parent
    main(["--module-first", "-o", str(docs_dest), str(package)])
    (docs_dest / "modules.rst").unlink()


def setup(app):
    app.connect("builder-inited", run_apidoc)
