import os
import pathlib
import shutil
import sys

sys.path.insert(0, os.path.abspath(".."))

import gflex

project = "gFlex"
author = "Andrew D. Wickert"
release = gflex.__version__
version = release
copyright = f"2010–2026, {author}"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "myst_parser",
]

autodoc_default_options = {
    "member-order": "bysource",
    "show-inheritance": True,
    "undoc-members": False,
}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy", None),
}

napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_use_param = True
napoleon_use_rtype = True

def _copy_changelog(app):
    shutil.copy(
        pathlib.Path(app.srcdir).parent / "CHANGES.md",
        pathlib.Path(app.srcdir) / "releases.md",
    )


def setup(app):
    app.connect("builder-inited", _copy_changelog)


html_theme = "furo"
html_static_path = ["_static"]
html_logo = "_static/logo.png"
