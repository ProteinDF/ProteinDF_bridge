"""Sphinx configuration for ProteinDF_bridge documentation."""

import os
import sys

sys.path.insert(0, os.path.abspath(".."))

# -- Project information -----------------------------------------------

project = "ProteinDF_bridge"
copyright = "2014, The ProteinDF development team"
author = "Toshiyuki HIRANO"

try:
    from importlib.metadata import version as _pkg_version

    release = _pkg_version("proteindf_bridge")
except Exception:
    release = "2024.3.0"
version = release

# -- General configuration -----------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx_copybutton",
    "sphinx_autodoc_typehints",
    "myst_parser",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

root_doc = "index"

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}
napoleon_google_docstring = True
napoleon_numpy_docstring = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

# -- i18n ------------------------------------------------------------------
# Base language for the docs source is English; Japanese is maintained as a
# gettext translation catalog under locale/ja/LC_MESSAGES via sphinx-intl.
# Build the Japanese site with: sphinx-build -D language=ja ...

language = "en"
locale_dirs = ["locale/"]
gettext_compact = False
gettext_uuid = True

# -- HTML output -------------------------------------------------------

html_theme = "furo"
html_static_path = ["_static"]
html_title = f"{project} {release}"

# -- LaTeX / PDF output -------------------------------------------------
# `make latexpdf` (English, pdflatex) and
# `sphinx-build -b latex -D language=ja . _build/latex/ja` (Japanese, uplatex)
# require a local TeX installation (e.g. `brew install --cask mactex-no-gui`,
# or texlive-lang-japanese on Linux, for the Japanese build).

latex_engine = "uplatex" if language == "ja" else "pdflatex"
latex_elements = {
    "papersize": "a4paper",
}
latex_documents = [
    (root_doc, "ProteinDF_bridge.tex", "ProteinDF\\_bridge Documentation", author, "manual"),
]
