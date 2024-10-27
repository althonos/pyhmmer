# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Imports -----------------------------------------------------------------

import datetime
import os
import sys
import re
import shutil
import semantic_version
import urllib.request

# -- Path setup --------------------------------------------------------------

docssrc_dir = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.dirname(docssrc_dir)

# When building on ReadTheDocs, we can't provide a local version of the Cython
# extensions, so we have to install the latest public version, and avoid
# patching the PYTHONPATH with the local development folder
if os.getenv("READTHEDOCS", "False") != "True":
    sys.path.insert(0, project_dir)

# Download the *See Also* cards from a centralized location so it can be kept
# up-to-date across all projects
with urllib.request.urlopen("https://gist.githubusercontent.com/althonos/5d6bf5a512d64dc951c42a91d5fc3fb3/raw/related.rst") as src:
    with open("related.rst", "wb") as dst:
        shutil.copyfileobj(src, dst)


# -- Project information -----------------------------------------------------

import pyhmmer

# extract the project metadata from the module itself
project = pyhmmer.__name__
author = re.match('(.*) <.*>', pyhmmer.__author__).group(1)
year = datetime.date.today().year
copyright = '{}, {}'.format("2020" if year==2020 else "2020-{}".format(year), author)

# extract the semantic version
semver = semantic_version.Version.coerce(pyhmmer.__version__)
version = str(semver.truncate(level="patch"))
release = str(semver)

# patch the docstring of `pyhmmer` so that we don't show the link to redirect
# to the docs (we don't want to see it when reading the docs already, duh!)
doc_lines = pyhmmer.__doc__.splitlines()
if "See Also:" in doc_lines:
    see_also = doc_lines.index("See Also:")
    pyhmmer.__doc__ = "\n".join(doc_lines[:see_also])

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.todo",
    "sphinx.ext.extlinks",
    "sphinx_design",
    "sphinxcontrib.jquery",
    "nbsphinx",
    "recommonmark",
    "IPython.sphinxext.ipython_console_highlighting",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    '_build',
    'Thumbs.db',
    '.DS_Store',
    '**.ipynb_checkpoints',
    'requirements.txt'
]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "monokailight"

# The name of the default role for inline references
default_role = "py:obj"


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'pydata_sphinx_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static/js', '_static/bibtex', '_static/css']
html_js_files = ["custom-icon.js"]
html_css_files = ["custom.css"]

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
html_theme_options = {
    "external_links": [
        {
            "url": "https://doi.org/10.1093/bioinformatics/btad214",
            "name": "Paper",
        },
    ],
    "show_toc_level": 2,
    "use_edit_page_button": True,
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/althonos/pyhmmer",
            "icon": "fa-brands fa-github",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/pyhmmer",
            "icon": "fa-custom fa-pypi",
        },
    ],
    "logo": {
        "text": "PyHMMER",
        "image_light": "_images/logo.png",
        "image_dark": "_images/logo.png",
    },
    "navbar_align": "left",
    "footer_start": ["copyright"],
    "footer_center": ["sphinx-version"],
}

html_context = {
    "github_user": "althonos",
    "github_repo": "pyhmmer",
    "github_version": "master",
    "doc_path": "docs",
}

html_favicon = '_images/favicon.ico'

# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = pyhmmer.__name__


# -- Extension configuration -------------------------------------------------

# -- Options for extlinks extension ------------------------------------------

extlinks = {
    'doi': ('https://doi.org/%s', 'doi:%s'),
    'pmid': ('https://pubmed.ncbi.nlm.nih.gov/%s', 'PMID:%s'),
    'isbn': ('https://www.worldcat.org/isbn/%s', 'ISBN:%s'),
}

# -- Options for imgmath extension -------------------------------------------

imgmath_image_format = "svg"

# -- Options for napoleon extension ------------------------------------------

napoleon_include_init_with_doc = True
napoleon_include_special_with_doc = True
napoleon_include_private_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = True
napoleon_use_rtype = False

# -- Options for autodoc extension -------------------------------------------

autoclass_content = "class"
autodoc_member_order = 'groupwise'
autosummary_generate = []
autodoc_typehints = 'none'

# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "psutil": ("https://psutil.readthedocs.io/en/latest/", None),
    "biopython": ("https://biopython.org/docs/latest/", None),
    "numpy": ("https://numpy.org/doc/stable/", None)
}

# -- Options for recommonmark extension --------------------------------------

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}

# -- Options for nbsphinx extension ------------------------------------------

nbsphinx_execute = 'auto'
nbsphinx_execute_arguments = [
    "--InlineBackend.figure_formats={'svg', 'pdf'}",
    "--InlineBackend.rc={'figure.dpi': 96}",
]
