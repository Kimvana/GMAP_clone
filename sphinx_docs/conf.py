import sys
from pathlib import Path

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'GMAP'
copyright = '2023, KE van Adrichem'
author = 'KE van Adrichem'
release = '0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc', 'numpydoc', 'sphinx_design',
              'sphinx.ext.mathjax']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

sys.path.append(str(Path(__file__).parent.parent.resolve()))

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'pydata_sphinx_theme'
html_theme_options = {
    "logo": {
        "image_light": "Figures/GMAP_banner_light.png",
        "image_dark": "Figures/GMAP_banner_dark.png",
    },
    "github_url": "https://github.com/lacourjansenlab/GMAP",
    "navbar_end": [
        "search-button",
        "theme-switcher",
        "navbar-icon-links",
    ],
    "navbar_persistent": [],
}
# html_static_path = ['_static']
html_favicon = "Figures/GMAP_ico.ico"

suppress_warnings = ["docutils"]

# -- Other options -----------------------------------------------------------
numpydoc_class_members_toctree = False
