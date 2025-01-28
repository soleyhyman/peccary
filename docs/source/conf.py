# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Needed imports ----------------------------------------------------------
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here.
# import pathlib
# import sys
import pathlib, sys
# sys.path.insert(0, pathlib.Path('source/conf.py').resolve().parents[2].as_posix()+'/peccary')
sys.path.insert(0, pathlib.Path(__file__).resolve().parents[2].as_posix()+'/peccary')

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'PECCARY'
copyright = '2024, Sóley Hyman'
author = 'Sóley Hyman'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.viewcode',
    'numpydoc',  # needs to be loaded *after* autodoc
    # 'sphinx.ext.napoleon',
    # 'autoapi.extension',
    'sphinx_automodapi.automodapi',
    'sphinx_automodapi.smart_resolver'
]

templates_path = ['_templates']
exclude_patterns = []

# autoapi_type = 'python'
# autoapi_dirs = ['../../peccary/']
# autoapi_dirs = ['C:/Users/soley/Documents/Arizona/AAA_GalaxyLab/peccary/peccary']
numpydoc_show_class_members = False
autodoc_typehints = "none"


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
# html_theme = 'furo'
html_theme = "pydata_sphinx_theme"
html_static_path = ['_static']
html_logo = "_static/peccary-logo-banner.png"
html_favicon = "_static/peccary-logo-icon.ico"
html_theme_options = {
  "show_nav_level": 2,
  "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/soleyhyman/peccary",
            "icon": "fa-brands fa-square-github",
            "type": "fontawesome",
        },
    ],
}
autosummary_generate = True