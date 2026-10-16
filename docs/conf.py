"""Configuration for the public CHRONOS documentation."""

project = "CHRONOS"
author = "Ariel Medero Borrell"
copyright = "2026, Ariel Medero Borrell"
language = "en"

extensions = ["myst_parser", "sphinx.ext.mathjax"]
source_suffix = {".md": "markdown"}
root_doc = "index"
myst_enable_extensions = ["dollarmath", "amsmath"]
highlight_language = "matlab"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
html_title = "CHRONOS documentation"
html_context = {
    "display_github": True,
    "github_user": "arielmb94",
    "github_repo": "CHRONOS-MPC",
    "github_version": "docs",
    "conf_py_path": "/docs/",
}
