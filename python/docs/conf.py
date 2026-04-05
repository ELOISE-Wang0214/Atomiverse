"""Sphinx configuration for Atomiverse Python documentation."""

# ── Project info ──────────────────────────────────────────────────────────────
project = "Atomiverse"
author = "Atomiverse Contributors"
copyright = "2025, Atomiverse Contributors"
release = "0.1.0"

# ── Extensions ────────────────────────────────────────────────────────────────
extensions = [
    "myst_parser",           # Markdown source files
    "sphinx.ext.autodoc",    # API docs from docstrings
    "sphinx.ext.napoleon",   # NumPy/Google docstring style
    "sphinx.ext.viewcode",   # Source links
    "sphinx.ext.intersphinx",
    "sphinx_copybutton",     # Copy button on code blocks
]

# ── Source ────────────────────────────────────────────────────────────────────
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
master_doc = "index"

# ── MyST ──────────────────────────────────────────────────────────────────────
myst_enable_extensions = ["colon_fence", "deflist"]

# ── Autodoc ───────────────────────────────────────────────────────────────────
autodoc_member_order = "bysource"
autodoc_typehints = "description"

# ── Intersphinx ───────────────────────────────────────────────────────────────
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
}

# ── HTML ──────────────────────────────────────────────────────────────────────
html_theme = "furo"
html_title = "Atomiverse"
html_theme_options = {
    "sidebar_hide_name": False,
    "light_css_variables": {
        "color-brand-primary": "#0066CC",
        "color-brand-content": "#0066CC",
    },
}
html_static_path = ["_static"]

# ── Copy-button ───────────────────────────────────────────────────────────────
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True
