"""Configuration for the public CHRONOS documentation."""

import atexit
import re
from pathlib import Path
from shutil import copyfile
from tempfile import TemporaryDirectory
from textwrap import dedent

project = "CHRONOS"
author = "Ariel Medero Borrell"
copyright = "2026, Ariel Medero Borrell"
language = "en"

extensions = [
    "myst_parser",
    "sphinx.ext.mathjax",
    "sphinx.ext.autodoc",
    "sphinxcontrib.matlab",
]
source_suffix = {".md": "markdown"}
root_doc = "index"
myst_enable_extensions = ["dollarmath", "amsmath"]
highlight_language = "matlab"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Give the recursive MATLAB parser only the selected public .m files.
# Never scan the repository's examples, live scripts, or autosave files.
_repository_root = Path(__file__).resolve().parent.parent
_matlab_sources = TemporaryDirectory(prefix="chronos-api-")
atexit.register(_matlab_sources.cleanup)
_api_source_paths = (
    "mpc_definition_api/init_func/init_mpc.m",
    "mpc_definition_api/init_func/init_mpc_dynamics.m",
    "mpc_definition_api/init_func/init_mpc_output.m",
    "mpc_definition_api/init_func/init_mpc_Tracking_cost.m",
    "mpc_definition_api/init_func/init_mpc_Control_cost.m",
    "mpc_definition_api/init_func/init_mpc_ControlRate_cost.m",
    "mpc_definition_api/init_func/init_mpc_Custom_cost.m",
    "mpc_definition_api/init_func/init_mpc_state_cnstr.m",
    "mpc_definition_api/init_func/init_mpc_control_cnstr.m",
    "mpc_definition_api/init_func/init_mpc_control_rate_cnstr.m",
    "mpc_definition_api/init_func/init_mpc_output_cnstr.m",
    "mpc_definition_api/init_func/init_mpc_custom_cnstr.m",
    "mpc_definition_api/init_func/init_mpc_ter_ingredients.m",
    "mpc_definition_api/init_func/init_mpc_ter_ingredients_dlqr.m",
    "mpc_definition_api/init_func/build_chronos_mpc.m",
    "mpc_definition_api/update_func/update_mpc_dynamics.m",
    "mpc_definition_api/update_func/update_mpc_output_vector.m",
    "mpc_definition_api/update_func/update_mpc_Tracking_cost.m",
    "mpc_definition_api/update_func/update_mpc_Control_cost.m",
    "mpc_definition_api/update_func/update_mpc_Control_Rate_cost.m",
    "mpc_definition_api/update_func/update_mpc_Custom_cost_vector.m",
    "mpc_definition_api/update_func/update_mpc_Custom_cost.m",
    "mpc_definition_api/update_func/update_mpc_state_cnstr_limits.m",
    "mpc_definition_api/update_func/update_mpc_control_cnstr_limits.m",
    "mpc_definition_api/update_func/update_mpc_control_rate_cnstr_limits.m",
    "mpc_definition_api/update_func/update_mpc_output_cnstr_limits.m",
    "mpc_definition_api/update_func/update_mpc_custom_cnstr_vector.m",
    "mpc_definition_api/update_func/update_mpc_custom_cnstr_limits.m",
    "mpc_definition_api/update_func/update_mpc_slack_cost.m",
    "mpc_definition_api/update_func/update_mpc_ter_ingredients.m",
    "solver/mpc_solve.m",
)
_api_function_names = {
    str(Path(path).with_suffix("")).replace("/", ".") for path in _api_source_paths
}
for _relative_path in _api_source_paths:
    _destination = Path(_matlab_sources.name) / _relative_path
    _destination.parent.mkdir(parents=True, exist_ok=True)
    copyfile(_repository_root / _relative_path, _destination)

matlab_src_dir = _matlab_sources.name
matlab_short_links = True

html_theme = "sphinx_rtd_theme"
html_title = "CHRONOS documentation"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_context = {
    "display_github": True,
    "github_user": "arielmb94",
    "github_repo": "CHRONOS-MPC",
    "github_version": "docs",
    "conf_py_path": "/docs/",
}


def _format_matlab_help(app, what, name, obj, options, lines):
    """Format existing MATLAB help as prose, argument entries, and examples."""
    # MATLAB help indentation is not reStructuredText: passing it through
    # can produce indentation warnings and unformatted examples. Adapt the
    # layout without changing descriptions or treating later prose as code.
    if name not in _api_function_names:
        return
    formatted = []
    in_arguments = False
    index = 0
    while index < len(lines):
        line = lines[index]
        text = line.strip()
        if text.startswith("Example"):
            in_arguments = False
            formatted.extend(["", f"**{text}**", ""])
        elif text in {"Input:", "Inputs:", "Output:", "Outputs:"}:
            in_arguments = True
            formatted.extend(["", f"**{text}**", ""])
        elif in_arguments and (argument := re.match(r"^(\w+)\s+-\s+(.*)$", text)):
            formatted.extend(["", f"``{argument[1]}`` — {argument[2]}"])
        elif not in_arguments and text and len(line) - len(line.lstrip()) >= 4:
            block = []
            while index < len(lines):
                item = lines[index]
                if item.strip() and len(item) - len(item.lstrip()) < 4:
                    break
                block.append(item)
                index += 1
            code = dedent("\n".join(block)).strip("\n")
            formatted.extend(["", ".. code-block:: matlab", ""])
            formatted.extend("    " + item for item in code.splitlines())
            formatted.append("")
            continue
        else:
            formatted.append(text)
        index += 1
    lines[:] = formatted + [""]


def setup(app):
    app.connect("autodoc-process-docstring", _format_matlab_help)
