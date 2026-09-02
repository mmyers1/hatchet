import argparse
from importlib.resources import path
from pathlib import Path
import webbrowser

import hatchet.data


def parsing_arguments(args=None):
    parser = argparse.ArgumentParser(
        prog='hatchet fixed-matrix-editor',
        description='Open the offline editor for HATCHet fixed matrix entry files.',
    )
    parser.add_argument(
        '--no-open',
        action='store_true',
        help='Print the editor path without opening a browser.',
    )
    return parser.parse_args(args)


def main(args=None):
    parsed = parsing_arguments(args)
    with path(hatchet.data, 'fixed_matrix_editor.html') as resource_path:
        editor_path = Path(resource_path).resolve()
        print(editor_path)
        if not parsed.no_open:
            webbrowser.open(editor_path.as_uri(), new=2)
        return editor_path
