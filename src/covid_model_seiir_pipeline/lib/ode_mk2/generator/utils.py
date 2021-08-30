import itertools
from pathlib import Path
from typing import Iterable, List, Union

import inflection

TEXTWIDTH = 118  # type: int
TAB = '    '  # type: str
SPACING = '\n\n'  # type: str


def make_module_docstring(description: str, file: Union[str, Path]) -> str:
    """Generates standard header with additional information from the description.

    Parameters
    ----------
    description
        Custom header text for the generated file.
    file
        Path to python module that generates the target module.

    Returns
    -------
        String representation of module doc string.

    """
    here = Path(file).resolve()
    out = f'"""{description}\n\n'
    out += f'This code is automatically generated by {Path(*here.parts[-2:])}\n\n'
    out += 'Any manual changes will be lost.\n"""\n'
    return out


def make_import(module_to_import: str, imports: Iterable[str] = ()) -> str:
    """Generates the necessary imports. Smart about importing modules or names.

    Parameters
    ----------
    module_to_import
        Name of the module.
    imports
        Named items to import. If empty import the module name.

    Returns
    -------
        Generated string for necessary imports.

    """
    if not imports:
        out = f"import {module_to_import}\n"
    else:
        out = hard_text_wrap(f'from {module_to_import} import (\n', imports, suffix=')\n')

    return out


def hard_text_wrap(start_string, items, suffix: str = ''):
    out = start_string
    for item in items:
        out += TAB + item + ',\n'
    out += suffix
    return out


def format_block(var1: Iterable[str], var2: Iterable[str],
                 indent_level: int = 1, max_len_var1: int = None) -> str:
    max_len_var1 = max([len(v) for v in var1]) if max_len_var1 is None else max_len_var1
    out = ''
    for v1 in var1:

        extra_space = ' ' * (max_len_var1 - len(v1) + 1)
        out += TAB * indent_level
        for v2 in var2:
            out += f"'{v1}_{v2}',{extra_space}"
        out += '\n'
    return out


def format_lines(content: Union[List[str], List[List[str]]], indent_level: int = 1) -> str:
    indent = TAB * indent_level
    if not content:
        return indent + "''\n"
    elif isinstance(content[0], str):
        return indent + ', '.join([f"'{c}'" for c in content]) + ',\n'
    else:
        out = ''
        col_widths = get_longest_item_by_column(content)
        for row in content:
            out += indent
            if len(row) == 1 and row[0] == '\n':
                out += '\n'
                continue
            for i, item in enumerate(row):
                extra_space = ' ' * (col_widths[i] - len(item) + 1)
                out += f"'{item}',{extra_space}"
            out += '\n'
        return out


def get_longest_item_by_column(content: List[List[str]]):
    lengths = []
    for row in content:
        for i, item in enumerate(row):
            if len(lengths) < i + 1:
                lengths.append(len(item))
            else:
                lengths[i] = max(lengths[i], len(item))
    return lengths


def make_named_tuple(name: str, content: Union[List[str], List[List[str]]], private: bool =True):
    prefix = '_' if private else ''
    return f"{prefix}{name} = namedtuple('{name}', [\n{format_lines(content)}])\n"


def make_content_array(rows: Union[List[str], List[List[str]]],
                       columns: Union[List[str], List[List[str]]],
                       row_prefix: bool = True) -> List[List[str]]:
    if not rows or not columns:
        raise ValueError
    if isinstance(rows[0], list):
        rows = ['_'.join(item) for item in itertools.product(*rows)]
    if isinstance(columns[0], list):
        columns = ['_'.join(item) for item in itertools.product(*columns)]

    template = '{row}_{column}' if row_prefix else '{column}_{row}'
    return [
        [template.format(row=r, column=c).rstrip('_') for c in columns] for r in rows
    ]


def make_index_map(map_name: str):
    return f'{inflection.underscore(map_name).upper()} = _{map_name}(*list(range(len(_{map_name}._fields))))\n'


def make_name_map(map_name: str, suffix: str = '_NAMES'):
    return f'{inflection.underscore(map_name).upper()}{suffix} = _{map_name}(*_{map_name}._fields)\n'



