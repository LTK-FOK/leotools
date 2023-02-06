"""Command line interface for the Lechner Earth Observation Toolset.

Examples:

    To list all the available functions, use the `list` argument.

        $ leotools list

    To print out help for a function, use `help` and a function name together.

        $ leotools help func_name
    
    To use a function, pass the function name and the args as positional 
    arguments and the kwargs as optional arguments.

        $ leotools func_name arg1 arg2 -kw3 arg3 -kw4 arg4
"""

import argparse
import re

### Other modules
import leotools.gistools as ltgt
import leotools.preproc as ltpp

# Catalog of outside functions
outside_catalog = {
    'profile': ltgt.print_profile,
    'tile': ltpp.reproj_tile,
    'datatake': ltpp.merge_datatake,
    'preproc': ltpp.preproc,
    'reformat': ltpp.reformat,
    'extras': ltpp.make_extras,
    'stack': ltgt.stack_images,
}

### Inner functions, not included in the catalogue
def arg_type(arg):
    """Casts string argument into correct int, float or string type."""

    if re.fullmatch(r'\d+', arg):
        return int(arg)
    elif re.fullmatch(r'\d*.\d+', arg):
        return float(arg)
    elif arg == 'None':
        return None
    elif arg == 'False':
        return False
    elif arg == 'True':
        return True
    else:
        return arg

def handle_kwargs(kwargs_list):
    """Converts optional command line arguments into a dictionary."""

    range_len = range(0, len(kwargs_list), 2) ### Skip every other number
    keys = [kwargs_list[i].lstrip('-') for i in range_len]
    vals = [arg_type(kwargs_list[i+1]) for i in range_len]
    return dict(zip(keys, vals))

### Helper functions, parts of the inner catalog
def help_func(func_name=None):
    """Prints help for the module or a function if one is given as an argument."""

    print() ### An empty line

    if func_name:
        if func_name in CATALOG:
            print(f"{func_name}\n")
            doc = CATALOG[func_name].__doc__
            print(doc)
        else:
            raise ValueError("Function name not found. Use the `list` command to get the function names.")

    else:
        print(__doc__)

def list_func():
    """Lists available functions."""

    print("\nAvailable functions:")

    for k in inside_catalog:
        doc_line = inside_catalog[k].__doc__.split('\n', 1)[0]
        print(f"{4*' '}{k}{(20-len(k))*' '}{doc_line}") ### Currently the list indent is 4 spaces

    print() ### An empty line

    for k in outside_catalog:
        doc_line = outside_catalog[k].__doc__.split('\n', 1)[0]
        print(f"{4*' '}{k}{(20-len(k))*' '}{doc_line}") ### Currently the list indent is 4 spaces

def info_func():
    """Information about common args.

    Inputs paths can be:
        - Path to a single file matching the pattern.
        - Path to a directory containing matching files.
        - Python list of paths to files and directories.
        - Path to a .txt file containing a list with each entry in a new line.

    Modes used for checking the output paths:
        0 - Raise error if path does not exist.
        1 - Create path and raise error if parents are missing.
        2 - Create path with all its missing parents.
    """

    print(info_func.__doc__)

### Catalog of inside functions
inside_catalog = {
    'help': help_func,
    'list': list_func,
    'info': info_func,
}

### Full functoin catalog
CATALOG = {**inside_catalog, **outside_catalog}

def main():
    """Calls the function and passes it the given arguments."""

    ### Argaparse
    parser = argparse.ArgumentParser(description="Use the `help` func to know more") ### --help or -h will print the description
    parser.add_argument('func', help="The name of the function you intend to run from the commandline.")
    parser.add_argument('args', nargs='*', default=[], help="Positional arguments passed to the function.")
    raw_args = parser.parse_known_args()
    
    args = [arg_type(i) for i in raw_args[0].args]
    kwargs = handle_kwargs(raw_args[1])
    CATALOG[raw_args[0].func](*args, **kwargs)

if __name__ == "__main__":
    main()