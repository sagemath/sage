r"""
Implementation of the command ``sage --fiximports``.

This file provides a tool to fix the modularization antipattern ``namespace_pkg_all_import``
reported by ``tox -e relint``. Sage library code should not import from ``sage.PAC.KAGE.all``
when ``sage.PAC.KAGE`` is an implicit namespace package.

AUTHORS:

- Alex Chandler
- Matthias Köppe

INSTRUCTIONS:

To fix the above issue for all Python and Cython source files (``.py``, ``.pyx``, ``.pxi``) files in the ``src/sage`` directory,
run the following from a terminal in ``SAGE_ROOT`` ::

    ./sage -python src/sage/misc/replace_dot_all.py

or ::

    ./sage --fiximports

Running replace_dot_all.py will call the function :func:`walkdir_replace_dot_all` which walks through all Python and Cython source files
and replaces certain ``from sage.PAC.KAGE.all import something`` (those matching the pattern above)
with the correct ``import`` statement by applying the function :func:`~sage.misc.dev_tools.import_statements`.

The user can also pass subdirectories of ``src/sage`` or specific files to fix. For example ::

    ./sage -python src/sage/misc/replace_dot_all.py src/sage/arith

will fix all files in ``src/sage/arith`` and ::

    ./sage -python src/sage/misc/replace_dot_all.py src/sage/arith/functions.pyx

will fix the file ``src/sage/arith/functions.py``. The file extension is necessary in the case of a specific file. The user can also
pass the verbose flag ``-v`` to print out the files being fixed. For example ::

    ./sage -python src/sage/misc/replace_dot_all.py -v src/sage/arith

will fix all files in ``src/sage/arith`` and print out the unusual examples of ``import`` statements it finds.

In some rare cases, such as ``import`` statements appearing in doctests, the program will not be able to fix the ``import`` statement. The program will
print out the location of the file and the line number of the exceptional ``import`` statement. The user can then manually fix the ``import`` statement.
The program will also (usually) print out the suggested replacement for the ``import`` statement. The user can then copy and paste this replacement
into the file. In the cases a suggested replacement is not printed out, the user should use the function :func:`~sage.misc.dev_tools.import_statements`
to find the correct ``import`` statement.
"""

# ****************************************************************************
#       Copyright (C) 2022-2023 Alex Chandler
#                     2023      Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

# Importing packages

from sage.misc.dev_tools import import_statements
import os
import re
import sys
import argparse

# We import this using __import__ so that "tox -e relint" does not complain about this source file.
__import__("sage.all", globals(), locals(), ["*"])


# Keep in sync with SAGE_ROOT/src/.relint.yml (namespace_pkg_all_import)

default_package_regex = (r"sage("
                         r"|[.](arith|categories|combinat|crypto|databases|data_structures|dynamics|ext|game_theory|games|geometry|graphs|groups|interfaces|manifolds|matrix|matroids|misc|modules|monoids|numerical|probability|quadratic_forms|quivers|rings|sat|schemes|sets|stats|tensor)[a-z0-9_.]*|[.]libs"
                         r")[.]all")


# Global variables

examples = list('ABCDEFGHIJ')  # controls how we print out interesting examples to the console
interesting_examples = dict(zip(examples, [0]*len(examples)))
log_messages = ''
number_examples_to_print = 100  # controls how many examples we print out to the console (100 effectively allows all unusual examples to be printed)
numberFiles, numberFilesMatchingRegex, numberFilesChanged, numberStatementsReplaced = 0, 0, 0, 0  # to print report on number of files changed


# Functions

def find_replacements(location, package_regex=None, verbose=False):
    r"""
    Locate the lines in the file at ``location`` which contain an ``import`` statement.

    INPUT:

    - ``location`` -- a file path
    - ``package_regex`` -- (default: :obj:`default_package_regex`) a regular expression matching
      the ``sage.PAC.KAGE.all`` package names from which we do not want to import.
    - ``verbose`` -- a parameter which if used will issue print statements when interesting examples are found

    OUTPUT:

    an array [row_index,import_index,replaced_commands,lines_spanned (optional)] with entries

    - ``row_index`` -- the row index (zero indexed) in the file which needs replacing
    - ``import_index`` -- the index of row row_index which marks the beginning of the word ``import`` in the line
    - ``replaced_commands`` -- the string which will replace the current line
    - ``lines_spanned`` -- the number of lines the original statement spans

    This output can be processed by the function :func:`process_line`.

    EXAMPLES::

        sage: # needs SAGE_SRC
        sage: from sage.misc.replace_dot_all import *
        sage: location = os.path.join(sage.env.SAGE_SRC, 'sage', 'plot', 'arc.py')
        sage: find_replacements(location, package_regex='sage[.]plot[.]all', verbose=True)
        [[..., ..., 'from sage.plot.graphics import Graphics']]
    """
    if package_regex is None:
        package_regex = default_package_regex
    regex = r"from\s+" + package_regex + r"\s+import"
    pattern = re.compile(regex)
    replacements = []
    global log_messages, interesting_examples
    with open(location) as fp:
        skip_line = False
        lines = fp.readlines()  # read all lines using readline()
        row_index = 0
        for row in lines:  # iterate over all lines of python file
            if pattern.search(row):  # (match the regex also do not want to mess with documentation)
                prefix = ''
                if '*' in row or 'SAGE_ROOT' in row:
                    if verbose and interesting_examples['J'] < number_examples_to_print:
                        interesting_examples['J'] += 1
                        log_messages += (f'J. Match but no changes made (import statement uses *) at {location}:{row_index + 1}. '
                                         f'Not applying any changes here.\n')
                    continue
                elif not (row.lstrip()[0:4] == 'from'):
                    skip_line = True
                    if '"' not in row and "'" not in row:
                        print(f'\n'
                              f'NEED TO CHANGE MANUALLY \n'
                              f'  Issue: line with import statement does not start with "from" \n'
                              f'  Location: at {location} \n'
                              f'  Line number: {row_index + 1}. \n'
                              f'  Giving correct import statements:\n')
                        leading_space = 0
                        while len(row) > 0 and row[leading_space] == ' ' and leading_space < len(row)-1:
                            leading_space += 1
                        prefix_space = leading_space
                        while row[prefix_space:prefix_space+4] != 'from':
                            prefix_space += 1
                        prefix = row[leading_space:prefix_space]
                        row = row[prefix_space:]
                    else:
                        print(f'\n'
                              f'NEED TO CHANGE MANUALLY \n'
                              f'  Issue: import statement does not start with "from" and contains quotation marks \n'
                              f'  Location: at {location}:{row_index + 1}. \n'
                              f'  Not able to suggest correct import statements. User must use the function import_statements().')
                        continue
                # find() method returns -1 if the value is not found, or if found it returns index of the first occurrence of the substring
                import_index = row.find('import ')
                modules = ''
                to_exec = row.strip()  # the import statement itself which we will clean and call to import modules e.g. "import (aa as a, bb, cc as c)"
                to_eval = row[import_index + 7:-1].strip()  # tuple of modules in import statement but e.g. "(aa as a, bb, cc as c)" is saved as "(a,bb,c)"
                to_eval_raw = row[import_index + 7:-1].strip()  # same as to_eval but we don't get rid of " as " parts e.g. "(aa as a, bb, cc as c)"
                span = 0  # keeps track of how many lines the import statement spans

                if '(' in row:  # for examples where we import a tuple of modules and the statement spans several lines
                    while ')' not in lines[row_index + span]:  # finding the line which closes the import statement
                        span += 1
                        to_exec += lines[row_index + span].strip()
                        to_eval += lines[row_index + span].strip()
                        to_eval_raw += lines[row_index + span].strip()
                    if span and verbose:  # useful to see these multiline examples for debugging
                        if " as " in to_eval_raw and interesting_examples['D'] < number_examples_to_print:
                            log_messages += f'D. Interesting example (spans multiple lines and has " as ") at {location}:{row_index + 1}\n'
                            interesting_examples['D'] += 1
                        elif interesting_examples['B'] < number_examples_to_print:
                            log_messages += f'B. Interesting example (spans multiple lines) at {location}:{row_index + 1}\n'
                            interesting_examples['B'] += 1

                # if there is an "as" statement inside to_eval, we want to keep only the new name for the module e.g. "(aa as a, bb, cc as c)" becomes "(a,bb,c)"
                while " as " in to_eval:
                    as_ind = to_eval.find(" as ")
                    j = as_ind - 1
                    while to_eval[j] not in [',', ' '] and j >= 0:
                        j -= 1
                    to_eval = to_eval[0:j+1] + to_eval[as_ind+4:]

                try:  # trying to execute the import statement so we can eval the modules and feed them to the function import_statements
                    to_exec = to_exec.replace("'", '').replace('"', '')
                    if (to_exec[-1] == ','):
                        to_exec = to_exec[:-1]
                    if sys.version_info.minor < 13:
                        exec(to_exec)
                    else:
                        loc = locals()
                        exec(to_exec, locals=loc)
                except ModuleNotFoundError as err:
                    print(f'ModuleNotFoundError: {err} found when trying to execute {to_exec}')
                except ImportError as err:
                    print(f'ImportError: {err} found when trying to execute {to_exec}')

                try:  # try to evaluate the list of module names to get a list of the modules themselves which we can call import_statements on
                    if sys.version_info.minor < 13:
                        modules = eval(to_eval)
                    else:
                        modules = eval(to_eval, locals=loc)
                except NameError as err:
                    print(f'NameError: {err} found when trying to evaluate {to_eval} at {location}:{row_index + 1}')
                except SyntaxError as err:
                    print(f'SyntaxError: {err} found when trying to evaluate {to_eval} at {location}:{row_index + 1}')

                # Need module to be a list of modules we are importing. If a single module was given, we make it a 1-element list.
                if not isinstance(modules, tuple):
                    modules = [modules]

                to_eval_list = to_eval.replace('(', '').replace(')', '').split(',')  # convert comma separated string to_eval to a list
                to_eval_list_raw = to_eval_raw.replace('(', '').replace(')', '').split(',')  # convert comma separated string to_eval_raw to a list
                to_eval_list_index = 0

                change_to = ''
                # constructs the appropriate replacement for the import statement and stores it (along with location data) in list replacements
                for mod in modules:
                    postfix = ''
                    as_index = -1
                    # saves the callable name of module in variable postfix (e.g. for module "b" called by "bb as b" we set postfix = " as b")
                    if " as " in to_eval_list_raw[to_eval_list_index]:
                        if verbose and interesting_examples['C'] < number_examples_to_print:
                            log_messages += f'C. Interesting example (" as " in tuple import) at {location}:{row_index + 1}\n'
                            interesting_examples['C'] += 1
                        as_index = to_eval_list_raw[to_eval_list_index].index(" as ")
                        postfix = to_eval_list_raw[to_eval_list_index][as_index:]
                    new_import_statement = import_statements(mod, answer_as_str=True, verbose=False)  # import statement for the current mod in the list module
                    import_index = new_import_statement.find('import')
                    new_mod_as_string = new_import_statement[import_index + 7:].strip()  # the name for the module given by the function import_statements
                    if as_index >= 0:
                        # the name for the module as originally called in the document (when there is an " as " statement)
                        original_mod_string = to_eval_list_raw[to_eval_list_index].strip()[:as_index]
                    else:
                        original_mod_string = to_eval_list[to_eval_list_index].strip()  # the name for the module as originally called in the document
                    if original_mod_string != new_mod_as_string:  # if the names differ, we use the original name as it was called in the document
                        if verbose and interesting_examples['A'] < number_examples_to_print:
                            log_messages += (f'A. Interesting example (module has multiple names) at {location}:{row_index + 1}. '
                                             f'Names: {original_mod_string}, {new_mod_as_string}. '
                                             f'Replacing new {new_mod_as_string} by original {original_mod_string}.\n')
                            interesting_examples['A'] += 1
                        new_import_statement = new_import_statement.replace(' ' + new_mod_as_string, ' ' + new_mod_as_string + ' as ' + original_mod_string)
                        if " as " in postfix and interesting_examples['G'] < number_examples_to_print:
                            log_messages += (f'G. Interesting example (module has multiple names) at {location}:{row_index + 1}. '
                                             f'Names: {original_mod_string}, {new_mod_as_string}. '
                                             f'Replacing new {new_mod_as_string} by original {original_mod_string}.\n')
                            interesting_examples['G'] += 1
                    if len(postfix.strip()) > 0:  # if module was called with " as " statement, we put that back in by adding the string "postfix"
                        # if " as " in new_import_statement locate the index of " as ", remove the end after this, and add the postfix there
                        if " as " in new_import_statement:
                            new_import_statement = new_import_statement[:new_import_statement.index(" as ")] + ' ' + postfix.strip()
                        else:
                            new_import_statement += (' ' + postfix.strip())
                    change_to += (prefix + new_import_statement + '\n')
                    to_eval_list_index += 1
                # [:-1] on change_to gets rid of the last '\n' we added which adds an unnecessary new line
                replacement = [row_index, import_index, change_to[:-1]].copy()
                if span:
                    # if original statement spanned multiple lines, we store that information to signal that we need to skip lines
                    # as we read the document in the function make_replacements_in_file
                    replacement.append(span)
                if not skip_line:
                    replacements.append(replacement)
            row_index += 1
            skip_line = False
    # keeping track of the numbers of files changed and statements replaced
    global numberStatementsReplaced, numberFilesChanged
    numberStatementsReplaced += len(replacements)
    if replacements:
        numberFilesChanged += 1
    return replacements


def process_line(location, line, replacements, row_index, verbose=False):
    r"""
    Modify a single source code ``line`` according to the given ``replacements``.

    INPUTS:

    - ``location`` -- a file path; only used for logging
    - ``line`` -- a source code line
    - ``replacements`` -- the array output from :func:`find_replacements`
    - ``row_index`` -- the line number where ``import`` appears
    - ``verbose`` -- if ``True``, issue print statements when interesting
      examples are found

    OUTPUT: an array ``[new_line, replacements]`` with entries

    - ``new_line`` -- the modified import statement (possibly now on several lines)
    - ``replacements`` -- just returns the original replacements with its index 0 element removed if ``replacements`` is nonempty

    EXAMPLES:

    Replacing the first line which needs a replacement in the source file with filepath ``src/sage/plot/arc.py``::

        sage: # needs SAGE_SRC
        sage: from sage.misc.replace_dot_all import *
        sage: location = os.path.join(sage.env.SAGE_SRC, 'sage', 'plot', 'arc.py')
        sage: replacements = find_replacements(location, package_regex='sage[.]plot[.]all', verbose=True); replacements
        [[477, 24, 'from sage.plot.graphics import Graphics']]
        sage: with open(location, "r") as file:
        ....:     lines = file.readlines()
        sage: row_index, col_number, *_ = replacements[0]
        sage: line = lines[row_index]
        sage: print(line.rstrip())
        from sage.plot.all import Graphics
        sage: new_line, replacements = process_line(location, line, replacements, row_index)
        sage: print(new_line)
        from sage.plot.graphics import Graphics
        sage: replacements
        []
    """
    line = line.rstrip()  # stripping line break
    new_line = ''
    global log_messages, interesting_examples
    if len(replacements) == 0:
        return line, replacements
    if row_index == replacements[0][0]:  # if line marked as containing .all
        replacement = replacements.pop(0)
        leading_space = 0
        while line and line[leading_space] == ' ' and leading_space < len(line)-1:
            leading_space += 1
        new_line = ' '*leading_space + replacement[2]  # adds leading space to first line (which may or may not start with 'from')
        new_line = new_line.replace('\n', '\n'+' '*leading_space)  # adds correct amount of indentation to the replacement at each line
        # new_line = replacement[2].replace('from ',' '*leading_space + 'from ') # adds correct amount of indentation to the replacement at each line
        if verbose and leading_space > 0:
            if len(replacement) == 4 and interesting_examples['F'] < number_examples_to_print:
                log_messages += f'F. Interesting example (has leading space and multiline) at {location}:{replacement[0] + 1}\n'
                interesting_examples['F'] += 1
            elif interesting_examples['E'] < number_examples_to_print:
                log_messages += f'E. Interesting example (has leading space) at {location}:{replacement[0] + 1}\n'
                interesting_examples['E'] += 1

    else:  # if line does not contain .all
        new_line = line
    return new_line, replacements


def make_replacements_in_file(location, package_regex=None, verbose=False, output=None):
    r"""
    Replace ``import`` statements in the file with filepath "location".

    INPUT:

    - ``location`` -- a file path
    - ``package_regex`` -- (default: :obj:`default_package_regex`) a regular expression matching
      the ``sage.PAC.KAGE.all`` package names from which we do not want to import.
    - ``verbose`` -- if ``True``, issue print statements when interesting examples are found
    - ``output`` -- a file path; if ``None``, overwrite the file given by ``location``

    EXAMPLES::

        sage: from sage.misc.replace_dot_all import *
        sage: import tempfile
        sage: with tempfile.TemporaryDirectory() as d:
        ....:     location = os.path.join(d, "input.py")
        ....:     with open(location, "w") as input:
        ....:         _ = input.write("from sage.plot.all import point2d\n")
        ....:         _ = input.write("from sage.plot.line import line\n")
        ....:     make_replacements_in_file(location, 'sage[.]plot[.]all', True)
        ....:     with open(location, "r") as output:
        ....:         for line in output:
        ....:             print(line.strip())
        from sage.plot.point import point2d
        from sage.plot.line import line
    """
    replacements = find_replacements(location, package_regex, verbose)
    with open(location) as file:
        lines = file.readlines()
    replaced_content = ""
    row_index = 0  # keeps track of the line number
    while row_index < len(lines):  # looping through the file
        line = lines[row_index]
        span = 0  # keeps track of number of lines import statement spans
        if replacements and row_index == replacements[0][0] and len(replacements[0]) == 4:
            span = replacements[0][3]  # if import statement spans span lines
        # returns the line if no replacements are needed and returns the processed line otherwise
        new_line, replacements = process_line(location, line, replacements, row_index, verbose=verbose)
        replaced_content += new_line + "\n"  # concatenate the new string and add an end-line break
        row_index += 1 + span
    if output is None:
        output = location
    with open(output, "w") as write_file:  # Open file in write mode
        write_file.write(replaced_content)  # overwriting the old file contents with the new/replaced content


def walkdir_replace_dot_all(dir, file_regex=r'.*[.](py|pyx|pxi)$', package_regex=None, verbose=False, *,
                            excluded_file_regex=r'auto-methods|replace_dot_all'):
    r"""
    Replace ``import`` statements in the files in directory ``dir`` matching the regex pattern ``file_regex``.

    INPUTS:

    - ``dir`` -- a directory path
    - ``file_regex`` -- a regular expression matching the file names to process
    - ``package_regex`` -- (default: :obj:`default_package_regex`) a regular expression matching
      the ``sage.PAC.KAGE.all`` package names from which we do not want to import.
    - ``verbose`` -- if ``True``, print statements when interesting examples are found
    - ``excluded_file_regex`` -- a regular expression matching the file names to exclude

    EXAMPLES::

        sage: # needs SAGE_SRC
        sage: from sage.misc.replace_dot_all import *
        sage: walkdir_replace_dot_all(os.path.join(sage.env.SAGE_SRC, 'sage'))  # not tested
    """
    global numberFiles, numberFilesMatchingRegex
    file_regex = re.compile(file_regex)
    excluded_file_regex = re.compile(excluded_file_regex)
    for root, dirs, files in os.walk(dir, topdown=False):
        for name in files:
            numberFiles += 1
            if file_regex.search(name) and not excluded_file_regex.search(name):
                numberFilesMatchingRegex += 1
                location = os.path.join(root, name)
                make_replacements_in_file(location, package_regex, verbose)


# ******************************************************** EXECUTES MAIN FUNCTION **********************************************************************
# this executes the main function in this file which writes over all import statements matching regex in files in specified location matching fileRegex:
if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser()
    # Optional arguments
    parser.add_argument(
        "location",
        metavar='files or directories',
        nargs='*',
        help=("Names of source directories or source files. "
              "If none given, walks through all files in src/sage."),
        type=str)
    parser.add_argument(
        "-v", "--verbose",
        help="Increase output verbosity. Shows locations of any unusual cases of import statements and the corresponding changes.",
        action="store_true")    # Parse arguments
    args = parser.parse_args()
    verbosity = args.verbose
    # Declare regular expressions
    file_regex = r'.*[.](py|pyx|pxi)$'
    package_regex = None
    # Execute the main function based on the specified location and verbosity
    if not args.location:
        from sage.env import SAGE_SRC

        args.location = [os.path.join(SAGE_SRC, 'sage')]
    try:
        for location in args.location:
            if not (location.endswith('.py') or location.endswith('.pxi')):
                # Assume directory
                walkdir_replace_dot_all(location, file_regex, package_regex, verbose=verbosity)
            else:
                # make replacements in file specified by location argument
                make_replacements_in_file(location, package_regex, verbose=verbosity)
    finally:
        # Print report also when interrupted
        if verbosity:
            log_messages_split = sorted(log_messages.rstrip().split('\n'))
            for i, message in enumerate(log_messages_split, start=1):
                # add index to each line
                print(f'{i}. {message.rstrip()}')
        report = 'REPORT:\n'
        report += f'Number of files checked: {numberFiles}\n'
        report += f'Number of files matching regex: {numberFilesMatchingRegex}\n'
        report += f'Number of files changed: {numberFilesChanged}\n'
        report += f'Number of import statements replaced: {numberStatementsReplaced}'
        print('*'*100 + '\n' + report + '\n' + '*'*100)

# ******************************************************************************************************************************************************
