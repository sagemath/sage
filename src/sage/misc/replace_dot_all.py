r"""
 This file provides a tool to fix the following issue:

 Modularization anti-patterns (from src/.relint.yml)
- name: 'namespace_pkg_all_import: import from .all of a namespace package'
  hint: |
    Sage library code should not import from sage.PAC.KAGE.all when sage.PAC.KAGE is an implicit
    Hint: namespace package. Type import_statements("SOME_IDENTIFIER") to find a more specific import.
  pattern: 'from\s+sage(|[.](arith|categories|combinat|ext|graphs(|[.]decompositions)|interfaces|libs|matrix|misc|numerical(|[.]backends)|rings|sets))[.]all\s+import'
  filePattern: '.*[.](py|pyx|pxi)$'
  error: false # Make this a warning instead of an error for now

AUTHORS:

- Alex Chandler
- Matthias Köppe

INSTRUCTIONS:

To fix the above issue for all .py, .pyx, and .pxi files in the src/sage directory,
run the following from a terminal in SAGE_ROOT ::

    ./sage -python src/sage/misc/replace_dot_all.py

Running replace_dot_all.py will call the function walkdir_replace_dot_all() which walks through all files in
src/sage matching the filePattern above and replaces certain 'from module.all import something' (those matching the pattern above)
with the correct import statement by applying the function import_statements from src/sage/misc/dev_tools.py.

The user can also pass arguments -l to specify a subdirectory of src/sage or even a specific file to fix. The root location is
automatically set to src/sage so you need only specify the path from there. For example ::

    ./sage -python src/sage/misc/replace_dot_all.py -l arith

will fix all files in src/sage/arith and ::

    ./sage -python src/sage/misc/replace_dot_all.py -l arith/functions.pyx

will fix the file src/sage/arith/functions.py. The file extension is necessary in the case of a specific file. The user can also
pass the verbose flag -v to print out the files being fixed. For example ::

    ./sage -python src/sage/misc/replace_dot_all.py -l arith -v

will fix all files in src/sage/arith and print out the unusual examples of import statements it finds.

In some rare cases, such as import statments appearing in doctests, the program will not be able to fix the import statement. The program will
print out the location of the file and the line number of the exceptional import statement. The user can then manually fix the import statement.
The program will also (usually) print out the suggested replacement for the import statement. The user can then copy and paste this replacement
into the file. In the cases a suggested replacement is not printed out, the user should use the function import_statements()
from src/sage/misc/dev_tools.py to find the correct import statement.
"""

# Importing packages

from sage.misc.dev_tools import import_statements
# import sage.all
from sage.all import *
import os
import re
import argparse


# to parse arguments passed to the script

def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser()
    # Optional arguments
    parser.add_argument(
        "-l", "--location", help="Location of directory or file (root set at src/sage so input path from here). If no argument given, walks through all files in src/sage.", type=str)
    parser.add_argument("-v", "--verbose", help="Increase output verbosity. Shows locations of any unusual cases of import statements and the corresponding changes.",
                        action="store_true")    # Parse arguments
    args = parser.parse_args()
    return args


# Global variables
optional_arguments = sys.argv
examples = list('ABCDEFGHIJ')  # controls how we print out interesting examples to the console
interesting_examples = dict(zip(examples, [0]*len(examples)))
log_messages = ''
number_examples_to_print = 100  # controls how many examples we print out to the console (100 effectively allows all unusual examples to be printed)
numberFiles, numberFilesMatchingRegex, numberFilesChanged, numberStatementsReplaced = 0, 0, 0, 0  # to print report on number of files changed


# Functions


def find_replacements(location, regex, verbose=False):
    r"""
    Locates the lines in the file at location which match the regex pattern.

    INPUT:

    - ``location`` -- a file path   file_to_change = 'schemes/elliptic_curves/ell_rational_field.py'
                                    location = cwd + file_to_change
    - ``regex`` -- a regular expression locating strings containing certain module.all import statements. The suggested value is ``regex = r"from\s+sage(|[.](arith|categories|combinat|ext|graphs(|[.]decompositions)|interfaces|libs|matrix|misc|numerical(|[.]backends)|rings|sets))[.]all\s+import"``.
    - ``verbose`` -- a parameter which if used will issue print statements when interesting examples are found

    OUTPUT:

    an array [row_index,import_index,replaced_commands,lines_spanned (optional)] with entries

    - ``row_index`` -- the row index (zero indexed) in the file which needs replacing
    - ``import_index`` -- the index of row row_index which marks the beginning of the word 'import' in the line
    - ``replaced_commands`` -- the string which will replace the current line
    - ``lines_spanned`` -- the number of lines the original statement spans

    this output is specifically designed to be fed into the function process_line(location,line,replacements,import_index) and
    called inside of the function "make_replacements_in_file(location)" to populate the variable "replacements"

    EXAMPLES::

        sage: from sage.misc.replace_dot_all import *
        sage: os.chdir(sage.env.SAGE_SRC + '/sage') # change to sage directory
        sage: cwd = os.getcwd() # Get the current working directory
        sage: file_to_change = 'structure/element.pyx'
        sage: location = cwd + file_to_change
        sage: regex = r"from\s+sage(|[.](arith|categories|combinat|ext|graphs(|[.]decompositions)|interfaces|libs|matrix|misc|numerical(|[.]backends)|rings|sets))[.]all\s+import"
        sage: replacements = find_replacements(location, regex, verbose=True)
    """

    pattern = re.compile(regex)
    replacements = []
    global log_messages, interesting_examples
    with open(location, "r+") as fp:
        skip_line = False
        lines = fp.readlines()  # read all lines using readline()
        row_index = 0
        for row in lines:  # iterate over all lines of python file
            if pattern.search(row):  # (match the regex also do not want to mess with documentation)
                prefix = ''
                if '*' in row or 'SAGE_ROOT' in row:
                    if verbose and interesting_examples['J'] < number_examples_to_print:
                        interesting_examples['J'] += 1
                        log_messages += f'J. Match but no changes made (import statement uses *) at {location} line number {row_index + 1}. Not applying any changes here.\n'
                    continue
                elif not (row.lstrip()[0:4] == 'from'):
                    skip_line = True
                    if '"' not in row and "'" not in row:
                        print(
                            f'\nNEED TO CHANGE MANUALLY \n  Issue: line with import statement does not start with "from" \n  Location: at {location} \n  Line number: {row_index + 1}. \n  Giving correct import statements:\n')
                        leading_space = 0
                        while len(row) > 0 and row[leading_space] == ' ' and leading_space < len(row)-1:
                            leading_space += 1
                        prefix_space = leading_space
                        while row[prefix_space:prefix_space+4] != 'from':
                            prefix_space += 1
                        prefix = row[leading_space:prefix_space]
                        row = row[prefix_space:]
                    else:
                        print(
                            f'\nNEED TO CHANGE MANUALLY \n  Issue: import statement does not start with "from" and contains quotation marks \n  Location: at {location} \n  Line number: {row_index + 1}. \n  Not able to suggest correct import statements. User must use the function import_statements().')
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
                    if span > 0 and verbose:  # useful to see these multiline examples for debugging
                        if " as " in to_eval_raw and interesting_examples['D'] < number_examples_to_print:
                            log_messages += f'D. Interesting example (spans multiple lines and has " as ") at {location} line number {row_index + 1}\n'
                            interesting_examples['D'] += 1
                        elif interesting_examples['B'] < number_examples_to_print:
                            log_messages += f'B. Interesting example (spans multiple lines) at {location} line number {row_index + 1}\n'
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
                    exec(to_exec)
                except ModuleNotFoundError as err:
                    print(f'ModuleNotFoundError: {err} found when trying to execute {to_exec}')
                except ImportError as err:
                    print(f'ImportError: {err} found when trying to execute {to_exec}')

                try:  # try to evaluate the list of module names to get a list of the modules themselves which we can call import_statements on
                    modules = eval(to_eval)
                except NameError as err:
                    print(f'NameError: {err} found when trying to evaluate {to_eval} at {location} line number {row_index + 1}')
                except SyntaxError as err:
                    print(f'SyntaxError: {err} found when trying to evaluate {to_eval} at {location} line number {row_index + 1}')

                # Need module to be a list of modules we are importing. If a single module was given, we make it a 1-element list.
                if not (type(modules) == tuple):
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
                            log_messages += f'C. Interesting example (" as " in tuple import) at {location} at line number {row_index + 1}\n'
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
                            log_messages += f'A. Interesting example (module has multiple names) at {location} line number {row_index + 1}. Names: {original_mod_string}, {new_mod_as_string}. Replacing new {new_mod_as_string} by original {original_mod_string}.\n'
                            interesting_examples['A'] += 1
                        new_import_statement = new_import_statement.replace(' ' + new_mod_as_string, ' ' + new_mod_as_string + ' as ' + original_mod_string)
                        if " as " in postfix and interesting_examples['G'] < number_examples_to_print:
                            log_messages += f'G. Interesting example (module has multiple names) at {location} line number {row_index + 1}. Names: {original_mod_string}, {new_mod_as_string}. Replacing new {new_mod_as_string} by original {original_mod_string}.\n'
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
                if span > 0:  # if original statement spanned multiple lines, we store that information to signal that we need to skip lines as we read the document in the function make_replacements_in_file
                    replacement.append(span)
                if skip_line is False:
                    replacements.append(replacement)
                else:
                    print(replacement[2])
            row_index += 1
            skip_line = False
    # keeping track of the numbers of files changed and statements replaced
    global numberStatementsReplaced, numberFilesChanged
    numberStatementsReplaced += len(replacements)
    if len(replacements) > 0:
        numberFilesChanged += 1
    return replacements


def process_line(location, line, replacements, import_index, verbose=False):
    r"""
    Designed specifically to be called inside of the function make_replacements_in_file(location) to process a single line while iterating over lines in a file

    INPUTS:

    - ``location`` -- a file path   file_to_change = 'schemes/elliptic_curves/ell_rational_field.py'
                                    location = cwd + file_to_change
    - ``line`` -- a line in a file
    - ``replacements`` -- the array output from find_replacements(location)
    - ``import_index`` -- the index in the line which locates 'import'
    - ``verbose`` -- a parameter which if used will issue print statements when interesting examples are found

    OUTPUT:

    an array ``[new_line, replacements]`` with entries

    - ``new_line`` -- the modified import statement (possibly now on several lines)
    - ``replacements`` -- just returns the original replacements with its index 0 element removed if replacements is nonempty

    EXAMPLES:

    - to replace the first line which needs a replacement at the file with filepath #'structure/element.pyx'::

        sage: from sage.misc.replace_dot_all import *
        sage: os.chdir(sage.env.SAGE_SRC + '/sage') # change to sage directory
        sage: cwd = os.getcwd() # Get the current working directory
        sage: file_to_change = 'structure/element.pyx'
        sage: location = cwd + file_to_change
        sage: regex = r"from\s+sage(|[.](arith|categories|combinat|ext|graphs(|[.]decompositions)|interfaces|libs|matrix|misc|numerical(|[.]backends)|rings|sets))[.]all\s+import"
        sage: replacements = find_replacements(location, regex, verbose=True)
        sage: file = open(location, "r")
        sage: lines = file.readlines()
        sage: row_index,import_index = replacements[0],replacements[1]
        sage: line = lines[row_index]
        sage: print(f'old line, old reps: {line,replacements}')
        sage: new_line,replacements = process_line(location,line,replacements,import_index)
        sage: print(f'new line, new reps: {new_line,replacements}')
    """

    line = line.rstrip()  # stripping line break
    new_line = ''
    global log_messages, interesting_examples
    if len(replacements) == 0:
        return line, replacements
    if import_index == replacements[0][0]:  # if line marked as containing .all
        replacement = replacements.pop(0)
        leading_space = 0
        while len(line) > 0 and line[leading_space] == ' ' and leading_space < len(line)-1:
            leading_space += 1
        new_line = ' '*leading_space + replacement[2]  # adds leading space to first line (which may or may not start with 'from')
        new_line = new_line.replace('\n', '\n'+' '*leading_space)  # adds correct amount of indentation to the replacement at each line
        # new_line = replacement[2].replace('from ',' '*leading_space + 'from ') # adds correct amount of indentation to the replacement at each line
        if verbose and leading_space > 0:
            if len(replacement) == 4 and interesting_examples['F'] < number_examples_to_print:
                log_messages += f'F. Interesting example (has leading space and multiline) at {location} row number {replacement[0] + 1}\n'
                interesting_examples['F'] += 1
            elif interesting_examples['E'] < number_examples_to_print:
                log_messages += f'E. Interesting example (has leading space) at {location} row number {replacement[0] + 1}\n'
                interesting_examples['E'] += 1

    else:  # if line does not contain .all
        new_line = line
    return new_line, replacements

# to make all replacements matching the regex in a single file with filepath "location"


def make_replacements_in_file(location, regex, verbose=False):
    r"""
    Writes over the file with filepath "location" making replacements for lines matching the regex pattern: 'from\s+sage(|[.](arith|categories|combinat|ext|graphs(|[.]decompositions)|interfaces|libs|matrix|misc|numerical(|[.]backends)|rings|sets))[.]all\s+import'

    INPUT:

    - ``location`` -- a file path   file_to_change = 'schemes/elliptic_curves/ell_rational_field.py'
                                    location = cwd + file_to_change
    - ``regex`` -- a regular expression locating strings containing certain module.all import statements. The suggested value is ``regex = r"from\s+sage(|[.](arith|categories|combinat|ext|graphs(|[.]decompositions)|interfaces|libs|matrix|misc|numerical(|[.]backends)|rings|sets))[.]all\s+import"``.
    - ``verbose`` -- a parameter which if used will issue print statements when interesting examples are found

    EXAMPLES::
        sage: from sage.misc.replace_dot_all import *
        sage: os.chdir(sage.env.SAGE_SRC + '/sage') # change to sage directory
        sage: cwd = os.getcwd() # Get the current working directory
        sage: file_to_change = 'structure/element.pyx'
        sage: location = cwd + file_to_change
        sage: regex = r"from\s+sage(|[.](arith|categories|combinat|ext|graphs(|[.]decompositions)|interfaces|libs|matrix|misc|numerical(|[.]backends)|rings|sets))[.]all\s+import"
        sage: make_replacements_in_file(location, regex)
    """

    replacements = find_replacements(location, regex, verbose)
    file = open(location, "r")
    replaced_content = ""
    row_index = 0  # keeps track of the line number
    lines = file.readlines()
    while row_index < len(lines):  # looping through the file
        line = lines[row_index]
        span = 0  # keeps track of number of lines import statement spans
        if len(replacements) > 0 and row_index == replacements[0][0] and len(replacements[0]) == 4:
            span = replacements[0][3]  # if import statement spans span lines
        # returns the line if no replacements are needed and returns the processed line otherwise
        new_line, replacements = process_line(location, line, replacements, row_index, verbose=verbose)
        replaced_content += new_line + "\n"  # concatenate the new string and add an end-line break
        row_index += 1 + span
    file.close()  # close the file
    write_file = open(location, "w")  # Open file in write mode
    write_file.write(replaced_content)  # overwriting the old file contents with the new/replaced content
    write_file.close()  # close the file

# to iterate over all files in src/sage matching the given regular expression


def walkdir_replace_dot_all(dir, fileRegex, regex, verbose=False):
    r"""
    Writes over the files in src/sage matching the regex pattern fileRegex making replacements to all lines in such files
    which match the regex pattern: 'from\s+sage(|[.](arith|categories|combinat|ext|graphs(|[.]decompositions)|interfaces|libs|matrix|misc|numerical(|[.]backends)|rings|sets))[.]all\s+import'

    INPUTS:

    - ``fileRegex`` -- a regular expression locating strings containing certain module.all import statements. The suggested value is ``fileRegex = r'.*[.](py|pyx|pxi)$'``.
    - ``regex`` -- a regular expression locating strings containing certain module.all import statements. The suggested value is ``regex = r"from\s+sage(|[.](arith|categories|combinat|ext|graphs(|[.]decompositions)|interfaces|libs|matrix|misc|numerical(|[.]backends)|rings|sets))[.]all\s+import"``.
    - ``verbose`` -- a parameter which if used will issue print statements when interesting examples are found

    EXAMPLES::

        sage: from sage.misc.replace_dot_all import *
        sage: os.chdir(sage.env.SAGE_SRC + '/sage') # change to sage directory
        sage: dir = os.getcwd() # Get the current working directory
        sage: fileRegex = r'.*[.](py|pyx|pxi)$'
        sage: regex = r"from\s+sage(|[.](arith|categories|combinat|ext|graphs(|[.]decompositions)|interfaces|libs|matrix|misc|numerical(|[.]backends)|rings|sets))[.]all\s+import"
        sage: walkdir_replace_dot_all(dir, fileRegex, regex)
    """
    global numberFiles, numberFilesMatchingRegex, numberFilesChanged, numberStatementsReplaced, log_messages
    pattern = re.compile(fileRegex)
    for root, dirs, files in os.walk(".", topdown=False):
        for name in files:
            numberFiles += 1
            if pattern.search(name):
                numberFilesMatchingRegex += 1
                location = os.path.join(root, name)[1:]
                if location.find('replace_dot_all') == -1:  # to avoid chaning anything in this file itself
                    make_replacements_in_file(dir + location, regex, verbose)
    # sort lines of log_messages by first character of each line(lines are separated by \n)
    sort_log_messages()
    if verbosity:
        print(log_messages)
    report = f'REPORT:\nNumber of files checked: {numberFiles}\nNumber of files matching regex: {numberFilesMatchingRegex}\nNumber of files changed: {numberFilesChanged}\nNumber of import statements replaced: {numberStatementsReplaced}'
    print('*'*100 + '\n' + report + '\n' + '*'*100)


def sort_log_messages():
    r"""
    If the user executes the function walkdir_replace_dot_all with the verbose parameter set to True, then the global variable log_messages will be a string containing all the log messages.
    This function sorts the lines of log_messages by the first character of each line(lines are separated by \n). This function is called at the end of walkdir_replace_dot_all.
    """
    global log_messages
    # split the log messages into a list of strings (each string is a line separated by a newline character)
    log_messages = log_messages.split('\n')
    # sort the list of strings
    log_messages.sort()
    # add index to each line
    for i in range(len(log_messages)):
        log_messages[i] = f'{i}. {log_messages[i]}'
    # join the list of strings into a single string separated by newline characters
    log_messages = '\n'.join(log_messages)[2:]


# ******************************************************** EXECUTES MAIN FUNCTION ****************************************************************************************
# this executes the main function in this file which writes over all import statements matching regex in files in specified location matching fileRegex:
if __name__ == "__main__":
    # Parse the arguments
    args = parseArguments()
    verbosity = args.verbose
    # Print arguments
    print("Executing replace_dot_all.py with arguments:")
    for a in args.__dict__:
        print('    ' + str(a) + ": " + str(args.__dict__[a]))
    # Declare regular expressions
    fileRegex = r'.*[.](py|pyx|pxi)$'
    regex = r"from\s+sage(|[.](arith|categories|combinat|ext|graphs(|[.]decompositions)|interfaces|libs|matrix|misc|numerical(|[.]backends)|rings|sets))[.]all\s+import"
    # Execute the main function based on the specified location and verbosity
    if args.location is None:
        os.chdir(sage.env.SAGE_SRC + '/sage')  # change to sage directory
        dir = os.getcwd()  # Get the current working directory
        walkdir_replace_dot_all(dir, fileRegex, regex, verbose=verbosity)
    elif args.location.find('.py') == -1 and args.location.find('.pxi') == -1:
        os.chdir(sage.env.SAGE_SRC + '/sage/' + args.location)  # change to directory specified by location argument
        dir = os.getcwd()  # Get the current working directory
        walkdir_replace_dot_all(dir, fileRegex, regex, verbose=verbosity)
    elif args.location.find('.py') != -1 or args.location.find('.pxi') != -1:
        # make replacements in file specified by location argument
        make_replacements_in_file(sage.env.SAGE_SRC + '/sage/' + args.location, regex, verbose=verbosity)
# ************************************************************************************************************************************************************************
