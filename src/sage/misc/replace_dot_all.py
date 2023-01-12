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

    ./sage src/sage/misc/replace_dot_all.py

Running replace_dot_all.py will call the function walkdir_replace_dot_all() which walks through all files in
src/sage matching the filePattern above and replaces certain 'from module.all import something' (those matching the pattern above)
with the correct import statement by applying the function import_statements from src/sage/misc/dev_tools.py.

ISSUES: (Note: In order to show examples of bad import statements and have them not be overwritten by this script, just include a * somewhere in the same line. Any lines with * are automatically skipped.)

- PROBLEM (fixed!): for import statements on multiple lines we need to preserve the original indentation after making changes with import_statements

- PROBLEM (fixed!): * from sage.arith.all import (hilbert_conductor_inverse, hilbert_conductor)
We need to deal with these separately i.e. rewrite this as a list of single string imports

-PROBLEM (fixed!): need to handle statements which import a module as a certain name

-PROBLEM (fixed!): combining the two above problems but multiline
e.g. * from sage.arith.all import (hilbert_conductor_inverse,
                                                  hilbert_conductor.
                                                      infinity as oo)

-PROBLEM (fixed!): Statements like
    from sage.rings.rational_field import QQ
can get replaced to
    from sage.rings.rational_field import Q
but then we have files calling QQ when it was imported as Q. This is an issue anytime an object has multiple names. We should have a
way of referencing the original way the object was called and replace it in the statement given from the function import_statements by hand.

-PROBLEM (fixed!): need to handle exceptions better in find_replacements

-PROBLEM: We currently only fix these issues in locations where 'from' are the first four characters of the line after any whitespace.
this excludes lines in documentation where an example is given of the form "sage: from module.submod.all import something". Current examples:

sage.env.SAGE_SRC/sage/categories/category.py line number 156.
*       sage: from sage.categories.all import Category
sage.env.SAGE_SRC/sage/categories/category.py line number 203.
*       sage: from sage.categories.all import Category
sage.env.SAGE_SRC/sage/misc/cachefunc.pyx line number 72.
*           ....: 'from sage.all import cached_method',
sage.env.SAGE_SRC/sage/misc/cachefunc.pyx line number 109.
*   sage: cython_code = ["from sage.all import cached_method, cached_in_parent_method, Category, Objects",
sage.env.SAGE_SRC/sage/misc/lazy_import.pyx line number 1043.
*       Importing my_Qp from here is deprecated; please use "from sage.all import Qp as my_Qp" instead.
sage.env.SAGE_SRC/sage/numerical/linear_functions.pyx line number 1067.
*           sage: from sage.rings.all import AA
sage.env.SAGE_SRC/sage/structure/element.pyx line number 2230.
*       sage: cython_code = ["from sage.all import cached_method, cached_in_parent_method, Category, Objects",
sage.env.SAGE_SRC/sage/structure/element.pyx line number 2328.
*            ....: from sage.all import cached_method, lazy_attribute, Category, Objects
sage.env.SAGE_SRC/sage/tests/books/computational-mathematics-with-sagemath/premierspas_doctest.py line number 120.
*  sage: from sage.all import pi
sage.env.SAGE_SRC/sage/env.py line number 14.
*   sage: cmd = "from sage.all import SAGE_ROOT, SAGE_LOCAL; print((SAGE_ROOT, SAGE_LOCAL))"


PROBLEM: *  from sage.arith.all import LCM       is valid but this gets replaced to   from sage.arith.functions import lcm      which span decided to replace by    from sage.arith.functions import LCM
* since in this file, they will be calling the function by the name LCM, the way they imported it. However, the statement     from sage.arith.functions import LCM      is not valid.
We get the error:

Traceback (most recent call last):
  File "sage.env.SAGE_SRC/sage/misc/replace_dot_all.py", line 89, in <module>
    from sage.all import * # which one is better to use?
  File "sage.env.SAGE_SRC/sage/all.py", line 135, in <module>
    from sage.rings.all      import *
  File "sage.env.SAGE_SRC/sage/rings/all.py", line 79, in <module>
    from .number_field.all import *
  File "sage.env.SAGE_SRC/sage/rings/number_field/all.py", line 2, in <module>
    from .number_field import (NumberField, NumberFieldTower, CyclotomicField, QuadraticField,
  File "sage.env.SAGE_SRC/sage/rings/number_field/number_field.py", line 137, in <module>
    from .unit_group import UnitGroup
  File "sage.env.SAGE_SRC/sage/rings/number_field/unit_group.py", line 163, in <module>
    from sage.groups.abelian_gps.values import AbelianGroupWithValues_class
  File "sage.env.SAGE_SRC/sage/groups/abelian_gps/values.py", line 76, in <module>
    from sage.groups.abelian_gps.abelian_group import AbelianGroup_class, _normalize
  File "sage.env.SAGE_SRC/sage/groups/abelian_gps/abelian_group.py", line 213, in <module>
    from sage.groups.abelian_gps.abelian_group_element import AbelianGroupElement
  File "sage.env.SAGE_SRC/sage/groups/abelian_gps/abelian_group_element.py", line 47, in <module>
    from sage.groups.abelian_gps.element_base import AbelianGroupElementBase
  File "sage.env.SAGE_SRC/sage/groups/abelian_gps/element_base.py", line 24, in <module>
    from sage.arith.functions import LCM
ImportError: cannot import name 'LCM' from 'sage.arith.functions' (sage.env.SAGE_SRC/sage/arith/functions.cpython-310-x86_64-linux-gnu.so)

* The issue is that the name was originally LCM and we want to preserve that so LCM can still be called in the given file. But somehow you can import LCM from
* sage.arith.all but not from sage.arith.functions
* we need to call from sage.arith.functions import lcm as LCM so instead of replacing
change_to_temp = change_to_temp.replace(new_mod_as_string,original_mod_string)
we should do
change_to_temp = change_to_temp.replace(new_mod_as_string,new_mod_as_string + ' as ' + original_mod_string) # originally thought it would be best to make these replacements but it causes strange issues


PROBLEM: for example in ell_rational_field.py we are now getting

from sage.rings.infinity import Infinity as infinity as oo

because we have two different cases where we use the as statement and in this case we are doing both.
"""

# Importing packages

from sage.misc.dev_tools import import_statements
# import sage.all
from sage.all import *  # which one is better to use?
import os
import re

# Global variables

examples = list('ABCDEFGHIJ')  # controls how we print out interesting examples to the console
interesting_examples = dict(zip(examples, [0]*len(examples)))
number_examples_to_print = 3
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

        python3: os.chdir('sage.env.SAGE_SRC/sage') # change to sage directory
        python3: cwd = os.getcwd() # Get the current working directory
        python3: file_to_change = 'structure/element.pyx'
        python3: location = cwd + file_to_change
        python3: replacements = find_replacements(location)
    """

    pattern = re.compile(regex)
    replacements = []
    with open(location, "r+") as fp:
        lines = fp.readlines()  # read all lines using readline()
        row_index = 0
        for row in lines:  # iterate over all lines of python file
            if pattern.search(row):  # (match the regex also do not want to mess with documentation)
                prefix = ''
                if '*' in row or 'SAGE_ROOT' in row:
                    if verbose:
                        print(
                            f'J. Match but no changes made (import statement uses *) at {location} line number {row_index + 1}. Not applying any changes here.')
                    continue
                elif not (row.lstrip()[0:4] == 'from'):
                    if '"' not in row and "'" not in row:
                        print(f'H. Interesting example (line with import statement does not start with "from") at {location} line number {row_index + 1}.')
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
                            f'I. Interesting example (import statement does not start with "from") at {location} line number {row_index + 1}. Not yet implemented...')
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
                            print(f'D. Interesting example (spans multiple lines and has " as ") at {location} line number {row_index + 1}')
                            interesting_examples['D'] += 1
                        elif interesting_examples['B'] < number_examples_to_print:
                            print(f'B. Interesting example (spans multiple lines) at {location} line number {row_index + 1}')
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
                        print(f'printing to_exec: {to_exec}')
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
                            print(f'C. Interesting example (" as " in tuple import) at {location} at line number {row_index + 1}')
                            interesting_examples['C'] += 1
                        as_index = to_eval_list_raw[to_eval_list_index].index(" as ")
                        postfix = to_eval_list_raw[to_eval_list_index][as_index:]
                        print(f'postfix at {location} at line number {row_index + 1}: {postfix}')
                    change_to_temp = import_statements(mod, answer_as_str=True, verbose=False)  # import statement for the current mod in the list module
                    import_index = change_to_temp.find('import')
                    new_mod_as_string = change_to_temp[import_index + 7:].strip()  # the name for the module given by the function import_statements
                    if as_index >= 0:
                        # the name for the module as originally called in the document (when there is an " as " statement)
                        original_mod_string = to_eval_list_raw[to_eval_list_index].strip()[:as_index]
                    else:
                        original_mod_string = to_eval_list[to_eval_list_index].strip()  # the name for the module as originally called in the document
                    if original_mod_string != new_mod_as_string:  # if the names differ, we use the original name as it was called in the document
                        if verbose and interesting_examples['A'] < number_examples_to_print:
                            print(
                                f'A. Interesting example (module has multiple names) at {location} line number {row_index + 1}. Names: {original_mod_string}, {new_mod_as_string}. Replacing new {new_mod_as_string} by original {original_mod_string}.')
                            interesting_examples['A'] += 1
                        rep = new_mod_as_string + ' as ' + original_mod_string
                        print(f'changing {new_mod_as_string} to {rep} in {change_to_temp} at {location} at line number {row_index + 1}')
                        change_to_temp = change_to_temp.replace(' '+new_mod_as_string, ' '+new_mod_as_string + ' as ' + original_mod_string)
                        print(f'the result is {change_to_temp}')
                        # change_to_temp = change_to_temp.replace(new_mod_as_string,original_mod_string) # originally this was the replacement but changed to the above to fix issues
                        if " as " in postfix and interesting_examples['G'] < number_examples_to_print:
                            print(
                                f'G. Interesting example (module has multiple names) at {location} line number {row_index + 1}. Names: {original_mod_string}, {new_mod_as_string}. Replacing new {new_mod_as_string} by original {original_mod_string}.')
                            interesting_examples['G'] += 1
                    if len(postfix.strip()) > 0:  # if module was called with " as " statement, we put that back in by adding the string "postfix"
                        # if " as " in change_to_temp locate the index of " as ", remove the end after this, and add the postfix there
                        if " as " in change_to_temp:
                            if verbose:
                                print(f'adding postfix {postfix} to {change_to_temp} after stripping existing "as" statement')
                            change_to_temp = change_to_temp[:change_to_temp.index(" as ")] + ' ' + postfix.strip()
                        else:
                            if verbose:
                                print(f'adding postfix {postfix} to {change_to_temp}')
                            change_to_temp += (' ' + postfix.strip())
                    change_to += (prefix + change_to_temp + '\n')
                    to_eval_list_index += 1
                # [:-1] on change_to gets rid of the last '\n' we added which adds an unnecessary new line
                replacement = [row_index, import_index, change_to[:-1]].copy()
                if span > 0:  # if original statement spanned multiple lines, we store that information to signal that we need to skip lines as we read the document in the function make_replacements_in_file
                    replacement.append(span)
                replacements.append(replacement)
            row_index += 1
    # to print a statement referencing each file which a change was made to
    # line_numbers = [replacement[0] for replacement in replacements]
    # if len(line_numbers)>0:
        # print(f'file {location} contains .all statements at lines {line_numbers}')
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

        python3: os.chdir('sage.env.SAGE_SRC/sage') # change to sage directory
        python3: cwd = os.getcwd() # Get the current working directory
        python3: file_to_change = 'structure/element.pyx'
        python3: location = cwd + file_to_change
        python3: replacements = find_replacements(location)
        python3: file = open(location, "r")
        python3: lines = file.readlines()
        python3: row_index,import_index = replacements[0],replacements[1]
        python3: line = lines[row_index]
        python3: print(f'old line, old reps: {line,replacements}')
        python3: new_line,replacements = process_line(location,line,replacements,row_index)
        python3: print(f'new line, new reps: {new_line,replacements}')
    """

    line = line.rstrip()  # stripping line break
    new_line = ''
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
                print(f'F. Interesting example (has leading space and multiline) at {location} row number {replacement[0] + 1}')
                interesting_examples['F'] += 1
            elif interesting_examples['E'] < number_examples_to_print:
                print(f'E. Interesting example (has leading space) at {location} row number {replacement[0] + 1}')
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
        python3: os.chdir('sage.env.SAGE_SRC/sage') # change to sage directory
        python3: cwd = os.getcwd() # Get the current working directory
        python3: file_to_change = 'structure/element.pyx'
        python3: location = cwd + file_to_change
        python3: make_replacements_in_file(location)
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
        python3: os.chdir('sage.env.SAGE_SRC/sage') # change to sage directory
        python3: cwd = os.getcwd() # Get the current working directory
        python3: walkdir_replace_dot_all()
    """
    global numberFiles, numberFilesMatchingRegex, numberFilesChanged, numberStatementsReplaced
    pattern = re.compile(fileRegex)
    for root, dirs, files in os.walk(".", topdown=False):
        for name in files:
            numberFiles += 1
            if pattern.search(name):
                numberFilesMatchingRegex += 1
                location = os.path.join(root, name)[1:]
                make_replacements_in_file(dir + location, regex, verbose)
    report = f'REPORT:\nNumber of files checked: {numberFiles}\nNumber of files matching regex: {numberFilesMatchingRegex}\nNumber of files changed: {numberFilesChanged}\nNumber of import statements replaced: {numberStatementsReplaced}'
    print('*'*100 + '\n' + report + '\n' + '*'*100)


# ******************************************************** EXECUTES MAIN FUNCTION ****************************************************************************************
# this executes the main function in this file which writes over all import statements matching regex in files in src/sage matching fileRegex:
if __name__ == "__main__":
    fileRegex = r'.*[.](py|pyx|pxi)$'
    regex = r"from\s+sage(|[.](arith|categories|combinat|ext|graphs(|[.]decompositions)|interfaces|libs|matrix|misc|numerical(|[.]backends)|rings|sets))[.]all\s+import"
    os.chdir(sage.env.SAGE_SRC + '/sage')  # change to sage directory
    os.chdir(sage.env.SAGE_SRC + '/sage/coding')  # change to a more specific sage directory if desired
    dir = os.getcwd()  # Get the current working directory
    walkdir_replace_dot_all(dir, fileRegex, regex, verbose=True)
# ************************************************************************************************************************************************************************
