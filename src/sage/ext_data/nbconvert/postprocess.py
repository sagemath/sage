#!/usr/bin/env sage-python

r"""
This postprocess script fixes some issues with .rst files returned by nbconvert:

- it fixes the rendering of the first title (that looks like markdown not rst)
- it removes the header that defines some LaTeX macros
- it removes the :math: role since this is the default one

AUTHORS:

- Thierry Monteil (2018): initial version.
"""

if __name__ == '__main__':

    import sys
    import re

    file_name = sys.argv[1]

    with open(file_name) as f:
            lines = f.readlines()

    # states of the parser
    wrong_math_def_started = False
    wrong_math_fixed = False
    wrong_title_fixed = False

    # processing
    new_file = ''
    for i, line in enumerate(lines):
        if line.startswith(' # ') and not wrong_title_fixed:
            new_file += re.sub('^ # ', '', line)
            new_file += '=' * (len(line) - 4) + '\n'
            wrong_title_fixed = True
        elif line.startswith('.. math::') and not wrong_math_fixed:
            pass
        elif line.startswith('   \\def') and not wrong_math_fixed:
            wrong_math_def_started = True
        elif wrong_math_def_started and not wrong_math_fixed:
            wrong_math_fixed = True
        else:
            new_file += re.sub(':math:', '', line)

    # write new file
    with open(file_name, 'w') as f:
        f.write(new_file)
