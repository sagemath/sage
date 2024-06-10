#!/usr/bin/env sage-python

########################################################################
# Originally based on a script by Andrew Dalke:
#    http://projects.scipy.org/pipermail/numpy-discussion/2008-July/035415.html
#
# 2012: Total rewrite by Volker Braun
########################################################################


########################################################################
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
########################################################################
import sys
import os
import time
import gc
import warnings

# Ignore collections.abc warnings, there are a lot of them but they are
# harmless. These warnings are also disabled in src/sage/all.py.
warnings.filterwarnings('ignore', category=DeprecationWarning,
                        message='.*collections[.]abc.*')

cmdline_args = sys.argv[2:]
have_cmdline_args = bool(cmdline_args)

direct_children_time = 0
import_counter = 0
parent = None
index_to_parent = dict()
all_modules = dict()

DEFAULT_LEVEL = 0


def new_import(name, globals={}, locals={}, fromlist=[], level=DEFAULT_LEVEL):
    """
    The new import function

    Note that ``name`` is not unique, it can be `sage.foo.bar` or `bar`.
    """
    global all_modules, import_counter, parent, direct_children_time
    old_direct_children_time = direct_children_time
    direct_children_time = 0
    old_parent = parent
    parent = this_import_counter = import_counter
    import_counter += 1
    t1 = time.time()
    module = old_import(name, globals, locals, fromlist, level)
    t2 = time.time()
    parent = old_parent
    elapsed_time = t2 - t1
    module_time = elapsed_time - direct_children_time
    direct_children_time = old_direct_children_time + elapsed_time
    index_to_parent[this_import_counter] = module
    data = all_modules.get(module, None)
    if data is not None:
        data['parents'].append(parent)
        data['import_names'].add(name)
        data['cumulative_time'] += elapsed_time
        data['time'] += module_time
        return module
    data = {'cumulative_time': elapsed_time,
            'time': module_time,
            'import_names': set([name]),
            'parents': [parent]}
    all_modules[module] = data
    return module


old_import = __builtins__.__import__
__builtins__.__import__ = new_import
gc.disable()
from sage.all import *
gc.enable()
__builtins__.__import__ = old_import

for data in all_modules.values():
    data['parents'] = set(index_to_parent.get(i, None)
                          for i in data['parents'])


module_by_speed = sorted(((data['time'], module, data)
                          for module, data in all_modules.items()),
                         key=lambda x: x[0])


def print_separator():
    print('=' * 72)


def print_headline(line):
    print('=={0:=<68}=='.format(' ' + line + ' '))


width = 10
fmt_header = '{0:>' + str(width) + '} {1:>' + str(width) + '} {2:>' + str(width) + '}  {3}'
fmt_number = '{0:>' + str(width) + '.3f} {1:>' + str(width) + '.3f} {2:>' + str(width) + '}  {3}'


def print_table(module_list, limit):
    global fmt_header, fmt_number
    print(fmt_header.format('exclude/ms', 'include/ms', '#parents', 'module name'))
    for t, module, data in module_list[-limit:]:
        print(fmt_number.format(1000 * t, 1000 * data['cumulative_time'],
                                len(data['parents']), module.__name__))


def guess_module_name(src):
    module = []
    src, ext = os.path.splitext(src)
    while src and src != '/':
        head, tail = os.path.split(os.path.abspath(src))
        if (tail == 'src' or any(os.path.exists(os.path.join(head, tail, f))
                                 for f in ('setup.py', 'pyproject.toml'))):
            return '.'.join(module)
        module.insert(0, tail)
        src = head
    return None


if not have_cmdline_args:
    print('== Slowest module imports (excluding / including children) ==')
    print_table(module_by_speed, 50)
    print('Total time (sum over exclusive time): {:.3f}ms'.format(1000 * sum(data[0] for data in module_by_speed)))
    print('Use sage -startuptime <module_name|file_name>... to get more details about specific modules.')
else:
    for module_arg in cmdline_args:
        matching_modules = [m for m in all_modules if m.__name__ == module_arg]
        if not matching_modules:
            if '/' in module_arg or any(module_arg.endswith(ext) for ext in ('.py', '.pyx')) or os.path.isdir(module_arg):
                file_name = module_arg
                module_arg = guess_module_name(file_name)
                if not module_arg:
                    print('Warning: "' + file_name + '" does not appear to be a Python module source file or package directory.')
                    continue
                else:
                    matching_modules = [m for m in all_modules if m.__name__.startswith(module_arg)]
            else:
                matching_modules = [m for m in all_modules if m.__name__.endswith(module_arg)]
        if not matching_modules:
            print('Warning: No modules loaded at startup correspond to {}'.format(module_arg))
        for module_name in matching_modules:
            parents = all_modules[module_name]['parents']
            print()
            print_separator()
            print_headline('Slowest modules importing {0}'.format(module_name.__name__))
            print_table([m for m in module_by_speed if m[1] in parents], 10)
            print()
            print_headline('Slowest modules imported by {0}'.format(module_name.__name__))
            print_table([m for m in module_by_speed if module_name in m[2]['parents']], 10)
            print()
            data = all_modules[module_name]
            print_headline('module ' + module_name.__name__)
            print('Time to import:  {0:.3f}ms'.format(1000 * data['time']))
            print('Cumulative time: {0:.3f}ms'.format(1000 * data['cumulative_time']))
            print('Names: {0}'.format(', '.join(data['import_names'])))
            print('File: {0}'.format(module_name.__file__))
