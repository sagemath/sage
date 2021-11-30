r"""
Keep track of attached files

TESTS::

    sage: attach('http://wstein.org/loadtest.py')
    Traceback (most recent call last):
    ...
    NotImplementedError: you cannot attach a URL

Check that no file clutter is produced::

    sage: dir = tmp_dir()
    sage: src = os.path.join(dir, 'foobar.sage')
    sage: with open(src, 'w') as f:
    ....:     _ = f.write('print("<output from attached file>")\n')
    sage: attach(src)
    <output from attached file>
    sage: os.listdir(dir)
    ['foobar.sage']
    sage: detach(src)

In debug mode backtraces contain code snippets. We need to manually
print the traceback because the python doctest module has special
support for exceptions and does not match them
character-by-character::

    sage: import traceback
    sage: with open(src, 'w') as f:
    ....:     _ = f.write('# first line\n')
    ....:     _ = f.write('# second line\n')
    ....:     _ = f.write('raise ValueError("third")   # this should appear in the source snippet\n')
    ....:     _ = f.write('# fourth line\n')

    sage: load_attach_mode(attach_debug=False)
    sage: try:
    ....:     attach(src)
    ....: except Exception:
    ....:     traceback.print_exc(file=sys.stdout)
    Traceback (most recent call last):
    ...
        exec(preparse_file(f.read()) + "\n", globals)
      File "<string>", line 3, in <module>
    ValueError: third
    sage: detach(src)

    sage: load_attach_mode(attach_debug=True)
    sage: try:
    ....:     attach(src)
    ....: except Exception:
    ....:     traceback.print_exc(file=sys.stdout)
    Traceback (most recent call last):
    ...
        exec(code, globals)
      File ".../foobar.sage....py", line ..., in <module>
        raise ValueError("third")   # this should appear in the source snippet
    ValueError: third
    sage: detach(src)
"""

# ****************************************************************************
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import os
import time
from IPython import get_ipython

from sage.repl.load import load, load_wrap
import sage.repl.inputhook
import sage.env

# The attached files as a dict of {filename:mtime}
attached = {}


load_debug_mode = False
attach_debug_mode = True


def load_attach_mode(load_debug=None, attach_debug=None):
    """
    Get or modify the current debug mode for the behavior of
    :func:`load` and :func:`attach` on ``.sage`` files.

    In debug mode, loaded or attached ``.sage`` files are preparsed
    through a file to make their tracebacks more informative. If not
    in debug mode, then ``.sage`` files are preparsed in memory only
    for performance.

    At startup, debug mode is ``True`` for attaching and ``False``
    for loading.

    .. NOTE::

        This function should really be deprecated and code executed
        from memory should raise proper tracebacks.

    INPUT:

    - ``load_debug`` -- boolean or ``None`` (default); if not
      ``None``, then set a new value for the debug mode for loading
      files.

    - ``attach_debug`` -- boolean or ``None`` (default); same as
      ``load_debug``, but for attaching files.

    OUTPUT:

    If all input values are ``None``, returns a tuple giving the
    current modes for loading and attaching.

    EXAMPLES::

        sage: load_attach_mode()
        (False, True)
        sage: load_attach_mode(attach_debug=False)
        sage: load_attach_mode()
        (False, False)
        sage: load_attach_mode(load_debug=True)
        sage: load_attach_mode()
        (True, False)
        sage: load_attach_mode(load_debug=False, attach_debug=True)
    """
    global load_debug_mode, attach_debug_mode
    if load_debug is None and attach_debug is None:
        return (load_debug_mode, attach_debug_mode)
    if load_debug is not None:
        load_debug_mode = load_debug
    if attach_debug is not None:
        attach_debug_mode = attach_debug


search_paths = []


def load_attach_path(path=None, replace=False):
    """
    Get or modify the current search path for :func:`load` and
    :func:`attach`.

    INPUT:

    - ``path`` -- string or list of strings (default: ``None``);
      path(s) to append to or replace the current path.

    - ``replace`` -- boolean (default: ``False``); if ``path`` is not
      ``None``, whether to *replace* the search path instead of
      *appending* to it.

    OUTPUT:

    ``None`` or a *reference* to the current search paths.

    EXAMPLES:

    First, we extend the example given in :func:`load`'s docstring::

        sage: sage.repl.attach.reset(); reset_load_attach_path()
        sage: load_attach_path()
        ['.']
        sage: t_dir = tmp_dir()
        sage: fullpath = os.path.join(t_dir, 'test.py')
        sage: with open(fullpath, 'w') as f:
        ....:     _ = f.write("print(37 * 3)")

    We put ``SAGE_TMP`` on the attach path for testing (otherwise this will
    load ``test.py`` from the current working directory if that happens
    to exist)::

        sage: load_attach_path(SAGE_TMP, replace=True)
        sage: attach('test.py')
        Traceback (most recent call last):
        ...
        OSError: did not find file 'test.py' to load or attach
        sage: load_attach_path(t_dir)
        sage: attach('test.py')
        111
        sage: attached_files() == [fullpath]
        True
        sage: sage.repl.attach.reset(); reset_load_attach_path()
        sage: load_attach_path() == ['.']
        True
        sage: load_attach_path(SAGE_TMP, replace=True)
        sage: load('test.py')
        Traceback (most recent call last):
        ...
        OSError: did not find file 'test.py' to load or attach

    The function returns a reference to the path list::

        sage: reset_load_attach_path(); load_attach_path()
        ['.']
        sage: load_attach_path('/path/to/my/sage/scripts'); load_attach_path()
        ['.', '/path/to/my/sage/scripts']
        sage: load_attach_path(['good', 'bad', 'ugly'], replace=True)
        sage: load_attach_path()
        ['good', 'bad', 'ugly']
        sage: p = load_attach_path(); p.pop()
        'ugly'
        sage: p[0] = 'weird'; load_attach_path()
        ['weird', 'bad']
        sage: reset_load_attach_path(); load_attach_path()
        ['.']
    """
    global search_paths
    if path is None:
        return search_paths
    else:
        if not isinstance(path, list):
            path = [path]
        if replace:
            search_paths = path
        else:
            for p in path:
                if not p:
                    continue
                if p not in search_paths:
                    search_paths.append(p)


def reset_load_attach_path():
    """
    Reset the current search path for :func:`load` and :func:`attach`.

    The default path is ``'.'`` plus any paths specified in the
    environment variable ``SAGE_LOAD_ATTACH_PATH``.

    EXAMPLES::

        sage: load_attach_path()
        ['.']
        sage: t_dir = tmp_dir()
        sage: load_attach_path(t_dir)
        sage: t_dir in load_attach_path()
        True
        sage: reset_load_attach_path(); load_attach_path()
        ['.']

    At startup, Sage adds colon-separated paths in the environment
    variable ``SAGE_LOAD_ATTACH_PATH``::

        sage: reset_load_attach_path(); load_attach_path()
        ['.']
        sage: os.environ['SAGE_LOAD_ATTACH_PATH'] = '/veni/vidi:vici:'
        sage: from importlib import reload
        sage: reload(sage.repl.attach)    # Simulate startup
        <module 'sage.repl.attach' from '...'>
        sage: load_attach_path()
        ['.', '/veni/vidi', 'vici']
        sage: del os.environ['SAGE_LOAD_ATTACH_PATH']
        sage: reload(sage.repl.preparse)    # Simulate startup
        <module 'sage.repl.preparse' from '...'>
        sage: reset_load_attach_path(); load_attach_path()
        ['.']
    """
    global search_paths
    search_paths = ['.']
    for path in os.environ.get('SAGE_LOAD_ATTACH_PATH', '').split(':'):
        load_attach_path(path=path)


# Set up the initial search path for loading and attaching files.  A
# user can modify the path with the function load_attach_path.
reset_load_attach_path()


def attach(*files):
    """
    Attach a file or files to a running instance of Sage and also load
    that file.

    .. NOTE::

        Attaching files uses the Python inputhook, which will conflict
        with other inputhook users. This generally includes GUI main loop
        integrations, for example tkinter. So you can only use tkinter or
        attach, but not both at the same time.

    INPUT:

    - ``files`` -- a list of filenames (strings) to attach.

    OUTPUT:

    Each file is read in and added to an internal list of watched files.
    The meaning of reading in a file depends on the file type:

    -  ``.py`` files are read in with no preparsing (so, e.g., ``2^3`` is 2
       bit-xor 3);

    -  ``.sage`` files are preparsed, then the result is read in;

    - ``.pyx`` files are *not* preparsed, but rather are compiled to a
       module ``m`` and then ``from m import *`` is executed.

    The contents of the file are then loaded, which means they are read
    into the running Sage session. For example, if ``foo.sage`` contains
    ``x=5``, after attaching ``foo.sage`` the variable ``x`` will be set
    to 5. Moreover, any time you change ``foo.sage``, before you execute
    a command, the attached file will be re-read automatically (with no
    intervention on your part).

    .. SEEALSO::

        :meth:`~sage.repl.load.load` is the same as :func:`attach`, but
        does not automatically reload a file when it changes.

    EXAMPLES:

    You attach a file, e.g., ``foo.sage`` or ``foo.py`` or
    ``foo.pyx``, to a running Sage session by typing::

        sage: attach('foo.sage')  # not tested

    Here we test attaching multiple files at once::

        sage: sage.repl.attach.reset()
        sage: t1 = tmp_filename(ext='.py')
        sage: with open(t1,'w') as f: _ = f.write("print('hello world')")
        sage: t2 = tmp_filename(ext='.py')
        sage: with open(t2,'w') as f: _ = f.write("print('hi there xxx')")
        sage: attach(t1, t2)
        hello world
        hi there xxx
        sage: set(attached_files()) == set([t1,t2])
        True

    .. SEEALSO::

        - :meth:`attached_files` returns a list of
          all currently attached files.

        - :meth:`detach` instructs Sage to remove a
          file from the internal list of watched files.

        - :meth:`load_attach_path` allows you to
          get or modify the current search path for loading and attaching
          files.
    """
    try:
        ipy = get_ipython()
    except NameError:
        ipy = None
    global attached
    for filename in files:
        if ipy:
            code = load_wrap(filename, attach=True)
            ipy.run_cell(code)
        else:
            load(filename, globals(), attach=True)


def add_attached_file(filename):
    """
    Add to the list of attached files

    This is a callback to be used from
    :func:`~sage.repl.load.load` after evaluating the attached
    file the first time.

    INPUT:

    - ``filename`` -- string, the fully qualified file name.

    EXAMPLES::

        sage: import sage.repl.attach as af
        sage: af.reset()
        sage: t = tmp_filename(ext='.py')
        sage: af.add_attached_file(t)
        sage: af.attached_files()
        ['/.../tmp_....py']
        sage: af.detach(t)
        sage: af.attached_files()
        []
    """
    sage.repl.inputhook.install()
    fpath = os.path.abspath(filename)
    attached[fpath] = os.path.getmtime(fpath)


def attached_files():
    """
    Returns a list of all files attached to the current session with
    :meth:`attach`.

    OUTPUT:

    The filenames in a sorted list of strings.

    EXAMPLES::

        sage: sage.repl.attach.reset()
        sage: t = tmp_filename(ext='.py')
        sage: with open(t,'w') as f: _ = f.write("print('hello world')")
        sage: attach(t)
        hello world
        sage: attached_files()
        ['/....py']
        sage: attached_files() == [t]
        True
    """
    global attached
    return sorted(attached)


def detach(filename):
    """
    Detach a file.

    This is the counterpart to :meth:`attach`.

    INPUT:

    - ``filename`` -- a string, or a list of strings, or a tuple of strings.

    EXAMPLES::

        sage: sage.repl.attach.reset()
        sage: t = tmp_filename(ext='.py')
        sage: with open(t,'w') as f: _ = f.write("print('hello world')")
        sage: attach(t)
        hello world
        sage: attached_files() == [t]
        True
        sage: detach(t)
        sage: attached_files()
        []

        sage: sage.repl.attach.reset(); reset_load_attach_path()
        sage: load_attach_path()
        ['.']
        sage: t_dir = tmp_dir()
        sage: fullpath = os.path.join(t_dir, 'test.py')
        sage: with open(fullpath, 'w') as f: _ = f.write("print(37 * 3)")
        sage: load_attach_path(t_dir, replace=True)
        sage: attach('test.py')
        111
        sage: attached_files() == [os.path.normpath(fullpath)]
        True
        sage: detach('test.py')
        sage: attached_files()
        []
        sage: attach('test.py')
        111
        sage: fullpath = os.path.join(t_dir, 'test2.py')
        sage: with open(fullpath, 'w') as f: _ = f.write("print(3)")
        sage: attach('test2.py')
        3
        sage: detach(attached_files())
        sage: attached_files()
        []

    TESTS::

        sage: detach('/dev/null/foobar.sage')
        Traceback (most recent call last):
        ...
        ValueError: file '/dev/null/foobar.sage' is not attached, see attached_files()
    """
    if isinstance(filename, str):
        filelist = [filename]
    else:
        filelist = [str(x) for x in filename]

    global attached
    for filename in filelist:
        fpath = os.path.expanduser(filename)
        if not os.path.isabs(fpath):
            for path in load_attach_path():
                epath = os.path.expanduser(path)
                fpath = os.path.join(epath, filename)
                fpath = os.path.abspath(fpath)
                if fpath in attached:
                    break
        if fpath in attached:
            attached.pop(fpath)
        else:
            raise ValueError("file '{0}' is not attached, see attached_files()".format(filename))
    if not attached:
        sage.repl.inputhook.uninstall()


def reset():
    """
    Remove all the attached files from the list of attached files.

    EXAMPLES::

        sage: sage.repl.attach.reset()
        sage: t = tmp_filename(ext='.py')
        sage: with open(t,'w') as f: _ = f.write("print('hello world')")
        sage: attach(t)
        hello world
        sage: attached_files() == [t]
        True
        sage: sage.repl.attach.reset()
        sage: attached_files()
        []
    """
    global attached
    attached = {}


def modified_file_iterator():
    """
    Iterate over the changed files

    As a side effect the stored time stamps are updated with the
    actual time stamps. So if you iterate over the attached files in
    order to reload them and you hit an error then the subsequent
    files are not marked as read.

    Files that are in the process of being saved are excluded.

    EXAMPLES::

        sage: sage.repl.attach.reset()
        sage: t = tmp_filename(ext='.py')
        sage: attach(t)
        sage: from sage.repl.attach import modified_file_iterator
        sage: list(modified_file_iterator())
        []
        sage: sleep(1)   # filesystem mtime granularity
        sage: with open(t, 'w') as f: _ = f.write('1')
        sage: list(modified_file_iterator())
        [('/.../tmp_....py', time.struct_time(...))]
    """
    global attached
    modified = {}
    for filename in list(attached):
        old_tm = attached[filename]
        if not os.path.exists(filename):
            print('### detaching file {0} because it does not exist (deleted?) ###'.format(filename))
            detach(filename)
            continue
        new_tm = os.path.getmtime(filename)
        if new_tm > old_tm:
            modified[filename] = new_tm

    if not modified:
        return
    time.sleep(0.1)  # sleep 100ms to give the editor time to finish saving

    for filename in list(modified):
        old_tm = modified[filename]
        new_tm = os.path.getmtime(filename)
        if new_tm == old_tm:
            # file was modified but did not change in the last 100ms
            attached[filename] = new_tm
            yield filename, time.gmtime(new_tm)


def reload_attached_files_if_modified():
    r"""
    Reload attached files that have been modified

    This is the internal implementation of the attach mechanism.

    EXAMPLES::

        sage: sage.repl.attach.reset()
        sage: from sage.repl.interpreter import get_test_shell
        sage: shell = get_test_shell()
        sage: tmp = tmp_filename(ext='.py')
        sage: with open(tmp, 'w') as f: _ = f.write('a = 2\n')
        sage: shell.run_cell('attach({0})'.format(repr(tmp)))
        sage: shell.run_cell('a')
        2
        sage: sleep(1)   # filesystem mtime granularity
        sage: with open(tmp, 'w') as f: _ = f.write('a = 3\n')

    Note that the doctests are never really at the command prompt
    where the automatic reload is triggered. So we have to do it
    manually::

        sage: shell.run_cell('from sage.repl.attach import reload_attached_files_if_modified')
        sage: shell.run_cell('reload_attached_files_if_modified()')
        ### reloading attached file tmp_....py modified at ... ###

        sage: shell.run_cell('a')
        3
        sage: shell.run_cell('detach({0})'.format(repr(tmp)))
        sage: shell.run_cell('attached_files()')
        []
        sage: shell.quit()
    """
    ip = get_ipython()
    for filename, mtime in modified_file_iterator():
        basename = os.path.basename(filename)
        timestr = time.strftime('%T', mtime)
        notice = '### reloading attached file {0} modified at {1} ###'.format(basename, timestr)
        if ip:
            print(notice)
            code = load_wrap(filename, attach=True)
            ip.run_cell(code)
        else:
            print(notice)
            load(filename, globals(), attach=True)
