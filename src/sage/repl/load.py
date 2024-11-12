# sage_setup: distribution = sagemath-repl
"""
Load Python, Sage, Cython, Fortran and Magma files in Sage
"""
# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import base64
from pathlib import Path

from sage.cpython.string import str_to_bytes, bytes_to_str, FS_ENCODING


def is_loadable_filename(filename):
    """
    Return whether a file can be loaded into Sage.

    This checks only whether its name ends in one of the supported
    extensions ``.py``, ``.pyx``, ``.sage``, ``.spyx``, ``.f``,
    ``.f90`` and ``.m``.

    .. NOTE:: :func:`load` assumes that `.m` signifies a Magma file.

    INPUT:

    - ``filename`` -- string or :class:`Path` object

    OUTPUT: boolean

    EXAMPLES::

        sage: sage.repl.load.is_loadable_filename('foo.bar')
        False
        sage: sage.repl.load.is_loadable_filename('foo.c')
        False
        sage: sage.repl.load.is_loadable_filename('foo.sage')
        True
        sage: sage.repl.load.is_loadable_filename('FOO.F90')
        True
        sage: sage.repl.load.is_loadable_filename('foo.m')
        True

        sage: from pathlib import Path
        sage: sage.repl.load.is_loadable_filename(Path('foo.py'))
        True
    """
    ext = Path(filename).suffix.lower()
    return ext in ('.py', '.pyx', '.sage', '.spyx', '.f', '.f90', '.m')


def load_cython(name):
    """
    Helper function to load a Cython file.

    INPUT:

    - ``name`` -- filename of the Cython file

    OUTPUT:

    - A string with Python code to import the names from the compiled
      module.
    """
    from sage.misc.cython import cython
    mod, dir = cython(str(name), compile_message=True, use_cache=True)
    import sys
    sys.path.append(dir)
    return f'from {mod} import *'


def load(filename, globals, attach=False):
    r"""
    Execute a file in the scope given by ``globals``. If the name starts with
    ``http://`` or ``https://``, it is treated as a URL and downloaded.

    .. NOTE::

        For Cython files, the situation is more complicated --
        the module is first compiled to a temporary module ``t`` and
        executed via::

            from t import *

    .. NOTE::

        The global ``load`` function is :func:`sage.misc.persist.load`,
        which delegates to this function for code file formats.

        ``%runfile`` magic can also be used, see
        :meth:`~sage.repl.ipython_extension.SageMagics.runfile`.

    INPUT:

    - ``filename`` -- string (denoting a filename or URL) or a :class:`Path` object

    - ``globals`` -- string:object dictionary; the context in which
      to execute the file contents

    - ``attach`` -- boolean (default: ``False``); whether to add the
      file to the list of attached files

    Loading an executable Sage script from the :ref:`command line <section-command-line>`
    will run whatever code is inside an

    ::

        if __name__ == "__main__":

    section, as the condition on ``__name__`` will hold true (code run from the
    command line is considered to be running in the ``__main__`` module.)

    EXAMPLES:

    Note that ``.py`` files are *not* preparsed::

        sage: t = tmp_filename(ext='.py')
        sage: with open(t, 'w') as f:
        ....:     _ = f.write("print(('hi', 2^3)); z = -2^7")
        sage: z = 1
        sage: sage.repl.load.load(t, globals())
        ('hi', 1)
        sage: z
        -7

    A ``.sage`` file *is* preparsed::

        sage: t = tmp_filename(ext='.sage')
        sage: with open(t, 'w') as f:
        ....:     _ = f.write("print(('hi', 2^3)); z = -2^7")
        sage: z = 1
        sage: sage.repl.load.load(t, globals())
        ('hi', 8)
        sage: z
        -128

    Cython files are *not* preparsed::

        sage: t = tmp_filename(ext='.pyx')
        sage: with open(t, 'w') as f:
        ....:     _ = f.write("print(('hi', 2^3)); z = -2^7")
        sage: z = 1
        sage: sage.repl.load.load(t, globals())                                         # needs sage.misc.cython
        Compiling ...
        ('hi', 1)
        sage: z
        -7

    If the file is not a Cython, Python, or Sage file, a :exc:`ValueError`
    is raised::

        sage: sage.repl.load.load(tmp_filename(ext='.foo'), globals())
        Traceback (most recent call last):
        ...
        ValueError: unknown file extension '.foo' for load or attach (supported extensions: .py, .pyx, .sage, .spyx, .f, .f90, .m)

    We load a file given at a remote URL (not tested for security reasons)::

        sage: sage.repl.load.load('https://www.sagemath.org/files/loadtest.py', globals())  # not tested
        hi from the net
        5

    We can load files using secure http (https)::

        sage: sage.repl.load.load('https://raw.githubusercontent.com/sagemath/sage-patchbot/3.0.0/sage_patchbot/util.py', globals())  # optional - internet

    We attach a file (note that :func:`~sage.repl.attach.attach`
    is equivalent, but available at the global scope by default)::

        sage: t = tmp_filename(ext='.py')
        sage: with open(t, 'w') as f:
        ....:     _ = f.write("print('hello world')")
        sage: sage.repl.load.load(t, globals(), attach=True)
        hello world
        sage: t in attached_files()
        True

    You cannot attach remote URLs (yet)::

        sage: sage.repl.load.load('https://www.sagemath.org/files/loadtest.py', globals(), attach=True)  # optional - internet
        Traceback (most recent call last):
        ...
        NotImplementedError: you cannot attach a URL

    The default search path for loading and attaching files is the
    current working directory, i.e., ``'.'``.  But you can modify the
    path with :func:`load_attach_path`::

        sage: import tempfile
        sage: sage.repl.attach.reset(); reset_load_attach_path()
        sage: load_attach_path()
        [PosixPath('.')]
        sage: with tempfile.TemporaryDirectory() as t_dir:
        ....:     fname = 'test.py'
        ....:     fullpath = os.path.join(t_dir, fname)
        ....:     with open(fullpath, 'w') as f:
        ....:         _ = f.write("print(37 * 3)")
        ....:     load_attach_path(t_dir, replace=True)
        ....:     attach(fname)
        111
        sage: sage.repl.attach.reset(); reset_load_attach_path() # clean up

    or by setting the environment variable ``SAGE_LOAD_ATTACH_PATH``
    to a colon-separated list before starting Sage::

        $ export SAGE_LOAD_ATTACH_PATH="/path/to/my/library:/path/to/utils"
        $ sage
        sage: load_attach_path()          # not tested
        ['.', '/path/to/my/library', '/path/to/utils']

    TESTS:

    Make sure that load handles filenames with spaces in the name or path::

        sage: t = tmp_filename(ext=' b.sage')
        sage: with open(t, 'w') as f:
        ....:     _ = f.write("print(2)")
        sage: sage.repl.load.load(t, globals())
        2

    Non-existing files with spaces give correct messages::

        sage: sage.repl.load.load("this file should not exist", globals())
        Traceback (most recent call last):
        ...
        OSError: did not find file 'this file should not exist' to load or attach
    """
    if attach:
        from sage.repl.attach import add_attached_file

    if isinstance(filename, bytes):
        # For Python 3 in particular, convert bytes filenames to str since the
        # rest of this functions operate on filename as a str
        filename = bytes_to_str(filename, FS_ENCODING, 'surrogateescape')

    if isinstance(filename, str) and filename.lower().startswith(('http://', 'https://')):
        if attach:
            # But see https://en.wikipedia.org/wiki/HTTP_ETag for how
            # we will do this.
            # https://diveintopython3.net/http-web-services.html#etags
            raise NotImplementedError("you cannot attach a URL")
        from sage.misc.remote_file import get_remote_file
        filename = get_remote_file(filename, verbose=False)

    filename = Path(filename).expanduser()

    from sage.repl.attach import load_attach_path
    for path in load_attach_path():
        fpath = (path / filename).expanduser()
        if fpath.is_file():
            break
    else:
        raise OSError('did not find file %r to load or attach' % str(filename))

    ext = fpath.suffix.lower()
    if ext == '.py':
        if attach:
            add_attached_file(fpath)
        with fpath.open() as f:
            code = compile(f.read(), fpath, 'exec')
            exec(code, globals)
    elif ext == '.sage':
        from sage.repl.attach import load_attach_mode
        from sage.repl.preparse import preparse_file_named, preparse_file
        load_debug_mode, attach_debug_mode = load_attach_mode()
        if (attach and attach_debug_mode) or ((not attach) and load_debug_mode):
            # Preparse to a file to enable tracebacks with
            # code snippets. Use preparse_file_named to make
            # the file name appear in the traceback as well.
            # See Issue 11812.
            if attach:
                add_attached_file(fpath)
            parsed_file = preparse_file_named(fpath)
            with parsed_file.open() as f:
                code = compile(f.read(), parsed_file, 'exec')
                exec(code, globals)
        else:
            # Preparse in memory only for speed.
            if attach:
                add_attached_file(fpath)
            with fpath.open() as f:
                exec(preparse_file(f.read()) + "\n", globals)
    elif ext in ['.spyx', '.pyx']:
        if attach:
            add_attached_file(fpath)
        exec(load_cython(fpath), globals)
    elif ext in ['.f', '.f90']:
        from sage.misc.inline_fortran import fortran
        with fpath.open() as f:
            fortran(f.read(), globals)
    elif ext == '.m':
        # Assume magma for now, though maybe .m is used by maple and
        # mathematica too, and we should really analyze the file
        # further.
        s = globals['magma'].load(fpath)
        i = s.find('\n')
        print(s[i + 1:])
    else:
        raise ValueError('unknown file extension %r for load or attach (supported extensions: .py, .pyx, .sage, .spyx, .f, .f90, .m)' % ext)


def load_wrap(filename, attach=False):
    """
    Encode a load or attach command as valid Python code.

    INPUT:

    - ``filename`` -- string or :class:`Path` object; the argument
      to the load or attach command

    - ``attach`` -- boolean (default: ``False``); whether to attach
      ``filename``, instead of loading it

    OUTPUT: string

    EXAMPLES::

        sage: sage.repl.load.load_wrap('foo.py', True)
        'sage.repl.load.load(sage.repl.load.base64.b64decode("Zm9vLnB5"),globals(),True)'
        sage: sage.repl.load.load_wrap('foo.sage')
        'sage.repl.load.load(sage.repl.load.base64.b64decode("Zm9vLnNhZ2U="),globals(),False)'
        sage: m = sage.repl.load.base64.b64decode("Zm9vLnNhZ2U=")
        sage: m == b'foo.sage'
        True
    """
    if isinstance(filename, Path):
        filename = str(filename)
    # Note: In Python 3, b64encode only accepts bytes, and returns bytes.
    b64 = base64.b64encode(str_to_bytes(filename, FS_ENCODING,
                                        "surrogateescape"))
    txt = 'sage.repl.load.load(sage.repl.load.base64.b64decode("{}"),globals(),{})'
    return txt.format(bytes_to_str(b64, 'ascii'), attach)
