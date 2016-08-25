"""
Temporary file handling

AUTHORS:

- Volker Braun, Jeroen Demeyer (2012-10-18): move these functions here
  from sage/misc/misc.py and make them secure, see :trac:`13579`.

- Jeroen Demeyer (2013-03-17): add :class:`atomic_write`,
  see :trac:`14292`.
"""

#*****************************************************************************
#       Copyright (C) 2012 Volker Braun <vbraun@stp.dias.ie>
#       Copyright (C) 2012 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

import os
import tempfile
import atexit

def delete_tmpfiles():
    """
    Remove the directory ``SAGE_TMP``.

    TESTS:

    This is automatically run when Sage exits, test this by running a
    separate session of Sage::

        sage: from sage.tests.cmdline import test_executable
        sage: child_SAGE_TMP, err, ret = test_executable(["sage", "-c", "print(SAGE_TMP)"])
        sage: err, ret
        ('', 0)
        sage: os.path.exists(child_SAGE_TMP)  # indirect doctest
        False

    The parent directory should exist::

        sage: parent_SAGE_TMP = os.path.normpath(child_SAGE_TMP + '/..')
        sage: os.path.isdir(parent_SAGE_TMP)
        True
    """
    import shutil
    from sage.misc.misc import SAGE_TMP
    shutil.rmtree(str(SAGE_TMP), ignore_errors=True)


# Run when Python shuts down
atexit.register(delete_tmpfiles)


#################################################################
# temporary directory
#################################################################

def tmp_dir(name="dir_", ext=""):
    r"""
    Create and return a temporary directory in
    ``$HOME/.sage/temp/hostname/pid/``

    The temporary directory is deleted automatically when Sage exits.

    INPUT:

    - ``name`` -- (default: ``"dir_"``) A prefix for the directory name.

    - ``ext`` -- (default: ``""``) A suffix for the directory name.

    OUTPUT:

    The absolute path of the temporary directory created, with a
    trailing slash (or whatever the path separator is on your OS).

    EXAMPLES::

        sage: d = tmp_dir('dir_testing_', '.extension')
        sage: d   # random output
        '/home/username/.sage/temp/hostname/7961/dir_testing_XgRu4p.extension/'
        sage: os.chdir(d)
        sage: _ = open('file_inside_d', 'w')

    Temporary directories are unaccessible by other users::

        sage: os.stat(d).st_mode & 0o077
        0
    """
    from sage.misc.misc import SAGE_TMP
    tmp = tempfile.mkdtemp(prefix=name, suffix=ext, dir=str(SAGE_TMP))
    name = os.path.abspath(tmp)
    return name + os.sep


#################################################################
# temporary filename
#################################################################

def tmp_filename(name="tmp_", ext=""):
    r"""
    Create and return a temporary file in
    ``$HOME/.sage/temp/hostname/pid/``

    The temporary file is deleted automatically when Sage exits.

    .. warning::

        If you need a particular file extension always use
        ``tmp_filename(ext=".foo")``, this will ensure that the file
        does not yet exist. If you were to use
        ``tmp_filename()+".foo"``, then you might overwrite an
        existing file!

    INPUT:

    - ``name`` -- (default: ``"tmp_"``) A prefix for the file name.

    - ``ext`` -- (default: ``""``) A suffix for the file name. If you
      want a filename extension in the usual sense, this should start
      with a dot.

    OUTPUT:

    The absolute path of the temporary file created.

    EXAMPLES::

        sage: fn = tmp_filename('just_for_testing_', '.extension')
        sage: fn  # random
        '/home/username/.sage/temp/hostname/8044/just_for_testing_tVVHsn.extension'
        sage: _ = open(fn, 'w')

    Temporary files are unaccessible by other users::

        sage: os.stat(fn).st_mode & 0o077
        0
    """
    from sage.misc.misc import SAGE_TMP
    handle, tmp = tempfile.mkstemp(prefix=name, suffix=ext, dir=str(SAGE_TMP))
    os.close(handle)
    name = os.path.abspath(tmp)
    return name


def graphics_filename(ext='.png'):
    """
    Deprecated SageNB graphics filename

    You should just use :meth:`tmp_filename`.

    When run from the Sage notebook, return the next available canonical
    filename for a plot/graphics file in the current working directory.
    Otherwise, return a temporary file inside ``SAGE_TMP``.

    INPUT:

    - ``ext`` -- (default: ``".png"``) A file extension (including the dot)
      for the filename.

    OUTPUT:

    The path of the temporary file created. In the notebook, this is
    a filename without path in the current directory. Otherwise, this
    an absolute path.

    EXAMPLES::

        sage: from sage.misc.temporary_file import graphics_filename
        sage: print(graphics_filename())  # random, typical filename for sagenb
        sage0.png

    TESTS:

    When doctesting, this returns instead a random temporary file.
    We check that it's a file inside ``SAGE_TMP`` and that the extension
    is correct::

        sage: fn = graphics_filename(ext=".jpeg")
        sage: fn.startswith(str(SAGE_TMP))
        True
        sage: fn.endswith('.jpeg')
        True

    Historically, it was also possible to omit the dot. This has been
    changed in :trac:`16640` but it will still work for now::

        sage: fn = graphics_filename("jpeg")
        doctest:...: DeprecationWarning: extension must now include the dot
        See http://trac.sagemath.org/16640 for details.
        sage: fn.endswith('.jpeg')
        True
    """
    import sage.plot.plot
    from sage.misc.superseded import deprecation
    if ext[0] not in '.-':
        deprecation(16640, "extension must now include the dot")
        ext = '.' + ext
    if sage.plot.plot.EMBEDDED_MODE:
        # Don't use this unsafe function except in the notebook, #15515
        i = 0
        while os.path.exists('sage%d%s'%(i,ext)):
            i += 1
        filename = 'sage%d%s'%(i,ext)
        return filename
    else:
        deprecation(17234,'use tmp_filename instead')
        return tmp_filename(ext=ext)


#################################################################
# write to a temporary file and move it in place
#################################################################
class atomic_write:
    """
    Write to a given file using a temporary file and then rename it
    to the target file. This renaming should be atomic on modern
    operating systems. Therefore, this class can be used to avoid race
    conditions when a file might be read while it is being written.
    It also avoids having partially written files due to exceptions
    or crashes.

    This is to be used in a ``with`` statement, where a temporary file
    is created when entering the ``with`` and is moved in place of the
    target file when exiting the ``with`` (if no exceptions occured).

    INPUT:

    - ``target_filename`` -- the name of the file to be written.
      Normally, the contents of this file will be overwritten.

    - ``append`` -- (boolean, default: False) if True and
      ``target_filename`` is an existing file, then copy the current
      contents of ``target_filename`` to the temporary file when
      entering the ``with`` statement. Otherwise, the temporary file is
      initially empty.

    - ``mode`` -- (default: ``0o666``) mode bits for the file. The
      temporary file is created with mode ``mode & ~umask`` and the
      resulting file will also have these permissions (unless the
      mode bits of the file were changed manually).

    EXAMPLES::

        sage: from sage.misc.temporary_file import atomic_write
        sage: target_file = tmp_filename()
        sage: open(target_file, "w").write("Old contents")
        sage: with atomic_write(target_file) as f:
        ....:     f.write("New contents")
        ....:     f.flush()
        ....:     open(target_file, "r").read()
        'Old contents'
        sage: open(target_file, "r").read()
        'New contents'

    The name of the temporary file can be accessed using ``f.name``.
    It is not a problem to close and re-open the temporary file::

        sage: from sage.misc.temporary_file import atomic_write
        sage: target_file = tmp_filename()
        sage: open(target_file, "w").write("Old contents")
        sage: with atomic_write(target_file) as f:
        ....:     f.close()
        ....:     open(f.name, "w").write("Newer contents")
        sage: open(target_file, "r").read()
        'Newer contents'

    If an exception occurs while writing the file, the target file is
    not touched::

        sage: with atomic_write(target_file) as f:
        ....:     f.write("Newest contents")
        ....:     raise RuntimeError
        Traceback (most recent call last):
        ...
        RuntimeError
        sage: open(target_file, "r").read()
        'Newer contents'

    Some examples of using the ``append`` option. Note that the file
    is never opened in "append" mode, it is possible to overwrite
    existing data::

        sage: target_file = tmp_filename()
        sage: with atomic_write(target_file, append=True) as f:
        ....:     f.write("Hello")
        sage: with atomic_write(target_file, append=True) as f:
        ....:     f.write(" World")
        sage: open(target_file, "r").read()
        'Hello World'
        sage: with atomic_write(target_file, append=True) as f:
        ....:     f.seek(0)
        ....:     f.write("HELLO")
        sage: open(target_file, "r").read()
        'HELLO World'

    If the target file is a symbolic link, the link is kept and the
    target of the link is written to::

        sage: link_to_target = os.path.join(tmp_dir(), "templink")
        sage: os.symlink(target_file, link_to_target)
        sage: with atomic_write(link_to_target) as f:
        ....:     f.write("Newest contents")
        sage: open(target_file, "r").read()
        'Newest contents'

    We check the permission bits of the new file. Note that the old
    permissions do not matter::

        sage: os.chmod(target_file, 0o600)
        sage: _ = os.umask(0o022)
        sage: with atomic_write(target_file) as f:
        ....:     pass
        sage: oct(os.stat(target_file).st_mode & 0o777)
        '644'
        sage: _ = os.umask(0o077)
        sage: with atomic_write(target_file, mode=0o777) as f:
        ....:     pass
        sage: oct(os.stat(target_file).st_mode & 0o777)
        '700'

    Test writing twice to the same target file. The outermost ``with``
    "wins"::

        sage: open(target_file, "w").write(">>> ")
        sage: with atomic_write(target_file, append=True) as f, \
        ....:          atomic_write(target_file, append=True) as g:
        ....:     f.write("AAA"); f.close()
        ....:     g.write("BBB"); g.close()
        sage: open(target_file, "r").read()
        '>>> AAA'
    """
    def __init__(self, target_filename, append=False, mode=0o666):
        """
        TESTS::

            sage: from sage.misc.temporary_file import atomic_write
            sage: link_to_target = os.path.join(tmp_dir(), "templink")
            sage: os.symlink("/foobar", link_to_target)
            sage: aw = atomic_write(link_to_target)
            sage: print(aw.target)
            /foobar
            sage: print(aw.tmpdir)
            /
        """
        self.target = os.path.realpath(target_filename)
        self.tmpdir = os.path.dirname(self.target)
        self.append = append
        # Remove umask bits from mode
        umask = os.umask(0); os.umask(umask)
        self.mode = mode & (~umask)

    def __enter__(self):
        """
        Create and return a temporary file in ``self.tmpdir`` (normally
        the same directory as the target file).

        If ``self.append``, then copy the current contents of
        ``self.target`` to the temporary file.

        OUTPUT: a file returned by :func:`tempfile.NamedTemporaryFile`.

        TESTS::

            sage: from sage.misc.temporary_file import atomic_write
            sage: aw = atomic_write(tmp_filename())
            sage: with aw as f:
            ....:     os.path.dirname(aw.target) == os.path.dirname(f.name)
            True
        """
        self.tempfile = tempfile.NamedTemporaryFile(dir=self.tmpdir, delete=False)
        self.tempname = self.tempfile.name
        os.chmod(self.tempname, self.mode)
        if self.append:
            try:
                r = open(self.target).read()
            except IOError:
                pass
            else:
                self.tempfile.write(r)
        return self.tempfile

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        If the ``with`` block was successful, move the temporary file
        to the target file. Otherwise, delete the temporary file.

        TESTS:

        Check that the temporary file is deleted if there was an
        exception::

            sage: from sage.misc.temporary_file import atomic_write
            sage: with atomic_write(tmp_filename()) as f:
            ....:     tempname = f.name
            ....:     raise RuntimeError
            Traceback (most recent call last):
            ...
            RuntimeError
            sage: os.path.exists(tempname)
            False
        """
        # Flush the file contents to disk (to be safe even if the
        # system crashes) and close the file.
        if not self.tempfile.closed:
            self.tempfile.flush()
            os.fsync(self.tempfile.fileno())
            self.tempfile.close()

        if exc_type is None:
            # Success: move temporary file to target file
            try:
                os.rename(self.tempname, self.target)
            except OSError:
                os.unlink(self.target)
                os.rename(self.tempname, self.target)
        else:
            # Failure: delete temporary file
            os.unlink(self.tempname)
