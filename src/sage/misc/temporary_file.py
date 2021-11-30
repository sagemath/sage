# -*- coding: utf-8 -*-
"""
Temporary file handling

AUTHORS:

- Volker Braun, Jeroen Demeyer (2012-10-18): move these functions here
  from sage/misc/misc.py and make them secure, see :trac:`13579`.

- Jeroen Demeyer (2013-03-17): add :class:`atomic_write`,
  see :trac:`14292`.

- Sebastian Oehms (2021-08-07): add :class:`atomic_dir`,
  see :trac:`32344`
"""
# ****************************************************************************
#       Copyright (C) 2012 Volker Braun <vbraun@stp.dias.ie>
#       Copyright (C) 2012 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import io
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
        sage: f = open('file_inside_d', 'w')

    Temporary directories are unaccessible by other users::

        sage: os.stat(d).st_mode & 0o077
        0
        sage: f.close()
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
        sage: f = open(fn, 'w')

    Temporary files are unaccessible by other users::

        sage: os.stat(fn).st_mode & 0o077
        0
        sage: f.close()
    """
    from sage.misc.misc import SAGE_TMP
    handle, tmp = tempfile.mkstemp(prefix=name, suffix=ext, dir=str(SAGE_TMP))
    os.close(handle)
    name = os.path.abspath(tmp)
    return name


#################################################################
# write to a temporary file and move it in place
#################################################################
class atomic_write(object):
    """
    Write to a given file using a temporary file and then rename it
    to the target file. This renaming should be atomic on modern
    operating systems. Therefore, this class can be used to avoid race
    conditions when a file might be read while it is being written.
    It also avoids having partially written files due to exceptions
    or crashes.

    This is to be used in a ``with`` statement, where a temporary file
    is created when entering the ``with`` and is moved in place of the
    target file when exiting the ``with`` (if no exceptions occurred).

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
      mode bits of the file were changed manually). (Not to be confused with
      the file opening mode.)

    - ``binary`` -- (boolean, default: True on Python 2, False on Python 3) the
      underlying file is opened in binary mode.  If False then it is opened in
      text mode and an encoding with which to write the file may be supplied.

    - ``**kwargs`` -- additional keyword arguments passed to the underlying
      `io.open` call.

    EXAMPLES::

        sage: from sage.misc.temporary_file import atomic_write
        sage: target_file = tmp_filename()
        sage: with open(target_file, 'w') as f:
        ....:     _ = f.write("Old contents")
        sage: with atomic_write(target_file) as f:
        ....:     _ = f.write("New contents")
        ....:     f.flush()
        ....:     with open(target_file, 'r') as f2:
        ....:         f2.read()
        'Old contents'
        sage: with open(target_file, 'r') as f:
        ....:     f.read()
        'New contents'

    The name of the temporary file can be accessed using ``f.name``.
    It is not a problem to close and re-open the temporary file::

        sage: from sage.misc.temporary_file import atomic_write
        sage: target_file = tmp_filename()
        sage: with open(target_file, 'w') as f:
        ....:     _ = f.write("Old contents")
        sage: with atomic_write(target_file) as f:
        ....:     f.close()
        ....:     with open(f.name, 'w') as f2:
        ....:         _ = f2.write("Newer contents")
        sage: with open(target_file, 'r') as f:
        ....:     f.read()
        'Newer contents'

    If an exception occurs while writing the file, the target file is
    not touched::

        sage: with atomic_write(target_file) as f:
        ....:     _ = f.write("Newest contents")
        ....:     raise RuntimeError
        Traceback (most recent call last):
        ...
        RuntimeError
        sage: with open(target_file, 'r') as f:
        ....:     f.read()
        'Newer contents'

    Some examples of using the ``append`` option. Note that the file
    is never opened in "append" mode, it is possible to overwrite
    existing data::

        sage: target_file = tmp_filename()
        sage: with atomic_write(target_file, append=True) as f:
        ....:     _ = f.write("Hello")
        sage: with atomic_write(target_file, append=True) as f:
        ....:     _ = f.write(" World")
        sage: with open(target_file, 'r') as f:
        ....:     f.read()
        'Hello World'
        sage: with atomic_write(target_file, append=True) as f:
        ....:     _ = f.seek(0)
        ....:     _ = f.write("HELLO")
        sage: with open(target_file, 'r') as f:
        ....:     f.read()
        'HELLO World'

    If the target file is a symbolic link, the link is kept and the
    target of the link is written to::

        sage: link_to_target = os.path.join(tmp_dir(), "templink")
        sage: os.symlink(target_file, link_to_target)
        sage: with atomic_write(link_to_target) as f:
        ....:     _ = f.write("Newest contents")
        sage: with open(target_file, 'r') as f:
        ....:     f.read()
        'Newest contents'

    We check the permission bits of the new file. Note that the old
    permissions do not matter::

        sage: os.chmod(target_file, 0o600)
        sage: _ = os.umask(0o022)
        sage: with atomic_write(target_file) as f:
        ....:     pass
        sage: '{:#o}'.format(os.stat(target_file).st_mode & 0o777)
        '0o644'
        sage: _ = os.umask(0o077)
        sage: with atomic_write(target_file, mode=0o777) as f:
        ....:     pass
        sage: '{:#o}'.format(os.stat(target_file).st_mode & 0o777)
        '0o700'

    Test writing twice to the same target file. The outermost ``with``
    "wins"::

        sage: with open(target_file, 'w') as f:
        ....:     _ = f.write('>>> ')
        sage: with atomic_write(target_file, append=True) as f, \
        ....:          atomic_write(target_file, append=True) as g:
        ....:     _ = f.write("AAA"); f.close()
        ....:     _ = g.write("BBB"); g.close()
        sage: with open(target_file, 'r') as f:
        ....:     f.read()
        '>>> AAA'

    Supplying an encoding means we're writing the file in "text mode" (in the
    same sense as `io.open`) and so we must write unicode strings::

        sage: target_file = tmp_filename()
        sage: with atomic_write(target_file, binary=False,
        ....:                   encoding='utf-8') as f:
        ....:     _ = f.write(u'Hélas')
        sage: import io
        sage: with io.open(target_file, encoding='utf-8') as f:
        ....:     print(f.read())
        Hélas

    Supplying an encoding in binary mode (or other arguments that don't
    make sense to `io.open` in binary mode) is an error::

        sage: writer = atomic_write(target_file, binary=True,
        ....:                       encoding='utf-8')
        sage: with writer as f:
        ....:     _ = f.write(u'Hello')
        Traceback (most recent call last):
        ...
        ValueError: binary mode doesn't take an encoding argument
        sage: os.path.exists(writer.tempname)
        False
    """
    def __init__(self, target_filename, append=False, mode=0o666,
                 binary=None, **kwargs):
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
        umask = os.umask(0)
        os.umask(umask)
        self.mode = mode & (~umask)

        # 'binary' mode is the default on Python 2, whereas 'text' mode is the
        # default on Python 3--this reflects consistent handling of the default
        # str type on the two platforms
        self.binary = False if binary is None else binary
        self.kwargs = kwargs

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

        fd, name = tempfile.mkstemp(dir=self.tmpdir)
        self.tempname = os.path.abspath(name)

        rmode = 'r' + ('b' if self.binary else '')
        wmode = 'w+' + ('b' if self.binary else '')

        try:
            self.tempfile = io.open(name, wmode, **self.kwargs)
        except Exception:
            # Some invalid arguments were passed to io.open
            os.unlink(name)
            raise
        finally:
            os.close(fd)

        os.chmod(name, self.mode)
        if self.append:
            try:
                with io.open(self.target, rmode, **self.kwargs) as f:
                    r = f.read()
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

#################################################################
# write to a temporary directory and move it in place
#################################################################
class atomic_dir(object):
    """
    Write to a given directory using a temporary directory and then rename it
    to the target directory. This is for creating a directory whose contents
    are determined uniquely by the directory name. If multiple threads or
    processes attempt to create it in parallel, then it does not matter which
    thread created it. Despite this assumption the contents of the directories
    differ in the examples for demonstration purpose.

    See also :class:`atomic_write`.

    INPUT:

    - ``target_directory`` -- the name of the directory to be written.
      If it exists then the previous contents will be kept.

    EXAMPLES::

        sage: from sage.misc.temporary_file import atomic_dir
        sage: target_dir = tmp_dir()
        sage: with atomic_dir(target_dir) as d:
        ....:     target_file = os.path.join(d.name, 'test')
        ....:     with open(target_file, 'w') as f:
        ....:        _ = f.write("First")
        ....:        f.flush()
        ....:     with atomic_dir(target_dir) as e:
        ....:         target_file2 = os.path.join(e.name, 'test')
        ....:         with open(target_file2, 'w') as g:
        ....:            _ = g.write("Second")
        ....:            g.flush()
        ....:     with open(target_file, 'r') as f:
        ....:         f.read()
        'First'
        sage: with atomic_dir(target_dir) as d:
        ....:     target_file = os.path.join(d.name, 'test')
        ....:     with open(target_file, 'w') as f:
        ....:        _ = f.write("Third")
        sage: target = os.path.join(target_dir, 'test')
        sage: with open(target, 'r') as h:
        ....:     h.read()
        'Second'
    """
    def __init__(self, target_directory):
        r"""
        TESTS::

            sage: from sage.misc.temporary_file import atomic_dir
            sage: link_to_target = os.path.join(tmp_dir(), "templink")
            sage: os.symlink("/foobar", link_to_target)
            sage: aw = atomic_dir(link_to_target)
            sage: print(aw.target)
            /foobar
            sage: print(aw.tmpdir)
            /
        """
        self.target = os.path.realpath(target_directory)
        self.tmpdir = os.path.dirname(self.target)

    def __enter__(self):
        r"""
        Create and return a temporary directory in ``self.tmpdir`` (normally
        the same directory as the target file).

        OUTPUT: a directory returned by :func:`tempfile.TemporaryDirectory`.

        TESTS::

            sage: from sage.misc.temporary_file import atomic_dir
            sage: aw = atomic_dir(tmp_dir())
            sage: with aw as d:
            ....:     os.path.dirname(aw.target) == os.path.dirname(d.name)
            True
        """
        tdir = tempfile.TemporaryDirectory(dir=self.tmpdir)
        self.tempname = os.path.abspath(tdir.name)
        return tdir

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        If the ``with`` block was successful, move the temporary directory
        to the target directory. Otherwise, delete the temporary directory.

        TESTS:

        Check that the temporary directory is deleted if there was an
        exception::

            sage: from sage.misc.temporary_file import atomic_dir
            sage: with atomic_dir(tmp_dir()) as d:
            ....:     tempname = d.name
            ....:     raise RuntimeError
            Traceback (most recent call last):
            ...
            RuntimeError
            sage: os.path.exists(tempname)
            False
        """
        import shutil
        if exc_type is None:
            # Success: move temporary file to target file
            try:
                os.rename(self.tempname, self.target)
            except OSError:
                # Race: Another thread or process must have created the directory
                pass
        else:
            # Failure: delete temporary file
            shutil.rmtree(self.tempname)
