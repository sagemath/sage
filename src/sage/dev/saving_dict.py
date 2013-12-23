r"""
A dictionary which is automatically stored in the file system

This module provides a special dictionary class for :class:`sagedev.SageDev`
which is automatically written to the file system whenever it changes.

AUTHORS:

- David Roe, Julian Rueth, Robert Bradshaw: initial version

"""
#*****************************************************************************
#       Copyright (C) 2013 David Roe <roed.math@gmail.com>
#                          Julian Rueth <julian.rueth@fsfe.org>
#                          Robert Bradshaw <robertwb@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import cPickle
from cStringIO import StringIO
import os
import collections

def _raise():
    r"""
    Helper method for :class:`SavingDict` which reraises an exception.

    TESTS::

        sage: from sage.dev.saving_dict import _raise
        sage: try:
        ....:     raise Exception("this is a test")
        ....: except Exception:
        ....:     s = "the exception was caught"
        ....:     _raise()
        Traceback (most recent call last):
        ...
        Exception: this is a test
        sage: print s
        the exception was caught
    """
    raise

class SavingDict(collections.MutableMapping):
    r"""
    Dictionary-like class that saves itself on the filesystem.

    INPUT:

    - ``filename`` -- a string, file to store SavingDict

    - ``values`` -- dictionary-like object that sets initial values or ``None``
      (default: ``None``); by default initial values will be loaded from
      filename, if it exists, otherwise the initial values will be empty.

    - ``default`` -- callabable object requiring no arguments that provides a
      default value for keys not in the :class:`SavingDict` (default:
      ``_raise``, which raises an exception if a key is not present)

    - ``paired`` -- another :class:`SavingDict` or ``None`` (default:
      ``None``), this dictionary will be paired as per :meth:`set_paired`

    EXAMPLES::

        sage: from sage.dev.saving_dict import SavingDict
        sage: tmp = tmp_filename()
        sage: sd = SavingDict(tmp)
        sage: sd['cat'] = 'meow'; sd['cat']
        'meow'
        sage: sd['cow'] = 'moo'; sd
        {'cow': 'moo', 'cat': 'meow'}
        sage: del sd; sd
        Traceback (most recent call last):
        ...
        NameError: name 'sd' is not defined

    Creating another instance (with the same file name) loads the
    cached values::

        sage: sd = SavingDict(tmp); sd
        {'cow': 'moo', 'cat': 'meow'}
        sage: sd._erase()
    """
    def __init__(self,
            filename,
            values=None,
            default=_raise,
            paired=None):
        """
        Initialization.

        EXAMPLES::

            sage: from sage.dev.saving_dict import SavingDict
            sage: sd = SavingDict(tmp_filename())
            sage: sd['cat'] = 'meow'; sd
            {'cat': 'meow'}
            sage: sd._erase()
        """
        if not callable(default):
            raise ValueError("default must be callable")

        self._filename = filename
        if values is None:
            self._dict = SavingDict.load_dict_from_file(filename)
        else:
            self._dict = dict(values) # explicitly make copy
        self._default = default
        self._pairing = None
        if paired is not None:
            self.set_paired(paired)
        if not os.path.exists(self._filename):
            self._write()
            
    def __repr__(self):
        r"""
        Return a printable representation of this object.

        TESTS::

            sage: from sage.dev.saving_dict import SavingDict
            sage: sd = SavingDict(tmp_filename(), {0:1, 1:2})
            sage: repr(sd)
            '{0: 1, 1: 2}'
            sage: sd._erase()
        """
        return repr(self._dict)

    def unset_pairing(self):
        r"""
        Unset any pairing that was constructed with :meth:`set_paired`.

        EXAMPLES::

            sage: from sage.dev.saving_dict import SavingDict
            sage: sd1 = SavingDict(tmp_filename())
            sage: sd2 = SavingDict(tmp_filename(), paired=sd1)
            sage: sd1[0] = 1; sd1[0]; sd2[1]
            1
            0
            sage: sd1.unset_pairing()
            sage: sd1._erase()
            sage: del sd1[0]; sd1[0]
            Traceback (most recent call last):
            ...
            KeyError: 0
            sage: sd2[1]
            0
            sage: sd2._erase()
        """
        if self._pairing:
            with self._pairing as paired:
                paired.unset_pairing()
        self._pairing = None

    def set_paired(self, other):
        r"""
        Set another class:`SavingDict` to be updated with the reverse of this
        one and vice versa.

        EXAMPLES::

            sage: from sage.dev.saving_dict import SavingDict
            sage: sd1 = SavingDict(tmp_filename())
            sage: sd2 = SavingDict(tmp_filename())
            sage: sd1.set_paired(sd2)
            sage: sd1[0] = 1; sd1[0]; sd2[1]
            1
            0
            sage: del sd1[0]
            sage: sd1[0]
            Traceback (most recent call last):
            ...
            KeyError: 0
            sage: sd2[1]
            Traceback (most recent call last):
            ...
            KeyError: 1
            sage: sd2[2] = 3; sd1[3]; sd2[2]
            2
            3
            sage: sd2[2] = 4; sd1[3]
            Traceback (most recent call last):
            ...
            KeyError: 3
            sage: sd1[4]; sd2[2]
            2
            4
            sage: sd1._erase(); sd2._erase()
        """
        if not isinstance(other, SavingDict):
            raise ValueError("other is not a SavingDict")

        self.unset_pairing()
        other.unset_pairing()

        class Pairing(object):
            def __init__(this, obj):
                this._obj = obj
            def __enter__(this):
                this._entered[0] += 1
                return other if this._obj is self else self
            def __exit__(this, type, value, tb):
                this._entered[0] -= 1
            def __nonzero__(this):
                return not this._entered[0]

        Pairing._entered = [0] # ints are immutable

        self._pairing, other._pairing = Pairing(self), Pairing(other)

    def _write(self):
        r"""
        Write this dictionary to disk.

        EXAMPLES::

            sage: from sage.dev.saving_dict import SavingDict
            sage: tmp = tmp_filename()
            sage: sd = SavingDict(tmp, {0:1, 1:2})
            sage: SavingDict(tmp)
            {}
            sage: sd._write()
            sage: SavingDict(tmp)
            {0: 1, 1: 2}
            sage: sd._erase()
        """
        from sage.doctest import DOCTEST_MODE
        if DOCTEST_MODE:
            from sage.misc.misc import SAGE_TMP
            SAGE_TMP = str(SAGE_TMP)
            error = "write attempt to a saving_dict in a doctest: "+self._filename
            assert os.path.abspath(self._filename).startswith(SAGE_TMP), error

        import tempfile
        dirname = os.path.dirname(os.path.abspath(self._filename))
        fd, tmpfile = tempfile.mkstemp(dir=dirname)
        s = cPickle.dumps(self._dict, protocol=2)
        with os.fdopen(fd, "wb") as F:
            F.write(s)
        try:
            # This move is atomic (the files are on the same filesystem)
            os.rename(tmpfile, self._filename)
        except (IOError, OSError):
            # Lesser operation systems cannot do atomic moves (looking at you, windows)
            os.unlink(self._filename)
            os.rename(tmpfile, self._filename)

    def _erase(self):
        """
        Erase the backing file.
        
        EXAMPLES::

            sage: from sage.dev.saving_dict import SavingDict
            sage: tmp = tmp_filename()
            sage: sd = SavingDict(tmp, {0:1, 1:2})
            sage: sd._write()
            sage: os.path.exists(tmp)
            True
            sage: sd._erase()
            sage: os.path.exists(tmp)
            False
        """
        os.unlink(self._filename)

    def __setitem__(self, key, value):
        r"""
        Set ``key`` to ``value``.

        TESTS::

            sage: from sage.dev.saving_dict import SavingDict
            sage: sd = SavingDict(tmp_filename())
            sage: sd['cow'] = 'moo'
            sage: sd
            {'cow': 'moo'}
            sage: sd._erase()
        """
        if self._pairing:
            with self._pairing as paired:
                if key in self:
                    del paired[self[key]]
                paired[value] = key
        self._dict[key] = value
        self._write()

    def __delitem__(self, key):
        r"""
        Remove ``key`` from this dictionary.

        TESTS::

            sage: from sage.dev.saving_dict import SavingDict
            sage: sd = SavingDict(tmp_filename(), {'cow': 'moo'}); sd
            {'cow': 'moo'}
            sage: del sd['cow']
            sage: del sd['cow']
            Traceback (most recent call last):
            ...
            KeyError: 'cow'
            sage: sd
            {}
            sage: sd._erase()
        """
        if self._pairing and key in self:
            with self._pairing as paired:
                del paired[self[key]]
        del self._dict[key]
        self._write()

    def __getitem__(self, key):
        r"""
        Return the value for ``key``.

        TESTS::

            sage: from sage.dev.saving_dict import SavingDict
            sage: sd = SavingDict(tmp_filename(), {'cow': 'moo'}); sd
            {'cow': 'moo'}
            sage: sd['cow']
            'moo'
            sage: sd['moo']
            Traceback (most recent call last):
            ...
            KeyError: 'moo'
            sage: sd._erase()
        """
        try:
            return self._dict[key]
        except KeyError:
            return self._default()

    def __contains__(self, key):
        r"""
        Return whether this dictionary contains ``key``.

        TESTS::

            sage: from sage.dev.saving_dict import SavingDict
            sage: sd = SavingDict(tmp_filename(), {'cow': 'moo'}); sd
            {'cow': 'moo'}
            sage: 'cow' in sd
            True
            sage: 'moo' in sd
            False
            sage: sd._erase()
        """
        return key in self._dict

    def __len__(self):
        r"""
        Return the number of keys in this dictionary.

        TESTS::

            sage: from sage.dev.saving_dict import SavingDict
            sage: sd = SavingDict(tmp_filename()); len(sd)
            0
            sage: sd['cow'] = 'moo'
            sage: len(sd)
            1
            sage: sd._erase()
        """
        return len(self._dict)

    def __iter__(self):
        r"""
        Return an iterator over the keys of this dictionary.

        TESTS::

            sage: from sage.dev.saving_dict import SavingDict
            sage: sd = SavingDict(tmp_filename(), {'cow':'moo', 0:1}); sd
            {0: 1, 'cow': 'moo'}
            sage: for key in sd:
            ...       print key, sd[key]
            0 1
            cow moo
            sage: sd._erase()
        """
        return iter(self._dict)

    @classmethod
    def load_dict_from_file(cls, filename):
        r"""
        Load a pickled dictionary from ``filename``, defaults to {} if the file
        does not exist.

        TESTS::

            sage: from sage.dev.saving_dict import SavingDict
            sage: d = SavingDict.load_dict_from_file(''); d
            {}
            sage: d['cow'] = 'moo'
            sage: import cPickle
            sage: tmp = tmp_filename()
            sage: with open(tmp, 'w') as f:
            ....:     f.write(cPickle.dumps(d, protocol=2))
            sage: SavingDict.load_dict_from_file(tmp)
            {'cow': 'moo'}
            sage: os.unlink(tmp)
        """
        if os.path.exists(filename):
            with open(filename) as F:
                s = F.read()
            if s:
                unpickler = cPickle.Unpickler(StringIO(s))
                try:
                    return unpickler.load()
                except:
                    # catch-all exception! Unpickling can cause all
                    # kinds of exceptions, e.g. AttributeError if the
                    # Sage source code changed.
                    pass
        return {}
