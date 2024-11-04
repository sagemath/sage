# sage_setup: distribution = sagemath-objects
# cython: old_style_globals=True
# The old_style_globals directive is important for load() to work correctly.
# However, this should be removed in favor of user_globals; see
# https://github.com/sagemath/sage/issues/18083

r"""
Object persistence

You can load and save most Sage object to disk using the load and
save member functions and commands.

.. NOTE::

   It is impossible to save certain Sage objects to disk. For example,
   if `x` is a MAGMA object, i.e., a wrapper around an object
   that is defined in MAGMA, there is no way to save `x` to
   disk, since MAGMA doesn't support saving of individual objects to
   disk.


-  Versions: Loading and saving of objects is guaranteed to work
   even if the version of Python changes. Saved objects can be loaded
   in future versions of Python. However, if the data structure that
   defines the object, e.g., in Sage code, changes drastically (or
   changes name or disappears), then the object might not load
   correctly or work correctly.

-  Objects are zlib compressed for space efficiency.
"""

import io
import os
import pickle
import sys

from textwrap import dedent

# change to import zlib to use zlib instead; but this
# slows down loading any data stored in the other format
import zlib
import bz2
comp = zlib
comp_other = bz2

from sage.misc.sage_unittest import TestSuite


# We define two global dictionaries `already_pickled` and
# `already_unpickled`, which are intended to help you to implement
# pickling when cyclic definitions could happen.
#
# You have the guarantee that the dictionary `already_pickled`
# (resp. `already_unpickled`) will be cleared after the complete
# pickling (resp. unpickling) process has been completed (and not
# before!).
# Apart from this, you are free to use these variables as you like.
#
# However, the standard utilisation is the following.
# The pickling method (namely `__reduce__`) checks if the id of the
# current element appears in the dictionary `already_pickled`. If it
# does not, the methods records that this element is about to be
# pickled by adding the entry { id: True } to `already_pickled`.
# In all cases, the pickling method pickles the element and includes
# the id in the tuple of arguments that will be passed in afterwards
# to the unpickling function.
# The unpickling function then receives all the necessary information
# to reconstruct the element, together with an id.
# If this id appears as a key of the dictionary `already_unpickled`, it
# returns the corresponding value.
# Otherwise, it builds the element, adds the entry { id: element } to
# `already_unpickled` and finally returns the element.
#
# For a working example, see sage.rings.padics.lazy_template.LazyElement_unknown
already_pickled = { }
already_unpickled = { }


cdef _normalize_filename(s):
    """
    Append the .sobj extension to a filename if it doesn't already have it.
    """
    if s[-5:] != '.sobj':
        return s + '.sobj'

    return s


def load(*filename, compress=True, verbose=True, **kwargs):
    r"""
    Load Sage object from the file with name filename, which will have
    an ``.sobj`` extension added if it doesn't have one.  Or, if the input
    is a filename ending in ``.py``, ``.pyx``, ``.sage``, ``.spyx``,
    ``.f``, ``.f90`` or ``.m``, load that file into the current running
    session.

    Loaded files are not loaded into their own namespace, i.e., this is
    much more like Python's ``execfile`` than Python's ``import``.

    This function also loads a ``.sobj`` file over a network by
    specifying the full URL.  (Setting ``verbose = False`` suppresses
    the loading progress indicator.)

    When a pickle created with Python 2 is unpickled in Python 3, Sage
    uses the default encoding ``latin1`` to unpickle data of type :class:`str`.

    Finally, if you give multiple positional input arguments, then all
    of those files are loaded, or all of the objects are loaded and a
    list of the corresponding loaded objects is returned.

    If ``compress`` is true (the default), then the data stored in the file
    are supposed to be compressed. If ``verbose`` is true (the default), then
    some logging is printed when accessing remote files. Further keyword
    arguments are passed to :func:`pickle.load`.

    EXAMPLES::

        sage: u = 'https://www.sagemath.org/files/test.sobj'
        sage: s = load(u)                                                  # optional - internet
        Attempting to load remote file: https://www.sagemath.org/files/test.sobj
        Loading started
        Loading ended
        sage: s                                                            # optional - internet
        'hello SageMath'

    We test loading a file or multiple files or even mixing loading files and objects::

        sage: t = tmp_filename(ext='.py')
        sage: with open(t, 'w') as f:
        ....:     _ = f.write("print('hello world')")
        sage: load(t)
        hello world
        sage: load(t,t)
        hello world
        hello world
        sage: t2 = tmp_filename(); save(2/3,t2)
        sage: load(t,t,t2)
        hello world
        hello world
        [None, None, 2/3]

    Files with a ``.sage`` extension are preparsed. Also note that we
    can access global variables::

        sage: t = tmp_filename(ext='.sage')
        sage: with open(t, 'w') as f:
        ....:     _ = f.write("a += Mod(2/3, 11)")  # This evaluates to Mod(8, 11)
        sage: a = -1
        sage: load(t)
        sage: a
        7

    We can load Fortran files::

        sage: code = '      subroutine hello\n         print *, "Hello World!"\n      end subroutine hello\n'
        sage: t = tmp_filename(ext='.F')
        sage: with open(t, 'w') as f:
        ....:     _ = f.write(code)
        sage: load(t)                                                                   # needs numpy
        sage: hello                                                                     # needs numpy
        <fortran ...>
    """
    import sage.repl.load
    if len(filename) != 1:
        v = [load(file, compress=compress, verbose=verbose, **kwargs) for file in filename]
        # Return v if one of the filenames refers to an object and not
        # a loadable filename.
        for file in filename:
            if not sage.repl.load.is_loadable_filename(file):
                return v
        return

    filename = filename[0]

    if sage.repl.load.is_loadable_filename(filename):
        sage.repl.load.load(filename, globals())
        return

    ## Check if filename starts with "http://" or "https://"
    lower = filename.lower()
    if lower.startswith("http://") or lower.startswith("https://"):
        from sage.misc.remote_file import get_remote_file
        filename = get_remote_file(filename, verbose=verbose)
        tmpfile_flag = True
    else:
        tmpfile_flag = False
        filename = _normalize_filename(filename)

    ## Load file by absolute filename
    with open(filename, 'rb') as fobj:
        X = loads(fobj.read(), compress=compress, **kwargs)
    try:
        X._default_filename = os.path.abspath(filename)
    except AttributeError:
        pass

    ## Delete the tempfile, if it exists
    if tmpfile_flag:
        os.unlink(filename)

    return X


def _base_save(obj, filename, compress=True):
    """
    Base implementation for dumping an object to a ``.sobj`` file.

    This is the implementation used by :meth:`SageObject.save` and by
    :func:`save` unless the object being saved has a custom ``.save()`` method,
    in which case that is tried first.

    Otherwise this is equivalent to :func:`_base_dumps` just with the resulting
    pickle data saved to a ``.sobj`` file.
    """

    filename = _normalize_filename(filename)

    with open(filename, 'wb') as fobj:
        fobj.write(_base_dumps(obj, compress=compress))

    return filename


def save(obj, filename, compress=True, **kwargs):
    """
    Save ``obj`` to the file with name ``filename``, which will have an
    ``.sobj`` extension added if it doesn't have one and if ``obj``
    doesn't have its own ``save()`` method, like e.g. Python tuples.

    For image objects and the like (which have their own ``save()`` method),
    you may have to specify a specific extension, e.g. ``.png``, if you
    don't want the object to be saved as a Sage object (or likewise, if
    ``filename`` could be interpreted as already having some extension).

    .. WARNING::

       This will *replace* the contents of the file if it already exists.

    EXAMPLES::

        sage: import tempfile
        sage: d = tempfile.TemporaryDirectory()
        sage: a = matrix(2, [1,2, 3,-5/2])                                              # needs sage.modules
        sage: objfile = os.path.join(d.name, 'test.sobj')
        sage: objfile_short = os.path.join(d.name, 'test')
        sage: save(a, objfile)                                                          # needs sage.modules
        sage: load(objfile_short)                                                       # needs sage.modules
        [   1    2]
        [   3 -5/2]

        sage: # needs sage.plot sage.schemes
        sage: E = EllipticCurve([-1,0])
        sage: P = plot(E)
        sage: save(P, objfile_short)   # saves the plot to "test.sobj"
        sage: save(P, filename=os.path.join(d.name, "sage.png"), xmin=-2)
        sage: save(P, os.path.join(d.name, "filename.with.some.wrong.ext"))
        Traceback (most recent call last):
        ...
        ValueError: allowed file extensions for images are
        '.eps', '.pdf', '.pgf', '.png', '.ps', '.sobj', '.svg'!
        sage: print(load(objfile))
        Graphics object consisting of 2 graphics primitives

        sage: save("A python string", os.path.join(d.name, 'test'))
        sage: load(objfile)
        'A python string'
        sage: load(objfile_short)
        'A python string'
        sage: d.cleanup()

    TESTS:

    Check that :issue:`11577` is fixed::

        sage: import tempfile
        sage: with tempfile.NamedTemporaryFile(suffix='.bar') as f:
        ....:     save((1,1), f.name)
        ....:     load(f.name)
        (1, 1)
    """

    if not os.path.splitext(filename)[1] or not hasattr(obj, 'save'):
        filename = _normalize_filename(filename)

    if filename.endswith('.sobj'):
        try:
            obj.save(filename=filename, compress=compress, **kwargs)
        except (AttributeError, RuntimeError, TypeError):
            _base_save(obj, filename, compress=compress)
    else:
        # Saving an object to an image file.
        obj.save(filename, **kwargs)


def _base_dumps(obj, compress=True):
    """
    Base implementation for dumping an object to a pickle in Sage.

    This is the implementation used by :meth:`SageObject.dumps` and by
    :func:`dumps` unless the object being dumped has a custom ``.dumps()``
    method, in which case that is tried first.
    """

    global already_pickled
    gherkin = SagePickler.dumps(obj)
    already_pickled = { }

    if compress:
        return comp.compress(gherkin)

    return gherkin


def dumps(obj, compress=True):
    """
    Dump obj to a string s.  To recover obj, use ``loads(s)``.

    .. SEEALSO:: :func:`loads`

    EXAMPLES::

        sage: a = 2/3
        sage: s = dumps(a)
        sage: a2 = loads(s)
        sage: type(a) is type(a2)
        True
        sage: a2
        2/3
    """
    global already_pickled
    if make_pickle_jar:
        picklejar(obj)
    try:
        ans = obj.dumps(compress)
    except (AttributeError, RuntimeError, TypeError):
        ans = _base_dumps(obj, compress=compress)
    already_pickled = { }
    return ans


# This is used below, and also by explain_pickle.py
unpickle_override = {}


def register_unpickle_override(module, name, callable, call_name=None):
    r"""
    Python pickles include the module and class name of classes.
    This means that rearranging the Sage source can invalidate old
    pickles.  To keep the old pickles working, you can call
    register_unpickle_override with an old module name and class name,
    and the Python callable (function, class with __call__ method, etc.)
    to use for unpickling.  (If this callable is a value in some module,
    you can specify the module name and class name, for the benefit of
    :func:`~sage.misc.explain_pickle.explain_pickle` when called with ``in_current_sage=True``).)

    EXAMPLES:

    Imagine that there used to be an ``old_integer`` module and old
    pickles essentially trying to do the following::

        sage: unpickle_global('sage.rings.old_integer', 'OldInteger')
        Traceback (most recent call last):
        ...
        ImportError: cannot import OldInteger from sage.rings.old_integer, call register_unpickle_override('sage.rings.old_integer', 'OldInteger', ...) to fix this

    After following the advice from the error message, unpickling
    works::

        sage: from sage.misc.persist import register_unpickle_override
        sage: register_unpickle_override('sage.rings.old_integer', 'OldInteger', Integer)
        sage: unpickle_global('sage.rings.old_integer', 'OldInteger')
        <... 'sage.rings.integer.Integer'>

    In many cases, unpickling problems for old pickles can be resolved with a
    simple call to ``register_unpickle_override``, as in the example above and
    in many of the ``sage`` source files.  However, if the underlying data
    structure has changed significantly then unpickling may fail and it
    will be necessary to explicitly implement unpickling methods for the
    associated objects. The python pickle protocol is described in detail on the
    web and, in particular, in the `python pickling documentation`_. For example, the
    following excerpt from this documentation shows that the unpickling of
    classes is controlled by their :meth:`__setstate__` method.

    ::

        object.__setstate__(state)

            Upon unpickling, if the class also defines the method :meth:`__setstate__`, it is
            called with the unpickled state. If there is no :meth:`__setstate__` method,
            the pickled state must be a dictionary and its items are assigned to the new
            instance's dictionary. If a class defines both :meth:`__getstate__` and
            :meth:`__setstate__`, the state object needn't be a dictionary and these methods
            can do what they want.

    .. _python pickling documentation: https://docs.python.org/library/pickle.html#pickle-protocol

    By implementing a :meth:`__setstate__` method for a class it should be
    possible to fix any unpickling problems for the class. As an example of what
    needs to be done, we show how to unpickle a :class:`CombinatorialObject`
    object using a class which also inherits from
    :class:`~sage.structure.element.Element`. This exact problem often arises
    when refactoring old code into the element framework. First we create a
    pickle to play with::

        sage: from sage.structure.element import Element
        sage: class SourPickle(CombinatorialObject): pass
        sage: class SweetPickle(CombinatorialObject, Element): pass
        sage: import __main__
        sage: __main__.SourPickle = SourPickle
        sage: __main__.SweetPickle = SweetPickle  # a hack to allow us to pickle command line classes
        sage: gherkin = dumps(SourPickle([1, 2, 3]))

    Using :func:`register_unpickle_override` we try to sweeten our pickle, but
    we are unable to eat it::

        sage: from sage.misc.persist import register_unpickle_override
        sage: register_unpickle_override('__main__', 'SourPickle', SweetPickle)
        sage: loads(gherkin)
        Traceback (most recent call last):
        ...
        KeyError: 0

    The problem is that the ``SweetPickle`` has inherited a
    :meth:`~sage.structure.element.Element.__setstate__` method from
    :class:`~sage.structure.element.Element` which is not compatible with
    unpickling for :class:`CombinatorialObject`. We can fix this by explicitly
    defining a new :meth:`__setstate__` method::

        sage: class SweeterPickle(CombinatorialObject, Element):
        ....:     def __setstate__(self, state):
        ....:         # a pickle from CombinatorialObject is just its instance
        ....:         # dictionary
        ....:         if isinstance(state, dict):
        ....:             # this is a fudge: we need an appropriate parent here
        ....:             self._set_parent(Tableaux())
        ....:             self.__dict__ = state
        ....:         else:
        ....:             P, D = state
        ....:             if P is not None:
        ....:                 self._set_parent(P)
        ....:             self.__dict__ = D
        sage: __main__.SweeterPickle = SweeterPickle
        sage: register_unpickle_override('__main__', 'SourPickle', SweeterPickle)
        sage: loads(gherkin)                                                            # needs sage.combinat
        [1, 2, 3]
        sage: loads(dumps(SweeterPickle([1, 2, 3])))  # check that pickles work for SweeterPickle
        [1, 2, 3]

    The ``state`` passed to :meth:`__setstate__` will usually be something like
    the instance dictionary of the pickled object, however, with some older
    classes such as :class:`CombinatorialObject` it will be a tuple. In general,
    the ``state`` can be any python object.  ``Sage`` provides a special tool,
    :func:`~sage.misc.explain_pickle.explain_pickle`, which can help in figuring
    out the contents of an old pickle. Here is a second example.

    ::

        sage: class A():
        ....:    def __init__(self,value):
        ....:        self.original_attribute = value
        ....:    def __repr__(self):
        ....:        return 'A(%s)' % self.original_attribute
        sage: class B():
        ....:    def __init__(self,value):
        ....:        self.new_attribute = value
        ....:    def __setstate__(self,state):
        ....:        try:
        ....:            self.new_attribute = state['new_attribute']
        ....:        except KeyError:      # an old pickle
        ....:            self.new_attribute = state['original_attribute']
        ....:    def __repr__(self):
        ....:        return 'B(%s)' % self.new_attribute
        sage: import __main__
        sage: # a hack to allow us to pickle command line classes
        sage: __main__.A = A
        sage: __main__.B = B
        sage: A(10)
        A(10)
        sage: loads(dumps(A(10)))
        A(10)
        sage: sage.misc.explain_pickle.explain_pickle(dumps(A(10)))
        pg_A = unpickle_global('__main__', 'A')
        si = unpickle_newobj(pg_A, ())
        pg_make_integer = unpickle_global('sage.rings.integer', 'make_integer')
        unpickle_build(si, {'original_attribute':pg_make_integer('a')})
        si
        sage: from sage.misc.persist import register_unpickle_override
        sage: register_unpickle_override('__main__', 'A', B)
        sage: loads(dumps(A(10)))
        B(10)
        sage: loads(dumps(B(10)))
        B(10)

    Pickling for python classes and extension classes, such as cython, is
    different -- again this is discussed in the `python pickling
    documentation`_. For the unpickling of extension classes you need to write
    a :meth:`__reduce__` method which typically returns a tuple ``(f,
    args,...)`` such that ``f(*args)`` returns (a copy of) the original object.
    The following code snippet is the
    :meth:`~sage.rings.integer.Integer.__reduce__` method from
    :class:`sage.rings.integer.Integer`.

    .. code-block:: cython

        def __reduce__(self):
            'Including the documentation properly causes a doc-test failure so we include it as a comment:'
            #* '''
            #* This is used when pickling integers.
            #*
            #* EXAMPLES::
            #*
            #*     sage: n = 5
            #*     sage: t = n.__reduce__(); t
            #*     (<built-in function make_integer>, ('5',))
            #*     sage: t[0](*t[1])
            #*     5
            #*     sage: loads(dumps(n)) == n
            #*     True
            #* '''
            # This single line below took me HOURS to figure out.
            # It is the *trick* needed to pickle Cython extension types.
            # The trick is that you must put a pure Python function
            # as the first argument, and that function must return
            # the result of unpickling with the argument in the second
            # tuple as input. All kinds of problems happen
            # if we don't do this.
            return sage.rings.integer.make_integer, (self.str(32),)
    """
    unpickle_override[(module, name)] = (callable, call_name)


def unpickle_global(module, name):
    r"""
    Given a module name and a name within that module (typically a class
    name), retrieve the corresponding object.  This normally just looks
    up the name in the module, but it can be overridden by
    register_unpickle_override.  This is used in the Sage unpickling
    mechanism, so if the Sage source code organization changes,
    register_unpickle_override can allow old pickles to continue to work.

    EXAMPLES::

        sage: from sage.misc.persist import unpickle_override, register_unpickle_override
        sage: unpickle_global('sage.rings.integer', 'Integer')
        <class 'sage.rings.integer.Integer'>

    Now we horribly break the pickling system::

        sage: register_unpickle_override('sage.rings.integer', 'Integer', Rational, call_name=('sage.rings.rational', 'Rational'))
        sage: unpickle_global('sage.rings.integer', 'Integer')
        <class 'sage.rings.rational.Rational'>

    and we reach into the internals and put it back::

        sage: del unpickle_override[('sage.rings.integer', 'Integer')]
        sage: unpickle_global('sage.rings.integer', 'Integer')
        <class 'sage.rings.integer.Integer'>

    A meaningful error message with resolution instructions is displayed for
    old pickles that accidentally got broken because a class or entire module
    was moved or renamed::

        sage: unpickle_global('sage.all', 'some_old_class')
        Traceback (most recent call last):
        ...
        ImportError: cannot import some_old_class from sage.all, call
        register_unpickle_override('sage.all', 'some_old_class', ...)
        to fix this

        sage: unpickle_global('sage.some_old_module', 'some_old_class')
        Traceback (most recent call last):
        ...
        ImportError: cannot import some_old_class from sage.some_old_module, call
        register_unpickle_override('sage.some_old_module', 'some_old_class', ...)
        to fix this

    TESTS:

    Test that :func:`register_unpickle_override` calls in lazily imported modules
    are respected::

        sage: unpickle_global('sage.combinat.root_system.type_A', 'ambient_space')      # needs sage.modules
        <class 'sage.combinat.root_system.type_A.AmbientSpace'>
    """
    unpickler = unpickle_override.get((module, name))
    if unpickler is not None:
        return unpickler[0]

    def error():
        raise ImportError("cannot import {1} from {0}, call "
            "register_unpickle_override({0!r}, {1!r}, ...) to fix this".format(
                module, name))

    mod = sys.modules.get(module)
    if mod is not None:
        try:
            return getattr(mod, name)
        except AttributeError:
            error()
    try:
        __import__(module)
    except ImportError:
        error()

    # Importing the module may have run register_unpickle_override.
    unpickler = unpickle_override.get((module, name))
    if unpickler is not None:
        return unpickler[0]

    mod = sys.modules[module]
    return getattr(mod, name)


class _BasePickler(pickle.Pickler):
    """
    Provides the Python 3 implementation for
    :class:`sage.misc.persist.SagePickler`.

    This is simpler than the Python 2 case since ``pickle.Pickler`` is a
    modern built-in type which can be easily subclassed to provide new
    functionality.

    See the documentation for that class for tests and examples.
    """

    def __init__(self, file_obj, protocol=None, persistent_id=None, *,
                    fix_imports=True):
        super(_BasePickler, self).__init__(file_obj, protocol,
                                            fix_imports=fix_imports)
        self._persistent_id = persistent_id

        def persistent_id(self, obj):
            """
            Implement persistence of external objects with the
            ``persistent_id`` function given at instantiation, if any.
            Otherwise returns ``None`` as in the base class.

            See the documentation for
            :class:`sage.misc.persist.SagePickler` for more details.
            """

            if self._persistent_id is not None:
                return self._persistent_id(obj)

            return super(_BasePickler, self).persistent_id(obj)


# Since Python 3.4, protocol version 4 is the "best" protocol, so we
# use that by default on Sage
_DEFAULT_PROTOCOL_VERSION = 4


class _BaseUnpickler(pickle.Unpickler):
    """
    Provides the Python 3 implementation for
    :class:`sage.misc.persist.SageUnpickler`.

    This is simpler than the Python 2 case since ``pickle.Unpickler`` is
    a modern built-in type which can be easily subclassed to provide new
    functionality.

    See the documentation for that class for tests and examples.
    """

    def __init__(self, file_obj, persistent_load=None, *, **kwargs):
        kwargs.setdefault('encoding', 'latin1')
        super(_BaseUnpickler, self).__init__(file_obj, **kwargs)
        self._persistent_load = persistent_load

    def persistent_load(self, pid):
        """
        Implement persistent loading with the ``persistent_load`` function
        given at instantiation, if any.  Otherwise raises a
        ``pickle.UnpicklingError`` as in the base class.

        See the documentation for :class:`sage.misc.persist.SageUnpickler`
        for more details.
        """

        if self._persistent_load is not None:
            return self._persistent_load(pid)

        return super(_BaseUnpickler, self).persistent_load(pid)

    def find_class(self, module, name):
        """
        The Unpickler uses this class to load module-level objects.
        Contrary to the name, it is used for functions as well as classes.

        (This is equivalent to what was previously called ``find_global``
        which seems like a better name, albeit still somewhat misleading).

        This is just a thin wrapper around
        :func:`sage.misc.persist.unpickle_global`
        """

        # First try using Sage's unpickle_global to go through the unpickle
        # override system
        try:
            return unpickle_global(module, name)
        except ImportError:
            # Failing that, go through the base class's find_class to give
            # it a try (this is necessary in particular to support
            # fix_imports)
            return super(_BaseUnpickler, self).find_class(module, name)


class SagePickler(_BasePickler):
    r"""
    Subclass ``pickle.Pickler`` with Sage-specific default options, and
    built-in support for external object persistence.

    INPUT:

    - ``file_obj`` -- a readable file-like object returning ``bytes`` from
      which the pickle data will be loaded

    - ``persistent_id`` -- callable or None; if given this callable takes a
      single object to be pickled, and returns an "ID" (a key with which to
      restore the object upon unpickling, which may itself be any pickleable
      object).  See the Python documentation on `pickling and unpickling
      external objects`_ for more details.

    - ``py2compat`` -- on Python 3 only, this creates pickles that have a
      better chance of being read on Python 2, by using protocol version 2
      (instead of 4) and fixing up imports of standard library modules and
      types whose names changed between Python 2 and 3.  This is enabled by
      default for the best chances of cross-Python compatibility.

    - Further arguments are passed to :func:`pickle.load`, where in Python-3
      Sage sets the default ``encoding='latin1'``. This is essential to make
      pickles readable in Python-3 that were created in Python-2. See
      :issue:`28444` for details.

    .. _pickling and unpickling external objects: https://docs.python.org/2.7/library/pickle.html#pickling-and-unpickling-external-objects

    EXAMPLES::

        sage: from sage.misc.persist import (
        ....:     unpickle_override, register_unpickle_override, SageUnpickler)
        sage: from sage.rings.integer import make_integer
        sage: from io import BytesIO
        sage: def fake_constructor(x):
        ....:     print("unpickling an Integer")
        ....:     return make_integer(x)
        sage: register_unpickle_override('sage.rings.integer', 'make_integer',
        ....:                            fake_constructor)
        sage: unp = SageUnpickler(BytesIO(dumps(1, compress=False)))
        sage: unp.load()
        unpickling an Integer
        1
        sage: del unpickle_override[('sage.rings.integer', 'make_integer')]

    The ``SagePickler`` can also be passed a ``persistent_id`` function::

        sage: table = {1: 'a', 2: 'b'}
        sage: # in practice this might be a database or something...
        sage: def load_object_from_table(obj_id):
        ....:     tag, obj_id
        ....:     return table[obj_id]

    TESTS:

    The following is an indirect doctest.
    ::

        sage: class Foo():
        ....:     def __init__(self, s):
        ....:         self.bar = s
        ....:     def __reduce__(self):
        ....:         return Foo, (self.bar,)
        ....:
        sage: import __main__
        sage: __main__.Foo = Foo

    The data that is passed to ``loads`` in the following line was created
    by ``dumps(Foo('\x80\x07')`` in Python-2. We demonstrate that it can
    be correctly unpickled in Python-3::

        sage: g = loads(b'x\x9ck`J\x8e\x8f\xcfM\xcc\xcc\x8b\x8f\xe7r\xcb\xcf\xe7*d\x0cej`/dj\r*d\xd6\x03\x00\x89\xc5\x08{')
        sage: type(g), g.bar
        (<class '__main__.Foo'>, '\x80\x07')

    The following line demonstrates what would happen without :issue:`28444`::

        sage: loads(b'x\x9ck`J\x8e\x8f\xcfM\xcc\xcc\x8b\x8f\xe7r\xcb\xcf\xe7*d\x0cej`/dj\r*d\xd6\x03\x00\x89\xc5\x08{', encoding='ASCII')
        Traceback (most recent call last):
        ...
        UnicodeDecodeError: 'ascii' codec can...t decode byte 0x80 in position 0: ordinal not in range(128)
    """

    def __init__(self, file_obj, persistent_id=None, py2compat=True):
        protocol = _DEFAULT_PROTOCOL_VERSION

        if py2compat:
            protocol = 2

        super(SagePickler, self).__init__(file_obj, protocol=protocol,
                                          persistent_id=persistent_id)

    @classmethod
    def dumps(cls, obj, **kwargs):
        """
        Equivalent to :func:`pickle.dumps` but using the
        :class:`sage.misc.persist.SagePickler`.

        INPUT:

        - ``obj`` -- the object to pickle

        - ``kwargs`` -- keyword arguments passed to the
          :class:`sage.misc.persist.SagePickler` constructor

        OUTPUT: ``pickle`` -- the pickled object as ``bytes``

        EXAMPLES::

            sage: import pickle
            sage: from sage.misc.persist import SagePickler
            sage: gherkin = SagePickler.dumps(1)
            sage: pickle.loads(gherkin)
            1
        """
        global already_pickled
        buf = io.BytesIO()
        pickler = cls(buf, **kwargs)
        pickler.dump(obj)
        already_pickled = { }
        return buf.getvalue()


class SageUnpickler(_BaseUnpickler):
    """
    Subclass ``pickle.Unpickler`` to control how certain objects get unpickled
    (registered overrides, specifically).

    This is only needed in Python 3 and up.  On Python 2 the behavior of the
    ``cPickle`` module is customized differently.

    This class simply overrides ``Unpickler.find_class`` to wrap
    ``sage.misc.persist.unpickle_global``.

    INPUT:

    - ``file_obj`` -- a readable file-like object returning ``bytes`` from
      which the pickle data will be loaded

    - ``persistent_load`` -- callable or None; if given this callable
      implements loading of persistent external objects.  The function
      should take a single argument, the persistent object ID. See the
      Python documentation on `pickling and unpickling external objects`_
      for more details.

    - ``kwargs`` -- additional keyword arguments passed to the
      ``pickle.Unpickler`` constructor

    .. _pickling and unpickling external objects: https://docs.python.org/2.7/library/pickle.html#pickling-and-unpickling-external-objects

    EXAMPLES::

        sage: from sage.misc.persist import (
        ....:     unpickle_override, register_unpickle_override, SageUnpickler)
        sage: from sage.rings.integer import make_integer
        sage: from io import BytesIO
        sage: def fake_constructor(x):
        ....:     print("unpickling an Integer")
        ....:     return make_integer(x)
        sage: register_unpickle_override('sage.rings.integer', 'make_integer',
        ....:                            fake_constructor)
        sage: unp = SageUnpickler(BytesIO(dumps(1, compress=False)))
        sage: unp.load()
        unpickling an Integer
        1
        sage: del unpickle_override[('sage.rings.integer', 'make_integer')]

    The ``SageUnpickler`` can also be passed a ``persistent_load`` function::

        sage: table = {1: 'a', 2: 'b'}
        sage: # in practice this might be a database or something...
        sage: def load_object_from_table(obj_id):
        ....:     tag, obj_id
        ....:     return table[obj_id]
    """

    @classmethod
    def loads(cls, data, **kwargs):
        """
        Equivalent to :func:`pickle.loads` but using the
        :class:`sage.misc.persist.SagePickler`.

        INPUT:

        - ``data`` -- the pickle data as ``bytes``

        - ``kwargs`` -- keyword arguments passed to the
          :class:`sage.misc.persist.SageUnpickler` constructor

        OUTPUT: ``obj`` -- the object that was serialized to the given pickle data

        EXAMPLES::

            sage: import pickle
            sage: from sage.misc.persist import SageUnpickler
            sage: gherkin = pickle.dumps(1)
            sage: SageUnpickler.loads(gherkin)
            1
        """

        return cls(io.BytesIO(data), **kwargs).load()


def loads(s, compress=True, **kwargs):
    r"""
    Recover an object x that has been dumped to a string s
    using ``s = dumps(x)``.

    .. SEEALSO:: :func:`dumps`

    EXAMPLES::

        sage: a = matrix(2, [1,2, 3,-4/3])                                              # needs sage.modules
        sage: s = dumps(a)                                                              # needs sage.modules
        sage: loads(s)                                                                  # needs sage.modules
        [   1    2]
        [   3 -4/3]

    If compress is ``True`` (the default), it will try to decompress
    the data with zlib and with bz2 (in turn); if neither succeeds,
    it will assume the data is actually uncompressed.  If ``compress==False``
    is explicitly specified, then no decompression is attempted.
    Further arguments are passed to python's :func:`pickle.load`.

    ::

        sage: v = [1..10]
        sage: loads(dumps(v, compress=False)) == v
        True
        sage: loads(dumps(v, compress=False), compress=True) == v
        True
        sage: loads(dumps(v, compress=True), compress=False)
        Traceback (most recent call last):
        ...
        UnpicklingError: invalid load key, 'x'.

    The next example demonstrates that Sage strives to avoid data loss
    in the transition from Python-2 to Python-3. The problem is that Python-3
    by default would not be able to unpickle a non-ASCII Python-2 string appearing
    in a pickle. See :issue:`28444` for details.
    ::

        sage: class Foo():
        ....:     def __init__(self, s):
        ....:         self.bar = s
        ....:     def __reduce__(self):
        ....:         return Foo, (self.bar,)
        ....:
        sage: import __main__
        sage: __main__.Foo = Foo

    The data that is passed to ``loads`` in the following line was created
    by ``dumps(Foo('\x80\x07')`` in Python-2.
    ::

        sage: g = loads(b'x\x9ck`J\x8e\x8f\xcfM\xcc\xcc\x8b\x8f\xe7r\xcb\xcf\xe7*d\x0cej`/dj\r*d\xd6\x03\x00\x89\xc5\x08{')
        sage: type(g), g.bar
        (<class '__main__.Foo'>, '\x80\x07')

    The following line demonstrates what would happen without :issue:`28444`::

        sage: loads(b'x\x9ck`J\x8e\x8f\xcfM\xcc\xcc\x8b\x8f\xe7r\xcb\xcf\xe7*d\x0cej`/dj\r*d\xd6\x03\x00\x89\xc5\x08{', encoding='ASCII')
        Traceback (most recent call last):
        ...
        UnicodeDecodeError: 'ascii' codec can...t decode byte 0x80 in position 0: ordinal not in range(128)
    """
    if not isinstance(s, bytes):
        raise TypeError("s must be bytes")

    if compress:
        try:
            s = comp.decompress(s)
        except Exception as msg1:
            try:
                s = comp_other.decompress(s)
            except Exception as msg2:
                # Maybe data is uncompressed?
                pass

    unpickler = SageUnpickler(io.BytesIO(s), **kwargs)
    global already_unpickled
    ans = unpickler.load()
    already_unpickled = { }
    return ans


cdef bint make_pickle_jar = 'SAGE_PICKLE_JAR' in os.environ


def picklejar(obj, dir=None):
    """
    Create pickled sobj of ``obj`` in ``dir``, with name the absolute
    value of the hash of the pickle of obj.  This is used in conjunction
    with :func:`unpickle_all`.

    To use this to test the whole Sage library right now, set the
    environment variable ``SAGE_PICKLE_JAR``, which will make it so
    :func:`dumps` will by default call :func:`picklejar` with the
    default dir.  Once you do that and doctest Sage, you'll find that
    the ``DOT_SAGE/pickle_jar`` directory contains a bunch of
    pickled objects along with corresponding txt descriptions of them.
    Use the :func:`unpickle_all` to see if they unpickle later.

    INPUT:

    - ``obj`` -- a pickleable object

    - ``dir`` -- string or ``None``; if ``None`` then ``dir`` defaults to
      ``DOT_SAGE/pickle_jar``

    EXAMPLES::

        sage: dir = tmp_dir()
        sage: sage.misc.persist.picklejar(1, dir)
        sage: sage.misc.persist.picklejar('test', dir)
        sage: len(os.listdir(dir))   # Two entries (sobj and txt) for each object
        4

    TESTS:

    Test an unaccessible directory::

        sage: import os, sys
        sage: s = os.stat(dir)
        sage: os.chmod(dir, 0o000)
        sage: try:
        ....:     uid = os.getuid()
        ....: except AttributeError:
        ....:     uid = -1
        sage: if uid == 0:
        ....:     print("OK (cannot test this as root)")
        ....: else:
        ....:     try:
        ....:         sage.misc.persist.picklejar(1, dir + '/noaccess')
        ....:     except PermissionError:
        ....:         print("OK (correctly raised PermissionError)")
        ....:     else:
        ....:         print("FAIL (did not raise an exception")
        OK...
        sage: os.chmod(dir, s.st_mode)
    """
    if dir is None:
        from sage.env import DOT_SAGE
        dir = os.path.join(DOT_SAGE, 'pickle_jar')
    try:
        os.makedirs(dir)
    except OSError as err:
        # It is not an error if the directory exists
        import errno
        if not err.errno == errno.EEXIST:
            raise

    global already_pickled
    s = comp.compress(SagePickler.dumps(obj))
    already_pickled = { }

    typ = str(type(obj))
    name = ''.join([x if (x.isalnum() or x == '_') else '_' for x in typ])
    base = os.path.join(dir, name)
    if os.path.exists(base):
        i = 0
        while os.path.exists(f'{base}-{i}'):
            i += 1
        base += f'-{i}'

    with open(base + '.sobj', 'wb') as fobj:
        fobj.write(s)

    import sage.version
    stamp = dedent("""\
        type(obj) = {typ}
        version = {ver}
        obj = {obj}
    """.format(typ=typ, ver=sage.version.version, obj=obj))

    with open(base + '.txt', 'w') as fobj:
        fobj.write(stamp)


def unpickle_all(target, debug=False, run_test_suite=False):
    """
    Unpickle all ``.sobj`` files in a directory or tar archive.

    INPUT:

    - ``target`` -- string; the name of a directory or of a (possibly
      compressed) tar archive that contains a single directory of
      ``.sobj`` files.  The tar archive can be in any format that
      python's ``tarfile`` module understands; for example,
      ``.tar.gz`` or ``.tar.bz2``.
    - ``debug`` -- boolean (default: ``False``)
      whether to report a stacktrace in case of failure
    - ``run_test_suite`` -- boolean (default: ``False``)
      whether to run ``TestSuite(x).run()`` on the unpickled objects

    OUTPUT:

    Typically, two lines are printed: the first reporting the number
    of successfully unpickled files, and the second reporting the
    number (zero) of failures. If there are failures, however, then a
    list of failed files will be printed before either of those lines,
    and the failure count will of course be nonzero.

    .. WARNING::

       You must only pass trusted data to this function, including tar
       archives. We use the "data" filter from PEP 706 if possible
       while extracting the archive, but even that is not a perfect
       solution, and it is only available since Python 3.11.4.

    EXAMPLES::

        sage: dir = tmp_dir()
        sage: sage.misc.persist.picklejar('hello', dir)
        sage: sage.misc.persist.unpickle_all(dir)
        Successfully unpickled 1 objects.
        Failed to unpickle 0 objects.
    """
    import os.path
    import tarfile

    ok_count = 0
    fail_count = 0
    failed = []
    tracebacks = []

    if os.path.isfile(target) and tarfile.is_tarfile(target):
        import tempfile
        with tempfile.TemporaryDirectory() as T:
            # Extract the tarball to a temporary directory. The "data"
            # filter only became available in python-3.11.4. See PEP
            # 706 for background.
            with tarfile.open(target) as tf:
                if hasattr(tarfile, "data_filter"):
                    tf.extractall(T, filter='data')
                else:
                    tf.extractall(T)

            # Ensure that the tarball contained exactly one thing, a
            # directory.
            bad_tarball_msg = "tar archive must contain only a single directory"
            contents = os.listdir(T)
            if len(contents) != 1:
                raise ValueError(bad_tarball_msg)

            dir = os.path.join(T, contents[0])
            if not os.path.isdir(dir):
                raise ValueError(bad_tarball_msg)

            # If everything looks OK, start this function over again
            # inside the extracted directory. Note: PEP 343 says the
            # temporary directory will be cleaned up even when the
            # "with" block is exited via a "return" statement. But
            # also note that "return" doesn't happen until the
            # recursive call to unpickle_all() has completed.
            return unpickle_all(dir, debug, run_test_suite)

    if not os.path.isdir(target):
        raise ValueError("target is neither a directory nor a tar archive")

    for A in sorted(os.listdir(target)):
        f = os.path.join(target, A)
        if os.path.isfile(f) and f.endswith('.sobj'):
            try:
                obj = load(f)
                if run_test_suite:
                    TestSuite(obj).run(catch = False)
                ok_count += 1
            except Exception:
                fail_count += 1
                if run_test_suite:
                    print(" * unpickle failure: TestSuite(load('%s')).run()" % f)
                else:
                    print(" * unpickle failure: load('%s')" % f)
                from traceback import print_exc
                print_exc()
                failed.append(A)
                if debug:
                    tracebacks.append(sys.exc_info())

    if failed:
        print("Failed:\n%s" % ('\n'.join(failed)))
    print("Successfully unpickled %s objects." % ok_count)
    print("Failed to unpickle %s objects." % fail_count)
    if debug:
        return tracebacks


def make_None(*args, **kwds):
    """
    Do nothing and return ``None``. Used for overriding pickles when
    that pickle is no longer needed.

    EXAMPLES::

        sage: from sage.misc.persist import make_None
        sage: print(make_None(42, pi, foo='bar'))                                       # needs sage.symbolic
        None
    """
    return None


def load_sage_object(cls, dic):   # not used
    X = cls.__new__(cls)
    try:
        X.__setstate__(dic)
    except AttributeError:
        X.__dict__ = dic
    return X


def load_sage_element(cls, parent, dic_pic):
    X = cls.__new__(cls)
    X._set_parent(parent)
    X.__dict__ = SageUnpickler.loads(dic_pic)
    return X


def db(name):
    r"""
    Load object with given name from the Sage database. Use x.db(name)
    or db_save(x, name) to save objects to the database.

    The database directory is ``$HOME/.sage/db``.
    """
    from sage.misc.misc import SAGE_DB
    return load('%s/%s' % (SAGE_DB, name))


def db_save(x, name=None):
    r"""
    Save x to the Sage database.

    The database directory is ``$HOME/.sage/db``.
    """
    try:
        x.db(name)
    except AttributeError:
        from sage.misc.misc import SAGE_DB
        save(x, '%s/%s' % (SAGE_DB, name))
