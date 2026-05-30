
import copyreg


class bar:
    pass


def metaclass(name, bases):
    """
    Create a new class in this metaclass.

    INPUT:

    - ``name`` -- string
    - ``bases`` -- tuple of classes

    EXAMPLES::

        sage: from sage.misc.test_class_pickling import metaclass, bar
        sage: c = metaclass("foo2", (bar, object))
        constructing class
        sage: c
        <class 'sage.misc.test_class_pickling.foo2'>
        sage: type(c)
        <class 'sage.misc.test_class_pickling.Metaclass'>
        sage: c.__bases__
        (<...sage.misc.test_class_pickling.bar...>, <... 'object'>)
    """
    print("constructing class")
    result = Metaclass(name, bases, {})
    result.reduce_args = (name, bases)
    return result


class Metaclass(type):
    """
    This metaclass illustrates the customization of how a class is pickled.
    It requires a slightly patched version of cPickle.

    See:

    - https://docs.python.org/3/library/copyreg.html#module-copyreg
    - http://groups.google.com/group/comp.lang.python/browse_thread/thread/66c282afc04aa39c/
    - http://groups.google.com/group/sage-devel/browse_thread/thread/583048dc7d373d6a/

    EXAMPLES::

        sage: from sage.misc.test_class_pickling import metaclass, bar
        sage: c = metaclass("foo", (bar, object))
        constructing class
        sage: import pickle
        sage: s = pickle.dumps(c)
        reducing a class
        sage: c2 = pickle.loads(s)
        constructing class
        sage: c == c2
        calling __eq__ defined in Metaclass
        True
    """
    def __eq__(self, other):
        print("calling __eq__ defined in Metaclass")
        return (type(self) is type(other)) and (self.reduce_args == other.reduce_args)

    def __reduce__(self):
        """
        Implement the pickle protocol for classes in this metaclass
        (not for the instances of this class!!!)

        EXAMPLES::

            sage: from sage.misc.test_class_pickling import metaclass, bar
            sage: c = metaclass("foo3", (bar, object))
            constructing class
            sage: c.__class__.__reduce__(c)
            reducing a class
            (<function metaclass at ...>,
             ('foo3', (<...sage.misc.test_class_pickling.bar...>, <...'object'>)))
        """
        print("reducing a class")
        return (metaclass, self.reduce_args)


copyreg.pickle(Metaclass, Metaclass.__reduce__)
