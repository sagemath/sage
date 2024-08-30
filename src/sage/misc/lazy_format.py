# sage_setup: distribution = sagemath-objects
"""
Lazy format strings
"""


class LazyFormat(str):
    """
    Lazy format strings.

    .. NOTE::

        We recommend to use :func:`sage.misc.lazy_string.lazy_string` instead,
        which is both faster and more flexible.

    An instance of :class:`LazyFormat` behaves like a usual format
    string, except that the evaluation of the ``__repr__`` method of
    the formatted arguments it postponed until actual printing.

    EXAMPLES:

    Under normal circumstances, :class:`Lazyformat` strings behave as usual::

        sage: from sage.misc.lazy_format import LazyFormat
        sage: LazyFormat("Got `%s`; expected a list")%3
        Got `3`; expected a list
        sage: LazyFormat("Got `%s`; expected %s")%(3, 2/3)
        Got `3`; expected 2/3

    To demonstrate the lazyness, let us build an object with a broken
    ``__repr__`` method::

        sage: class IDontLikeBeingPrinted():
        ....:     def __repr__(self):
        ....:         raise ValueError("do not ever try to print me")

    There is no error when binding a lazy format with the broken object::

        sage: lf = LazyFormat("<%s>")%IDontLikeBeingPrinted()

    The error only occurs upon printing::

        sage: lf
        <repr(<sage.misc.lazy_format.LazyFormat at 0x...>) failed: ValueError: do not ever try to print me>

    .. rubric:: Common use case:

    Most of the time, ``__repr__`` methods are only called during user
    interaction, and therefore need not be fast; and indeed there are
    objects ``x`` in Sage such ``x.__repr__()`` is time consuming.

    There are however some uses cases where many format strings are
    constructed but not actually printed. This includes error handling
    messages in :mod:`unittest` or :class:`TestSuite` executions::

        sage: QQ._tester().assertIn(0, QQ,
        ....:                "%s doesn't contain 0"%QQ)

    In the above ``QQ.__repr__()`` has been called, and the result
    immediately discarded. To demonstrate this we replace ``QQ`` in
    the format string argument with our broken object::

        sage: QQ._tester().assertTrue(True,
        ....:                "%s doesn't contain 0"%IDontLikeBeingPrinted())
        Traceback (most recent call last):
        ...
        ValueError: do not ever try to print me

    This behavior can induce major performance penalties when testing.
    Note that this issue does not impact the usual assert::

        sage: assert True, "%s is wrong"%IDontLikeBeingPrinted()

    We now check that :class:`LazyFormat` indeed solves the assertion problem::

        sage: QQ._tester().assertTrue(True,
        ....:               LazyFormat("%s is wrong")%IDontLikeBeingPrinted())
        sage: QQ._tester().assertTrue(False,
        ....:               LazyFormat("%s is wrong")%IDontLikeBeingPrinted())
        Traceback (most recent call last):
        ...
        AssertionError: ...
    """

    def __mod__(self, args):
        """
        Bind the lazy format string with its parameters.

        EXAMPLES::

            sage: from sage.misc.lazy_format import LazyFormat
            sage: form = LazyFormat("<%s>")
            sage: form
            unbound LazyFormat("<%s>")
            sage: form%"params"
            <params>
        """
        if hasattr(self, "_args"):  # self is already bound...
            self = LazyFormat(""+self)
        self._args = args
        return self

    def __repr__(self):
        """
        TESTS::

            sage: from sage.misc.lazy_format import LazyFormat
            sage: form = LazyFormat("<%s>")
            sage: form
            unbound LazyFormat("<%s>")
            sage: print(form)
            unbound LazyFormat("<%s>")
            sage: form%"toto"
            <toto>
            sage: print(form % "toto")
            <toto>
            sage: print(str(form % "toto"))
            <toto>
            sage: print((form % "toto").__repr__())
            <toto>
        """
        try:
            args = self._args
        except AttributeError:
            return "unbound LazyFormat(\""+self+"\")"
        else:
            return str.__mod__(self, args)

    __str__ = __repr__
