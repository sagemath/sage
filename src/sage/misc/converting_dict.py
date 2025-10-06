r"""
Converting Dictionary

At the moment, the only class contained in this module is a key
converting dictionary, which applies some function (e.g. type
conversion function) to all arguments used as keys.

.. It is conceivable that a other dicts might be added later on.

AUTHORS:

- Martin von Gagern (2015-01-31): initial version

EXAMPLES:

A ``KeyConvertingDict`` will apply a conversion function to all method
arguments which are keys::

    sage: from sage.misc.converting_dict import KeyConvertingDict
    sage: d = KeyConvertingDict(int)
    sage: d["3"] = 42
    sage: list(d.items())
    [(3, 42)]

This is used e.g. in the result of a variety, to allow access to the
result no matter how a generator is identified::

    sage: # needs sage.libs.singular sage.rings.number_field
    sage: K.<x,y> = QQ[]
    sage: I = ideal([x^2 + 2*y - 5, x + y + 3])
    sage: V = sorted(I.variety(AA), key=str)
    sage: v = V[0]
    sage: v['x'], v['y']
    (-2.464101615137755?, -0.535898384862246?)
    sage: list(v)[0].parent()
    Multivariate Polynomial Ring in x, y over Algebraic Real Field
"""
# ****************************************************************************
#       Copyright (C) 2015 Martin von Gagern <Martin.vGagern@gmx.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from collections.abc import Mapping


class KeyConvertingDict(dict):
    r"""
    A dictionary which automatically applies a conversions to its keys.

    The most common application is the case where the conversion
    function is the object representing some category, so that key
    conversion means a type conversion to adapt keys to that
    category. This allows different representations for keys which in
    turn makes accessing the correct element easier.

    INPUT:

    - ``key_conversion_function`` -- a function which will be
      applied to all method arguments which represent keys
    - ``data`` -- (optional) dictionary or sequence of key-value pairs
      to initialize this mapping

    EXAMPLES::

        sage: from sage.misc.converting_dict import KeyConvertingDict
        sage: d = KeyConvertingDict(int)
        sage: d["3"] = 42
        sage: list(d.items())
        [(3, 42)]
        sage: d[5.0] = 64
        sage: d["05"]
        64
    """

    def __init__(self, key_conversion_function, data=None):
        r"""
        Construct a dictionary with a given conversion function.

        EXAMPLES::

            sage: from sage.misc.converting_dict import KeyConvertingDict
            sage: d = KeyConvertingDict(int)
            sage: d["3"] = 42
            sage: list(d.items())
            [(3, 42)]
            sage: list(KeyConvertingDict(int, {"5": 7}).items())
            [(5, 7)]
            sage: list(KeyConvertingDict(int, [("9", 99)]).items())
            [(9, 99)]
        """
        super().__init__()
        self.key_conversion_function = key_conversion_function
        if data:
            self.update(data)

    def __getitem__(self, key):
        r"""
        Retrieve an element from the dictionary.

        INPUT:

        - ``key`` -- a value identifying the element, will be converted

        EXAMPLES::

            sage: from sage.misc.converting_dict import KeyConvertingDict
            sage: d = KeyConvertingDict(int)
            sage: d[3] = 42
            sage: d["3"]
            42
        """
        key = self.key_conversion_function(key)
        return super().__getitem__(key)

    def __setitem__(self, key, value):
        r"""
        Assign an element in the dictionary.

        INPUT:

        - ``key`` -- a value identifying the element, will be converted
        - ``value`` -- the associated value, will be left unmodified

        EXAMPLES::

            sage: from sage.misc.converting_dict import KeyConvertingDict
            sage: d = KeyConvertingDict(int)
            sage: d["3"] = 42
            sage: list(d.items())
            [(3, 42)]
        """
        key = self.key_conversion_function(key)
        return super().__setitem__(key, value)

    def __delitem__(self, key):
        r"""
        Remove a mapping from the dictionary.

        INPUT:

        - ``key`` -- a value identifying the element, will be converted

        EXAMPLES::

            sage: from sage.misc.converting_dict import KeyConvertingDict
            sage: d = KeyConvertingDict(int)
            sage: d[3] = 42
            sage: del d["3"]
            sage: len(d)
            0
        """
        key = self.key_conversion_function(key)
        return super().__delitem__(key)

    def __contains__(self, key):
        r"""
        Test whether a given key is contained in the mapping.

        INPUT:

        - ``key`` -- a value identifying the element, will be converted

        EXAMPLES::

            sage: from sage.misc.converting_dict import KeyConvertingDict
            sage: d = KeyConvertingDict(int)
            sage: d[3] = 42
            sage: "3" in d
            True
            sage: 4 in d
            False
        """
        key = self.key_conversion_function(key)
        return super().__contains__(key)

    def pop(self, key, *args):
        r"""
        Remove and retrieve a given element from the dictionary.

        INPUT:

        - ``key`` -- a value identifying the element, will be converted
        - ``default`` -- the value to return if the element is not mapped, optional

        EXAMPLES::

            sage: from sage.misc.converting_dict import KeyConvertingDict
            sage: d = KeyConvertingDict(int)
            sage: d[3] = 42
            sage: d.pop("3")
            42
            sage: d.pop("3", 33)
            33
            sage: d.pop("3")
            Traceback (most recent call last):
            ...
            KeyError: ...
        """
        key = self.key_conversion_function(key)
        return super().pop(key, *args)

    def setdefault(self, key, default=None):
        r"""
        Create a given mapping unless there already exists a mapping
        for that key.

        INPUT:

        - ``key`` -- a value identifying the element, will be converted
        - ``default`` -- the value to associate with the key

        EXAMPLES::

            sage: from sage.misc.converting_dict import KeyConvertingDict
            sage: d = KeyConvertingDict(int)
            sage: d.setdefault("3")
            sage: list(d.items())
            [(3, None)]
        """
        key = self.key_conversion_function(key)
        return super().setdefault(key, default)

    def update(self, *args, **kwds):
        r"""
        Update the dictionary with key-value pairs from another dictionary,
        sequence of key-value pairs, or keyword arguments.

        INPUT:

        - ``key`` -- a value identifying the element, will be converted
        - ``args`` -- a single dict or sequence of pairs
        - ``kwds`` -- named elements require that the conversion
          function accept strings

        EXAMPLES::

            sage: from sage.misc.converting_dict import KeyConvertingDict
            sage: d = KeyConvertingDict(int)
            sage: d.update([("3",1),(4,2)])
            sage: d[3]
            1
            sage: d.update({"5": 7, "9": 12})
            sage: d[9]
            12
            sage: d = KeyConvertingDict(QQ['x'])
            sage: d.update(x=42)
            sage: d
            {x: 42}
        """
        f = self.key_conversion_function
        u = super().update
        if args:
            if len(args) != 1:
                raise TypeError("update expected at most 1 argument")
            arg = args[0]
            if isinstance(arg, Mapping):
                seq = ((f(k), arg[k]) for k in arg)
            else:
                seq = ((f(k), v) for k, v in arg)
            u(seq)
        if kwds:
            seq = ((f(k), v) for k, v in kwds.items())
            u(seq)

    def _repr_pretty_(self, p, cycle):
        """
        For pretty printing in the Sage command prompt.

        Since ``KeyConvertingDict`` inherits from ``dict``, we just use IPython's
        built-in ``dict`` pretty printer.
        When :issue:`36801` is fixed, this function will be redundant.

        EXAMPLES::

            sage: from sage.misc.converting_dict import KeyConvertingDict
            sage: d = KeyConvertingDict(int)
            sage: d["3"] = 4
            sage: d["1"] = 2
            sage: repr(d)    # dictionaries are insertion ordered since Python 3.6
            '{3: 4, 1: 2}'
            sage: d          # indirect doctest
            {1: 2, 3: 4}

        The last example output will be ``{3: 4, 1: 2}`` outside of doctesting,
        see :func:`sage.doctest.forker.init_sage`.
        """
        p.pretty(dict(self))
