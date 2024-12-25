# sage_setup: distribution = sagemath-objects
r"""
Global options

The :class:`GlobalOptions` class provides a generic mechanism for
setting and accessing **global** options for parents in one or several
related classes, typically for customizing the representation of their
elements. This class will eventually also support setting options on a
parent by parent basis.

These options should be "attached" to one or more classes as an options method.

.. SEEALSO::

    For good examples of :class:`GlobalOptions` in action see
    :obj:`sage.combinat.partition.Partitions.options` and
    :obj:`sage.combinat.tableau.Tableaux.options`.

.. _construction_section:

Construction of options classes
-------------------------------

The general setup for creating a set of global options is::

    sage: from sage.structure.global_options import GlobalOptions
    sage: class MyOptions(GlobalOptions):
    ....:     '''
    ....:     Nice options
    ....:
    ....:     @OPTIONS@
    ....:     '''
    ....:     NAME = 'option name'
    ....:     module = 'sage.some_module.some_file'
    ....:     option_class = 'name_of_class_controlled_by_options'
    ....:     first_option = dict(default='with_bells',
    ....:                       description='Changes the functionality of _repr_',
    ....:                       values=dict(with_bells='causes _repr_ to print with bells',
    ....:                                   with_whistles='causes _repr_ to print with whistles'),
    ....:                       alias=dict(bells='option1', whistles='option2'))
    ....:     # second_option = dict(...)
    ....:     # third_option = dict(...)

Note the syntax using the ``class`` keyword. However, because of some
metaclass magic, the resulting ``MyOptions`` object becomes an instance
of ``GlobalOptions`` instead of a subclass. So, despite the ``class``
syntax, ``MyOptions`` is not a class.

The options constructed by :class:`GlobalOptions` have to be explicitly
associated to the class that they control using the following arguments:

- ``NAME`` -- a descriptive name for the options class; this is
  optional. The default is the name of the constructed class.

- ``module`` -- the sage module containing the options class (optional)

- ``option_class`` -- the name of the options class; this is optional and
  defaults to ``NAME`` if not explicitly set

It is only possible to pickle a :class:`GlobalOptions` class if the
corresponding module is specified *and* if the options are explicitly
attached to the corresponding class as a *options* method.

Each option is specified as a dictionary which describes the possible
values for the option and its documentation. The possible entries in this
dictionary are:

- ``alias`` -- allows for several option values to do the same thing

- ``alt_name`` -- an alternative name for this option

- ``checker`` -- a validation function which returns whether a user
  supplied value is valid or not. This is typically useful for large
  lists of legal values such as :class:`~sage.rings.semirings.non_negative_integer_semiring.NN`.

- ``default`` -- gives the default value for the option

- ``description`` -- a one line description of the option

- ``link_to`` -- links this option to another one in another set of
  global options. This is used for example to allow
  :class:`Partitions` and :class:`Tableaux` to share the same
  ``convention`` option.

- ``setter`` -- a function which is called **after** the value of the
  option is changed

- ``values`` -- dictionary assigning each valid value for the option
  to a short description of what it does

- ``case_sensitive`` -- boolean (default: ``True``); depending on
  whether the values of the option are case sensitive

For each option, either a complete list of possible values, via ``values``, or a
validation function, via ``checker``, must be given. The values can be quite
arbitrary, including user-defined functions which customize the default
behaviour of the classes such as the output of ``_repr_`` or :func:`latex`. See
:ref:`dispatcher` below, and :meth:`~GlobalOptions._dispatcher`, for more
information.

The documentation for the options is automatically constructed from
the docstring of the class by replacing the magic word ``@OPTIONS@``
with a description of each option.

The basic structure for defining a :class:`GlobalOptions` class is best
illustrated by an example::

    sage: from sage.structure.global_options import GlobalOptions
    sage: class Menu():
    ....:     class options(GlobalOptions):
    ....:         '''
    ....:         Fancy documentation
    ....:         -------------------
    ....:
    ....:         @OPTIONS@
    ....:
    ....:         The END!
    ....:         '''
    ....:         NAME = 'menu'
    ....:         entree = dict(default='soup',
    ....:                     description='The first course of a meal',
    ....:                     values=dict(soup='soup of the day', bread='oven baked'),
    ....:                     alias=dict(rye='bread'))
    ....:         appetizer = dict(alt_name='entree')
    ....:         main = dict(default='pizza', description='Main meal',
    ....:                   values=dict(pizza='thick crust', pasta='penne arrabiata'),
    ....:                   case_sensitive=False)
    ....:         dessert = dict(default='espresso', description='Dessert',
    ....:                      values=dict(espresso='life begins again',
    ....:                                  cake='waist begins again',
    ....:                                  cream='fluffy, white stuff'))
    ....:         tip = dict(default=10, description='Reward for good service',
    ....:                    checker = lambda tip: tip in range(0,20))
    sage: Menu.options
    Current options for menu
      - dessert: espresso
      - entree:  soup
      - main:    pizza
      - tip:     10

In the examples above, the options are constructed when the ``options``
object is created. However, it is also possible to construct the options
dynamically using the :meth:`GlobalOptions._add_to_options` methods.

For more details see :class:`GlobalOptions`.

Accessing and setting option values
-----------------------------------

All options and their values, when they are strings, are forced to be lower
case. The values of an options class can be set and accessed by calling the
class or by treating the class as an array.

Continuing the example from :ref:`construction_section`::

    sage: Menu.options
    Current options for menu
      - dessert: espresso
      - entree:  soup
      - main:    pizza
      - tip:     10
    sage: Menu.options.dessert
    espresso
    sage: Menu.options.dessert = 'cake'
    sage: Menu.options.dessert
    cake

Note that, provided there is no ambiguity, options and their values can be
abbreviated::

    sage: Menu.options('d')
    'cake'
    sage: Menu.options('m','t',des='esp', ent='sou')  # get and set several values at once
    ['pizza', 10]
    sage: Menu.options(t=15)
    sage: Menu.options('tip')
    15
    sage: Menu.options.tip
    15
    sage: Menu.options(e='s', m='Pi'); Menu.options()
    Current options for menu
      - dessert: cake
      - entree:  soup
      - main:    pizza
      - tip:     15
    sage: Menu.options(m='P')
    Traceback (most recent call last):
    ...
    ValueError: P is not a valid value for main in the options for menu


Setter functions
----------------

Each option of a :class:`GlobalOptions` can be equipped with an optional setter
function which is called **after** the value of the option is changed. In the
following example, setting the option 'add' changes the state of the class by
setting an attribute in this class using a :func:`classmethod`. Note that the
options object is inserted after the creation of the class in order to access
the :func:`classmethod` as ``A.setter``::

    sage: from sage.structure.global_options import GlobalOptions
    sage: class A(SageObject):
    ....:     state = 0
    ....:     @classmethod
    ....:     def setter(cls, option, val):
    ....:         cls.state += int(val)
    sage: class options(GlobalOptions):
    ....:     NAME = "A"
    ....:     add = dict(default=1,
    ....:                checker=lambda v: int(v)>0,
    ....:                description='An option with a setter',
    ....:                setter=A.setter)
    sage: A.options = options
    sage: A.options
    Current options for A
    - add: 1
    sage: a = A(); a.state
    1
    sage: a.options()
    Current options for A
    - add: 1
    sage: a.options(add=4)
    sage: a.state
    5
    sage: a.options()
    Current options for A
    - add: 4

Documentation for options
-------------------------

The documentation for a :class:`GlobalOptions` is automatically generated from
the supplied options. For example, the generated documentation for the options
``menu`` defined in :ref:`construction_section` is the following:

.. CODE-BLOCK:: text

    Fancy documentation
    -------------------

    OPTIONS:

    - ``appetizer`` -- alternative name for ``entree``
    - ``dessert`` -- (default: ``espresso``)
      Dessert

      - ``cake``     -- waist begins again
      - ``cream``    -- fluffy, white stuff
      - ``espresso`` -- life begins again

    - ``entree`` -- (default: ``soup``)
      The first course of a meal

      - ``bread`` -- oven baked
      - ``rye``   -- alias for ``bread``
      - ``soup``  -- soup of the day

    - ``main`` -- (default: ``pizza``)
      Main meal

      - ``pasta`` -- penne arrabiata
      - ``pizza`` -- thick crust

    - ``tip`` -- (default: ``10``)
      Reward for good service


    The END!

    See :class:`~sage.structure.global_options.GlobalOptions` for more features of these options.

In addition, help on each option, and its list of possible values, can be
obtained by (trying to) set the option equal to '?'::

    sage: Menu.options.dessert?                # not tested
    - ``dessert`` -- (default: ``espresso``)
      Dessert

      - ``cake``     -- waist begins again
      - ``cream``    -- fluffy, white stuff
      - ``espresso`` -- life begins again

.. _dispatcher:

Dispatchers
-----------

The whole idea of a :class:`GlobalOptions` class is that the options change the
default behaviour of the associated classes. This can be done either by simply
checking what the current value of the relevant option is. Another possibility
is to use the options class as a dispatcher to associated methods. To use the
dispatcher feature of a :class:`GlobalOptions` class it is necessary to implement
separate methods for each value of the option where the naming convention for
these methods is that they start with a common prefix and finish with the value
of the option.

If the value of a dispatchable option is set equal to a (user defined) function
then this function is called instead of a class method.

For example, the options ``MyOptions`` can be used to dispatch the ``_repr_``
method of the associated class ``MyClass`` as follows:

.. CODE-BLOCK:: python

    class MyClass(...):
        def _repr_(self):
            return self.options._dispatch(self,'_repr_','first_option')
        def _repr_with_bells(self):
            print('Bell!')
        def _repr_with_whistles(self):
            print('Whistles!')
    class MyOptions(GlobalOptions):
        ...

In this example, ``first_option`` is an option of ``MyOptions`` which takes
values ``bells``, ``whistles``, and so on. Note that it is necessary to make
``self``, which is an instance of ``MyClass``, an argument of the dispatcher
because :meth:`~GlobalOptions._dispatch()` is a method of :class:`GlobalOptions`
and not a method of ``MyClass``. Apart from ``MyOptions``, as it is a method of
this class, the arguments are the attached class (here ``MyClass``), the prefix
of the method of ``MyClass`` being dispatched, the option of ``MyOptions``
which controls the dispatching. All other arguments are passed through to the
corresponding methods of ``MyClass``. In general, a dispatcher is invoked as:

.. CODE-BLOCK:: python

    self.options._dispatch(self, dispatch_to, option, *args, **kargs)

Usually this will result in the method
``dispatch_to + '_' + MyOptions(options)`` of ``self`` being called with
arguments ``*args`` and ``**kargs`` (if ``dispatch_to[-1] == '_'`` then the
method ``dispatch_to + MyOptions(options)`` is called).

If ``MyOptions(options)`` is itself a function then the dispatcher will call
this function instead. In this way, it is possible to allow the user to
customise the default behaviour of this method. See
:meth:`~GlobalOptions._dispatch` for an example of how this can be achieved.

The dispatching capabilities of :class:`GlobalOptions` allows options to be
applied automatically without needing to parse different values of the option
(the cost is that there must be a method for each value). The dispatching
capabilities can also be used to make one option control several methods:

.. CODE-BLOCK:: python

    def __le__(self, other):
        return self.options._dispatch(self, '_le_','cmp', other)
    def __ge__(self, other):
        return self.options._dispatch(self, '_ge_','cmp', other)
    def _le_option_a(self, other):
        return ...
    def _ge_option_a(self, other):
        return ...
    def _le_option_b(self, other):
        return ...
    def _ge_option_b(self, other):
        return ...

See :meth:`~GlobalOptions._dispatch` for more details.

Doc testing
-----------

All of the options and their effects should be doc-tested. However, in order
not to break other tests, all options should be returned to their default state
at the end of each test. To make this easier, every :class:`GlobalOptions` class has
a :meth:`~GlobalOptions._reset()` method for doing exactly this.


Pickling
--------

Options classes can only be pickled if they are the options for some standard
sage class. In this case the class is specified using the arguments to
:class:`GlobalOptions`. For example
:meth:`~sage.combinat.partition.Partitions.options` is defined as:

.. CODE-BLOCK:: python

     class Partitions(UniqueRepresentation, Parent):
         ...
         class options(GlobalOptions):
             NAME = 'Partitions'
             module = 'sage.combinat.partition'
             ...

Here is an example to test the pickling of a :class:`GlobalOptions` instance::

    sage: TestSuite(Partitions.options).run()                                           # needs sage.combinat

TESTS:

Check that the old call syntax still works::

    sage: class Menu():
    ....:     options = GlobalOptions('menu',
    ....:         doc='Fancy documentation\n'+'-'*19, end_doc='The END!',
    ....:         entree=dict(default='soup',
    ....:                     description='The first course of a meal',
    ....:                     values=dict(soup='soup of the day', bread='oven baked'),
    ....:                     alias=dict(rye='bread')),
    ....:         appetizer=dict(alt_name='entree'),
    ....:         main=dict(default='pizza', description='Main meal',
    ....:                   values=dict(pizza='thick crust', pasta='penne arrabiata'),
    ....:                   case_sensitive=False),
    ....:         dessert=dict(default='espresso', description='Dessert',
    ....:                      values=dict(espresso='life begins again',
    ....:                                  cake='waist begins again',
    ....:                                  cream='fluffy, white stuff')),
    ....:         tip=dict(default=10, description='Reward for good service',
    ....:         checker=lambda tip: tip in range(0,20))
    ....:     )
    sage: Menu.options
    Current options for menu
      - dessert: espresso
      - entree:  soup
      - main:    pizza
      - tip:     10

We can have a ``name`` option::

    sage: class MyOptions(GlobalOptions):
    ....:     name = dict(default='alice', values={'alice': "A", 'bob': "B"})
    sage: MyOptions
    Current options for MyOptions
      - name: alice

Check that the ``name`` and ``NAME`` keywords are both supported with
this syntax::

    sage: GlobalOptions(name='menu')
    Current options for menu
    sage: GlobalOptions(NAME='menu')
    Current options for menu
    sage: GlobalOptions()
    Traceback (most recent call last):
    ...
    TypeError: GlobalOptions() is missing keyword argument 'name'

Test the documentation examples above::

    sage: print(Menu.options.__doc__)
    Fancy documentation
    -------------------
    <BLANKLINE>
    OPTIONS:
    <BLANKLINE>
    - ``appetizer`` -- alternative name for ``entree``
    - ``dessert`` -- (default: ``espresso``)
      Dessert
    <BLANKLINE>
      - ``cake``     -- waist begins again
      - ``cream``    -- fluffy, white stuff
      - ``espresso`` -- life begins again
    <BLANKLINE>
    - ``entree`` -- (default: ``soup``)
      The first course of a meal
    <BLANKLINE>
      - ``bread`` -- oven baked
      - ``rye``   -- alias for ``bread``
      - ``soup``  -- soup of the day
    <BLANKLINE>
    - ``main`` -- (default: ``pizza``)
      Main meal
    <BLANKLINE>
      - ``pasta`` -- penne arrabiata
      - ``pizza`` -- thick crust
    <BLANKLINE>
    - ``tip`` -- (default: ``10``)
      Reward for good service
    <BLANKLINE>
    <BLANKLINE>
    <BLANKLINE>
    <BLANKLINE>
    The END!
    See :class:`~sage.structure.global_options.GlobalOptions` for more features of these options.

::

    sage: print(Menu.options.dessert.__doc__)
    - ``dessert`` -- (default: ``espresso``)
      Dessert
    <BLANKLINE>
      - ``cake``     -- waist begins again
      - ``cream``    -- fluffy, white stuff
      - ``espresso`` -- life begins again
    <BLANKLINE>

AUTHORS:

- Andrew Mathas (2013): initial version
- Andrew Mathas (2016): overhaul making the options attributes, enabling
                        pickling and attaching the options to a class.
- Jeroen Demeyer (2017): use subclassing to create instances
"""

# ****************************************************************************
#       Copyright (C) 2013,2016 Andrew Mathas <andrew dot mathas at sydney dot edu dot au>
#       Copyright (C) 2017 Jeroen Demeyer <J.Demeyer@UGent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from importlib import import_module
from pickle import PicklingError
from textwrap import dedent

from sage.misc.instancedoc import instancedoc


class Option:
    r"""
    An option.

    Each option for an options class is an instance of this class which
    implements the magic that allows the options to the attributes of the
    options class that can be looked up, set and called.

    By way of example, this class implements the following functionality.

    EXAMPLES::

        sage: # needs sage.combinat
        sage: Partitions.options.display
        list
        sage: Partitions.options.display = 'compact'
        sage: Partitions.options.display('list')
        sage: Partitions.options._reset()

    TESTS::

        sage: TestSuite(Partitions.options.display).run()                               # needs sage.combinat
    """
    __name__ = 'Option class'

    def __init__(self, options, name):
        r"""
        Initialise an option by settings its ``name``, "parent" option class
        ``options`` and doc-string.

        EXAMPLES::

            sage: type(Partitions.options.display)    # indirect doctest                # needs sage.combinat
            <class 'sage.structure.global_options.Option'>
        """
        self._name = name
        self._options = options
        self.__doc__ = options._doc[name]
        super().__init__()

    def __repr__(self):
        r"""
        Return a string representation for this collection of options.

        EXAMPLES::

            sage: Partitions.options.display  # indirect doctest                        # needs sage.combinat
            list
        """
        # NOTE: we intentionally use str() instead of repr()
        return str(self._options[self._name])

    def __add__(self, other):
        r"""
        Return the object obtained by adding ``self`` and ``other``, where
        ``self`` behaves like its value.

        EXAMPLES::

            sage: Tableaux.options.convention + ' is good'                              # needs sage.combinat
            'English is good'
        """
        return self._options[self._name] + other

    def __radd__(self, other):
        r"""
        Return the object obtained by adding ``other`` and ``self``, where
        ``self`` behaves like its value.

        EXAMPLES::

            sage: 'Good ' + Tableaux.options.convention                                 # needs sage.combinat
            'Good English'
        """
        return other + self._options[self._name]

    def __mul__(self, other):
        r"""
        Return the object obtained by adding ``self`` and ``other``, where
        ``self`` behaves like its value.

        EXAMPLES::

            sage: Tableaux.options.convention + ' is good'                              # needs sage.combinat
            'English is good'
        """
        return self._options[self._name] * other

    def __rmul__(self, other):
        r"""
        Return the object obtained by r-adding ``other`` and ``self``, where
        ``self`` behaves like its value.

        EXAMPLES::

            sage: 'Good ' + Tableaux.options.convention                                 # needs sage.combinat
            'Good English'
        """
        return other * self._options[self._name]

    def __bool__(self) -> bool:
        r"""
        Return the value of this option interpreted as a boolean.

        EXAMPLES::

            sage: # needs sage.combinat sage.graphs sage.modules
            sage: RiggedConfigurations.options.half_width_boxes_type_B
            True
            sage: 'yes' if RiggedConfigurations.options.half_width_boxes_type_B else 'no'
            'yes'
            sage: RiggedConfigurations.options.half_width_boxes_type_B = False
            sage: 'yes' if RiggedConfigurations.options.half_width_boxes_type_B else 'no'
            'no'
            sage: RiggedConfigurations.options._reset()
        """
        return bool(self._options[self._name])

    def __call__(self, *args, **kwds):
        r"""
        Get or set value of the option ``self``.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: Partitions.options.display()       # indirect doctest
            'list'
            sage: Partitions.options.display('exp')  # indirect doctest
            sage: Partitions.options.display()       # indirect doctest
            'exp_low'
            sage: Partitions.options._reset()

        TESTS:

        Check that values can be set to ``None`` (:issue:`30763`)::

            sage: from sage.structure.global_options import GlobalOptions
            sage: class config(GlobalOptions):
            ....:     size = dict(default=42,
            ....:                 description='integer or None',
            ....:                 checker=lambda val: val is None or val >= 0)
            sage: config.size()
            42
            sage: config.size(None)
            sage: config.size() is None
            True
            sage: config._reset()

        Check the deprecation::

            sage: config.size(value=None)
            doctest:...: DeprecationWarning: keyword argument "value" should be replaced by positional argument
            See https://github.com/sagemath/sage/issues/30763 for details.
            sage: config.size() is None
            True
            sage: config.size(1, 2)
            Traceback (most recent call last):
            ...
            TypeError: option takes at most one argument "value"
            sage: config.size(unknown=3)
            Traceback (most recent call last):
            ...
            TypeError: option takes at most one argument "value"
            sage: config.size(4, value=5)
            Traceback (most recent call last):
            ...
            TypeError: option takes at most one argument "value"
        """
        if not args and not kwds:
            return self._options[self._name]
        if 'value' in kwds:
            from sage.misc.superseded import deprecation
            deprecation(30763, 'keyword argument "value" should be replaced '
                               'by positional argument')
            args += (kwds.pop('value'),)
        if len(args) > 1 or kwds:
            raise TypeError('option takes at most one argument "value"')
        self._options[self._name] = args[0]

    def __eq__(self, other):
        r"""
        Equality testing for an option in based on the value of the attribute.

        EXAMPLES::

            sage: Tableaux.options.convention                                           # needs sage.combinat
            English
            sage: Tableaux.options.convention == "English"                              # needs sage.combinat
            True
            sage: Tableaux.options.convention == "French"                               # needs sage.combinat
            False
        """
        return self._options[self._name] == other

    def __ne__(self, other):
        r"""
        Inequality testing for an option in based on the value of
        the attribute.

        EXAMPLES::

            sage: Tableaux.options.convention                                           # needs sage.combinat
            English
            sage: Tableaux.options.convention != "English"                              # needs sage.combinat
            False
            sage: Tableaux.options.convention != "French"                               # needs sage.combinat
            True
        """
        return self._options[self._name] != other

    def __hash__(self):
        r"""
        Return the hash of ``self``, which is the hash of the corresponding
        value.

        EXAMPLES::

            sage: hash(Tableaux.options.convention) == hash(Tableaux.options('convention'))         # needs sage.combinat
            True
        """
        return hash(self._options[self._name])

    def __str__(self):
        r"""
        Return the string representation of ``self``, which is the string of
        the corresponding value.

        EXAMPLES::

            sage: str(Tableaux.options.convention)                                      # needs sage.combinat
            'English'
        """
        return str(self._options[self._name])


class GlobalOptionsMetaMeta(type):
    def __call__(self, name, bases, dict):
        """
        Called when subclassing an instance of ``self``.

        Python translates ``class C(B): ...`` to
        ``meta = type(B); C = meta("C", (B,), ...)``.
        If we want to intercept this call ``meta(...)``, we need to
        customize ``__call__`` in the metaclass of ``meta``, which is
        this metametaclass.

        EXAMPLES::

            sage: from sage.structure.global_options import GlobalOptions
            sage: type(GlobalOptions)
            <class 'sage.structure.global_options.GlobalOptionsMeta'>
            sage: type(type(GlobalOptions))
            <class 'sage.structure.global_options.GlobalOptionsMetaMeta'>
            sage: class G(GlobalOptions): pass
            sage: type(G)
            <class 'sage.structure.global_options.GlobalOptions'>

        Since ``G`` is constructed using ``class`` syntax, the object
        gets a ``__module__`` attribute, different from the
        ``__module__`` attribute of its type (:class:`GlobalOptions`)::

            sage: G.__module__
            '__main__'
            sage: type(G).__module__
            'sage.structure.global_options'

        Multiple base classes are not allowed::

            sage: class G(GlobalOptions, object): pass
            Traceback (most recent call last):
            ...
            TypeError: GlobalOptions must be the only base class
        """
        # Allow only a single base class (which is GlobalOptions in
        # practice).
        # If we ever find a reasonable meaning for multiple base
        # classes, note that Python 2 and Python 3 have different
        # semantics for determining the metaclass when multiple base
        # classes are involved.
        #
        # Note: On Python 3 bases is empty if the class was declared
        # without any explicted bases:
        if not bases:
            bases = (object,)

        if len(bases) != 1:
            raise TypeError("GlobalOptions must be the only base class")

        base = bases[0]
        if not isinstance(base, self):
            # For the creation of the initial GlobalOptions class, fall
            # back to the usual class creation
            return self.__new__(self, name, bases, dict)

        # Create an instance of the base class (as opposed to an
        # instance of self, which would be typical for a metaclass)
        instance = base.__new__(base)

        # Split dict in options for instance.__init__ and attributes to
        # insert in the class __dict__
        kwds = {"NAME": name}
        for key, value in dict.items():
            if key.startswith("__"):
                instance.__dict__[key] = value
            else:
                kwds[key] = value

        instance.__init__(**kwds)
        return instance


class GlobalOptionsMeta(type, metaclass=GlobalOptionsMetaMeta):
    """
    Metaclass for :class:`GlobalOptions`.

    This class is itself an instance of :class:`GlobalOptionsMetaMeta`,
    which implements the subclass magic.
    """


@instancedoc
class GlobalOptions(metaclass=GlobalOptionsMeta):
    r"""
    The :class:`GlobalOptions` class is a generic class for setting and
    accessing global options for Sage objects.

    While it is possible to create instances of :class:`GlobalOptions`
    the usual way, the recommended syntax is to subclass from
    ``GlobalOptions``. Thanks to some metaclass magic, this actually
    creates an instance of ``GlobalOptions`` instead of a subclass.

    INPUT (as "attributes" of the class):

    - ``NAME`` -- specifies a name for the options class (optional;
      default: class name)

    - ``module`` -- gives the module that contains the associated options class

    - ``option_class`` -- gives the name of the associated module class
      (default: ``NAME``)

    - option = ``dict(...)`` -- dictionary specifying an option

    The options are specified by keyword arguments with their values
    being a dictionary which describes the option. The
    allowed/expected keys in the dictionary are:

    - ``alias`` -- defines alias/synonym for option values
    - ``alt_name`` -- alternative name for an option
    - ``checker`` -- a function for checking whether a particular value for
      the option is valid
    - ``default`` -- the default value of the option
    - ``description`` -- documentation string
    - ``link_to`` -- links to an option for this set of options to an
      option in another :class:`GlobalOptions`
    - ``setter`` -- a function (class method) which is called whenever this
      option changes
    - ``values`` -- dictionary of the legal values for this option (this
      automatically defines the corresponding ``checker``); this dictionary
      gives the possible options, as keys, together with a brief description
      of them
    - ``case_sensitive`` -- boolean (default: ``True``); depending on whether
      the values of the option are case sensitive

    Options and their values can be abbreviated provided that this
    abbreviation is a prefix of a unique option.

    EXAMPLES::

        sage: from sage.structure.global_options import GlobalOptions
        sage: class Menu():
        ....:     class options(GlobalOptions):
        ....:         '''
        ....:         Fancy documentation
        ....:         -------------------
        ....:
        ....:         @OPTIONS@
        ....:
        ....:         End of Fancy documentation
        ....:         '''
        ....:         NAME = 'menu'
        ....:         entree = dict(default='soup',
        ....:                     description='The first course of a meal',
        ....:                     values=dict(soup='soup of the day', bread='oven baked'),
        ....:                     alias=dict(rye='bread'))
        ....:         appetizer = dict(alt_name='entree')
        ....:         main = dict(default='pizza', description='Main meal',
        ....:                   values=dict(pizza='thick crust', pasta='penne arrabiata'),
        ....:                   case_sensitive=False)
        ....:         dessert = dict(default='espresso', description='Dessert',
        ....:                      values=dict(espresso='life begins again',
        ....:                                  cake='waist begins again',
        ....:                                  cream='fluffy white stuff'))
        ....:         tip = dict(default=10, description='Reward for good service',
        ....:                  checker=lambda tip: tip in range(0,20))
        sage: Menu.options
        Current options for menu
          - dessert: espresso
          - entree:  soup
          - main:    pizza
          - tip:     10
        sage: Menu.options(entree='s')         # unambiguous abbreviations are allowed
        sage: Menu.options(t=15)
        sage: (Menu.options['tip'], Menu.options('t'))
        (15, 15)
        sage: Menu.options()
        Current options for menu
          - dessert: espresso
          - entree:  soup
          - main:    pizza
          - tip:     15
        sage: Menu.options._reset(); Menu.options()
        Current options for menu
          - dessert: espresso
          - entree:  soup
          - main:    pizza
          - tip:     10
        sage: Menu.options.tip=40
        Traceback (most recent call last):
        ...
        ValueError: 40 is not a valid value for tip in the options for menu
        sage: Menu.options(m='p')           # ambiguous abbreviations are not allowed
        Traceback (most recent call last):
        ...
        ValueError: p is not a valid value for main in the options for menu

    The documentation for the options class is automatically generated from the
    information which specifies the options:

    .. CODE-BLOCK:: text

        Fancy documentation
        -------------------

        OPTIONS:

        - dessert:  (default: espresso)
          Dessert

          - ``cake``     -- waist begins again
          - ``cream``    -- fluffy white stuff
          - ``espresso`` -- life begins again

        - entree:  (default: soup)
          The first course of a meal

          - ``bread`` -- oven baked
          - ``rye``   -- alias for bread
          - ``soup``  -- soup of the day

        - main:  (default: pizza)
          Main meal

          - ``pasta`` -- penne arrabiata
          - ``pizza`` -- thick crust

        - tip:  (default: 10)
          Reward for good service

        End of Fancy documentation

        See :class:`~sage.structure.global_options.GlobalOptions` for more features of these options.

    The possible values for an individual option can be obtained by
    (trying to) set it equal to '?'::

        sage: Menu.options(des='?')
        - ``dessert`` -- (default: ``espresso``)
          Dessert
        <BLANKLINE>
          - ``cake``     -- waist begins again
          - ``cream``    -- fluffy white stuff
          - ``espresso`` -- life begins again
        <BLANKLINE>
        Current value: espresso
    """
    __name__ = 'options'

    def __init__(self, NAME=None, module='', option_class='', doc='', end_doc='', **options):
        r"""
        Initialize ``self``.

        TESTS::

            sage: from sage.structure.global_options import GlobalOptions
            sage: class menu(GlobalOptions):
            ....:     entree = dict(default='soup',
            ....:                 description='The first course of a meal',
            ....:                 values=dict(soup='soup of the day', bread='oven baked'),
            ....:                 alias=dict(rye='bread'))
            ....:     appetizer = dict(alt_name='entree')
            ....:     main = dict(default='pizza', description='Main meal',
            ....:               values=dict(pizza='thick crust', pasta='penne arrabiata'),
            ....:               case_sensitive=False)
            ....:     dessert = dict(default='espresso', description='Dessert',
            ....:                  values=dict(espresso='life begins again',
            ....:                              cake='waist begins again',
            ....:                              cream='fluffy white stuff'))
            ....:     tip = dict(default=10, description='Reward for good service',
            ....:              checker=lambda tip: tip in range(0,20))
            sage: menu._name  # Default name is class name
            'menu'
            sage: class specials(GlobalOptions):
            ....:     entree = dict(link_to=(menu, 'entree'))
            ....:     main_specials = dict(default='salmon', description='main course specials',
            ....:                   values=dict(salmon='a fish', crab='Sebastian'))
            sage: specials['entree'] = 'rye'
            sage: menu['entree']
            'bread'

            sage: class alias_test(GlobalOptions):
            ....:       "Test aliases with case sensitivity"
            ....:       test_opt = dict(default='Upper',
            ....:           description = 'Starts with an uppercase',
            ....:           values = dict(Upper="Starts with uppercase",
            ....:                         lower="only lowercase"),
            ....:           case_sensitive = False,
            ....:           alias = dict(UpperAlias='Upper', lower_alias='lower'))
            sage: alias_test['test_opt'] = 'Lower_Alias'
            sage: alias_test['test_opt']
            'lower'
            sage: alias_test['test_opt'] = 'upperalias'
            sage: alias_test['test_opt']
            'Upper'
        """
        if NAME is None:
            # Require a "name" keyword in **options
            try:
                NAME = options.pop("name")
            except KeyError:
                raise TypeError("GlobalOptions() is missing keyword argument 'name'")
        self._name = NAME

        # initialise the various dictionaries used by GlobalOptions
        self._alias = {}           # a dictionary of alias for the values of some options
        self._alt_names = {}       # a dictionary of alternative names for some options
        self._case_sensitive = {}  # a dictionary of booleans indicating to check case sensitivity
        self._checker = {}         # a dictionary of validity checkers for each option
        self.__default_value = {}  # a dictionary of the default options
        self._display_values = {}  # a dictionary of the output of the values
        self._doc = {}             # a dictionary of doc strings, forced by the linked options
        self._legal_values = {}    # a dictionary of lists of the legal values for each option
        self._linked_value = {}    # a dictionary of linked to other global options as (link, linked_option)
        self._setter = {}          # a dictionary of the list of setters
        self._value = {}           # a dictionary of the current options
        for option in options:
            self._add_option(option, options[option])

        self._options_module = module
        self._option_class = option_class or self._name

        # If the instance dict has a __doc__ attribute, use that as
        # docstring.
        try:
            self._docstring = dedent(self.__dict__.pop("__doc__"))
        except KeyError:
            self._docstring = "@OPTIONS@"

        # Finally, we build the doc string for the options
        # First we strip common white space off the front of doc and end_doc
        if doc:
            self._docstring = dedent(doc) + "\n\n" + self._docstring

        if end_doc:
            self._docstring = self._docstring + "\n\n" + dedent(end_doc)

        # Add docstring footer
        self._docstring += "\nSee :class:`~sage.structure.global_options.GlobalOptions` for more features of these options."

    def __repr__(self):
        r"""
        Return a string representation for this collection of options.

        EXAMPLES::

            sage: from sage.structure.global_options import GlobalOptions
            sage: class FoodOptions(GlobalOptions):
            ....:     NAME = 'daily meal'
            ....:     food = dict(default='apple', values=dict(apple='a fruit', pair='of what?'))
            ....:     drink = dict(default='water', values=dict(water='a drink', coffee='a lifestyle'))
            sage: FoodOptions
            Current options for daily meal
              - drink: water
              - food:  apple
        """
        options = list(self._value) + list(self._linked_value)
        for x in self._alt_names:
            options.remove(x)
        if not options:
            return 'Current options for {}'.format(self._name)

        options.sort()
        width = 1 + max(len(option) for option in options)
        txt = '\n'.join('  - {:{}} {}'.format(option + ':', width, self[option])
                        for option in options)
        return 'Current options for {}\n{}'.format(self._name, txt)

    def __call__(self, *get_value, **set_value):
        r"""
        Get or set value of the option ``option``.

        EXAMPLES::

            sage: from sage.structure.global_options import GlobalOptions
            sage: class FoodOptions(GlobalOptions):
            ....:     NAME = 'daily meal'
            ....:     food = dict(default='apple', values=dict(apple='a fruit', pair='of what?'))
            ....:     drink = dict(default='water', values=dict(water='a drink', coffee='a lifestyle'))
            ....:     beverage = dict(alt_name='drink')
            sage: FoodOptions()
            Current options for daily meal
              - drink: water
              - food:  apple
            sage: FoodOptions('food')
            'apple'
            sage: FoodOptions(food='pair'); FoodOptions()
            Current options for daily meal
              - drink: water
              - food:  pair
            sage: FoodOptions('beverage')
            'water'
            sage: FoodOptions(food='apple', drink='coffee'); FoodOptions()
            Current options for daily meal
              - drink: coffee
              - food:  apple
        """
        if not get_value and not set_value:
            return self

        if get_value:
            # use __getitem__ to return these options
            if len(get_value) == 1:
                return self[get_value[0]]
            else:
                return [self[option] for option in get_value]

        # use __setitem__ to set these options
        if set_value:
            for option in set_value:
                self[option] = set_value[option]

    def __getitem__(self, option):
        r"""
        Return the current value of the option ``option``.

        EXAMPLES::

            sage: from sage.structure.global_options import GlobalOptions
            sage: class FoodOptions(GlobalOptions):
            ....:     NAME = 'daily meal'
            ....:     food = dict(default='apple', values=dict(apple='a fruit', pair='of what?'))
            ....:     drink = dict(default='water', values=dict(water='a drink', coffee='a lifestyle'))
            sage: FoodOptions['drink']
            'water'
            sage: FoodOptions['d']
            'water'
        """
        option = self._match_option(option)
        if option in self._linked_value:
            link, linked_opt = self._linked_value[option]
            return link[linked_opt]
        elif option in self._value:
            if option in self._display_values:
                return self._display_values[option][self._value[option]]
            return self._value[option]

    def __setitem__(self, option, value):
        r"""
        The ``__setitem__`` method is used to change the current values of the
        options. It also checks that the supplied options are valid and changes
        any alias to its generic value.

        EXAMPLES::

            sage: from sage.structure.global_options import GlobalOptions
            sage: class FoodOptions(GlobalOptions):
            ....:     NAME = 'daily meal'
            ....:     food = dict(default='apple', values=dict(apple='a fruit', pair='of what?'))
            ....:     drink = dict(default='water', values=dict(water='a drink', coffee='a lifestyle'))
            sage: FoodOptions['drink']='coffee'; FoodOptions()
            Current options for daily meal
              - drink: coffee
              - food:  apple
            sage: FoodOptions(drink='w'); FoodOptions()
            Current options for daily meal
              - drink: water
              - food:  apple
            sage: FoodOptions(drink='?')
            - ``drink`` -- (default: ``water``)
            <BLANKLINE>
              - ``coffee`` -- a lifestyle
              - ``water``  -- a drink
            <BLANKLINE>
            Current value: water
        """
        option = self._match_option(option)
        if not callable(value):
            value = self._match_value(option, value)

        if value == '?':  # return help
            print('%s\nCurrent value: %s' % (self._doc[option], self[option]))
            return      # we do not want to call the setter below

        elif option in self._linked_value:
            link, linked_opt = self._linked_value[option]
            link[linked_opt] = value

        else:
            self._value[option] = value

        if option in self._setter:
            # if a setter function exists then call it with the associated
            # class, option and value
            self._setter[option](option, value)

    def _instancedoc_(self):
        r"""
        Return the docstring for the options class ``self``.

        EXAMPLES::

            sage: print(Partitions.options.__doc__)                                     # needs sage.combinat
            <BLANKLINE>
            Set and display the global options for elements of the partition,
            skew partition, and partition tuple classes.  If no parameters are
            set, then the function returns a copy of the options dictionary.
            <BLANKLINE>
            The ``options`` to partitions can be accessed as the method
            :obj:`Partitions.options` of :class:`Partitions` and
            related parent classes.
            <BLANKLINE>
            <BLANKLINE>
            OPTIONS:
            <BLANKLINE>
            - ``convention`` -- (default: ``English``)
              Sets the convention used for displaying tableaux and partitions
            <BLANKLINE>
              - ``English`` -- use the English convention
              - ``French``  -- use the French convention
            ...

        TESTS::

            sage: from sage.structure.global_options import GlobalOptions
            sage: print(GlobalOptions.__doc__)
            <BLANKLINE>
                The :class:`GlobalOptions` class is a generic class for setting...
        """
        options = 'OPTIONS:\n\n{}\n\n'.format('\n'.join(self._doc[opt] for opt in sorted(self._doc)))
        return self._docstring.replace("@OPTIONS@", options)

    def __setattr__(self, name, value=None):
        r"""
        Set the attribute ``name`` of the option class ``self`` equal to
        ``value``, if the attribute ``name`` exists.

        As the attributes of an option class are the actual options we need
        to be able to "trap" invalid options in a sensible way. We do this
        by sending any "non-standard" to :meth:`__setitem__` for processing.

        EXAMPLES::

            sage: Partitions.options.display = 'exp'                                    # needs sage.combinat
            sage: Partitions.options.dispplay = 'list'                                  # needs sage.combinat
            Traceback (most recent call last):
            ...
            ValueError: dispplay is not an option for Partitions
            sage: Partitions.options._reset()                                           # needs sage.combinat
        """
        # Underscore, and "special", attributes are set using type.__setattr__
        if name[0] == '_' or name in ['reset', 'dispatch', 'default_value']:
            return super().__setattr__(name, value)

        # General case: redirect to __setitem__
        self[name] = value

    def __setstate__(self, state):
        r"""
        This is a custom :meth:`__setstate__` method for unpickling instances of
        the :class:`GlobalOptions` class.

        The :meth:`__getstate__` method returns a dictionary with an
        ``options_class`` key which identifies the "parent" class for the options.
        This is then used to unpickle the options class.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: Partitions.options()
            Current options for Partitions
              - convention:        English
              - diagram_str:       *
              - display:           list
              - latex:             young_diagram
              - latex_diagram_str: \ast
            sage: Partitions.options.convention = "French"
            sage: pickle = dumps(Partitions.options)
            sage: Partitions.options._reset()        # reset options
            sage: loads(pickle)                      # indirect doctest
            Current options for Partitions
              - convention:        French
              - diagram_str:       *
              - display:           list
              - latex:             young_diagram
              - latex_diagram_str: \ast
            sage: Partitions.options._reset()
        """
        # open the options for the corresponding "parent" and copy all of
        # the data from its options class into unpickle
        options_class = getattr(import_module(state['options_module']), state['option_class'])
        unpickle = options_class.options
        state.pop('option_class')
        state.pop('options_module')
        for setting in unpickle.__dict__:
            self.__dict__[setting] = unpickle.__dict__[setting]

        # reset the options in ``self`` to their defaults
        self._reset()
        # apply the options stored in state
        for opt in state:
            self[opt] = state[opt]
        options_class.options = self

    def __getstate__(self):
        r"""
        Return a dictionary that can be used to pickle an instance of a
        :class:`GlobalOptions` class.

        This is called by :func:`pickle_GlobalOptions`.

        Typically instances of :class:`GlobalOptions` are
        complicated to create, so they do no pickle. If the options are
        associated to  "parent" class then it is possible to create the
        default version of the class and then add the non-default
        settings on top of this. This method returns a dictionary of the
        non-default options that can then be used to unpickle an
        instance of the option using :func:`unpickle_GlobalOptions`.

        EXAMPLES::

            sage: Partitions.options._reset()                                           # needs sage.combinat
            sage: Partitions.options.__getstate__()                                     # needs sage.combinat
             {'convention': 'English',
             'option_class': 'Partitions',
             'options_module': 'sage.combinat.partition'}
        """

        # options classes can be pickled only if they are the options for an
        # associated "parent" class that lives in self._module

        if not self._options_module:
            pickleable = False
        else:
            opt_mod = import_module(self._options_module)
            pickleable = hasattr(opt_mod, self._name) and hasattr(getattr(opt_mod, self._name), 'options')

        if not pickleable:
            raise PicklingError('%s cannot be pickled because it is not associated with a class' % self)

        pickle = {'option_class': self._option_class, 'options_module': self._options_module}
        for opt in self._value:
            if opt not in self._alt_names and self[opt] != self.__default_value[opt]:
                pickle[opt] = self[opt]
        for opt in self._linked_value:
            link, linked_opt = self._linked_value[opt]
            if opt not in self._alt_names and link[opt] != link.__default_value[opt]:
                pickle[opt] = self[opt]
        return pickle

    def __eq__(self, other):
        r"""
        Two options classes are equal if they return the same :meth:`__getstate__`.

        EXAMPLES::

            sage: Partitions.options == PartitionsGreatestLE.options  # indirect doctest            # needs sage.combinat
            True
            sage: Partitions.options == Tableaux.options                                # needs sage.combinat
            False
        """
        return isinstance(other, GlobalOptions) and self.__getstate__() == other.__getstate__()

    def _add_option(self, option, specifications):
        r"""
        Add an option.

        INPUT:

        - ``option`` -- string
        - ``specifications`` -- dictionary

        .. SEEALSO::

            :class:`GlobalOptions` for a description of the required
            ``specifications``.
        """
        if not isinstance(specifications, dict):
            raise TypeError("expected dict as specification of %r, got %r" % (option, specifications))

        doc = {}  # will be used to build the doc string
        self._case_sensitive[option] = True    # ``True`` by default
        self._legal_values[option] = []
        for spec in sorted(specifications):   # NB: options processed alphabetically!
            if spec == 'alias':
                self._alias[option] = specifications[spec]
                self._legal_values[option] += list(specifications[spec])
                for opt in specifications[spec]:
                    doc[opt] = 'alias for ``%s``' % specifications[spec][opt]
            elif spec == 'alt_name':
                self._alt_names[option] = specifications[spec]
                self._linked_value[option] = (self, specifications[spec])
                doc = '- ``%s`` -- alternative name for ``%s``' % (option, specifications[spec].lower())
            elif spec == 'case_sensitive':
                if not specifications[spec]:
                    for opt in self._legal_values:
                        self._display_values[option] = {val.lower(): val for val in self._legal_values[option]}
                        self._legal_values[option] = [val.lower() for val in self._legal_values[option]]
                    if option in self._alias:
                        self._alias[option] = {k.lower(): v.lower()
                                               for k, v in self._alias[option].items()}
                self._case_sensitive[option] = bool(specifications[spec])
            elif spec == 'checker':
                if not callable(specifications[spec]):
                    raise ValueError('the checker for %s must be callable' % option)
                self._checker[option] = specifications[spec]
            elif spec == 'default':
                self.__default_value[option] = specifications[spec]
            elif spec == 'link_to':
                if (isinstance(specifications[spec], tuple) and
                        len(specifications[spec]) == 2 and
                        isinstance(specifications[spec][0], GlobalOptions)):
                    link, linked_opt = specifications['link_to']  # for sanity
                    if linked_opt in link._value:
                        self._linked_value[option] = specifications['link_to']
                        link, linked_opt = specifications['link_to']  # for sanity
                        doc = link._doc[linked_opt]
                    elif linked_opt in link._linked_value:
                        self._linked_value[option] = link._linked_value[linked_opt]
                        doc = link._doc[linked_opt]
                    else:
                        raise ValueError("could not find link to {1} in {0}".format(*specifications[spec]))
                else:
                    raise ValueError("linked options must be specified as a string: 'linked_option' or a tuple: (link,linked_option)")
            elif spec == 'setter':
                if callable(specifications[spec]):
                    self._setter[option] = specifications[spec]
                else:
                    raise ValueError('the setter for %s must be a function' % option)
            elif spec == 'values':
                for val in specifications[spec]:
                    doc[val] = specifications[spec][val]
                doc.update(specifications[spec])
                if self._case_sensitive[option]:
                    self._legal_values[option] += list(specifications[spec])
                    self._display_values[option] = {val: val for val in specifications[spec]}
                else:
                    self._legal_values[option] += [val.lower() for val in specifications[spec]]
                    self._display_values[option] = {val.lower(): val for val in specifications[spec]}
            elif spec != 'description':
                raise ValueError('Initialization error in Global options for %s: %s not recognized!' % (self._name, spec))

        # now build the doc string for this option
        if doc == {} and 'description' not in specifications:
            raise ValueError('no documentation specified for %s in the options for %s' % (option, self._name))

        # first a necessary hack to initialise the option in self._doc because __setitem__ calls _match_option
        self._doc[option] = ''
        if option in self._linked_value:
            self._doc[option] = doc
        else:
            width = max(len(v) for v in doc) + 4 if doc != {} else 4
            if len(doc) > 0:
                self._doc[option] = '- ``{}`` -- (default: ``{}``)\n{}\n{}\n'.format(
                    option, self._default_value(option),
                    '  %s\n' % specifications['description'] if 'description' in specifications else '',
                    '\n'.join('  - {:{}} -- {}'.format('``' + val + '``', width, doc[val])
                              for val in sorted(doc)))
            else:
                self._doc[option] = '- ``{}`` -- (default: ``{}``)\n{}'.format(
                    option, self._default_value(option),
                    '  %s\n' % specifications['description'] if 'description' in specifications else '')

        # sanity check for non-linked options
        if option not in self._linked_value:
            if 'default' not in specifications:
                raise ValueError('a default value for %s must be given' % option)

            if not (option in self._checker or option in self._legal_values):
                raise ValueError('a value checker or a list of valid values for %s must be given' % option)

            # finally, set, check and process the default value using  __setitem__
            self[option] = self.__default_value[option]
            self.__default_value[option] = self._value[option]  # in case the default is an alias

        # Build getters and setters for this option. As we have
        # overridden __setattr__, we call object.__setattr__ directly
        super().__setattr__(option, Option(self, option))

    def _match_option(self, option):
        r"""
        Check whether ``option`` is the name of an option of ``self``, or a
        prefix thereof.

        INPUT:

        - ``option`` -- string

        EXAMPLES::

            sage: from sage.structure.global_options import GlobalOptions
            sage: class FoodOptions(GlobalOptions):
            ....:     NAME = 'daily meal'
            ....:     food = dict(default='apple', values=dict(apple='a fruit', pair='of what?'))
            ....:     drink = dict(default='water', values=dict(water='a drink', coffee='a lifestyle'))
            sage: FoodOptions('food') # indirect doctest
            'apple'
            sage: FoodOptions('f')
            'apple'
        """
        if option in self._doc:
            return option

        # a lower case version of the option
        loption = option.lower()

        # as it is not an option try and match it with a prefix to an option,
        # without checking case using the fact that the keys of self._doc is a
        # list of the options, both normal and linked
        matches = [opt for opt in self._doc if opt.lower().startswith(loption)]
        if matches and all(m.startswith(matches[0]) for m in matches):
            return matches[0]
        elif len(matches) > 1:
            # as there is more than one match check case as well
            matches = [mat for mat in matches if mat.startswith(option)]
            if matches and all(m.startswith(matches[0]) for m in matches):
                return matches[0]
            else:
                raise ValueError('%s is an ambiguous option for %s' % (option, self._name))

        # if we are still here this is not a good option!
        raise ValueError('%s is not an option for %s' % (option, self._name))

    def _match_value(self, option, value):
        r"""
        Check whether ``value`` is a valid value for ``option``.

        INPUT:

        - ``option`` -- string: the name of an option, or prefix thereof
        - ``value`` -- a value or ``'?'``

        EXAMPLES::

            sage: from sage.structure.global_options import GlobalOptions
            sage: class FoodOptions(GlobalOptions):
            ....:     NAME = 'daily meal'
            ....:     food = dict(default='apple', values=dict(apple='a fruit', pair='of what?'))
            ....:     drink = dict(default='water', values=dict(water='a drink', wine='a lifestyle'))
            sage: FoodOptions(f='a') # indirect doctest
            sage: FoodOptions('f')
            'apple'
            sage: FoodOptions(d='wi'); FoodOptions('f')
            'apple'
            sage: FoodOptions(d='w')
            Traceback (most recent call last):
            ...
            ValueError: w is not a valid value for drink in the options for daily meal
        """
        if value == "?":
            return value   # help on this value

        if option in self._linked_value:
            link, linked_opt = self._linked_value[option]
            return link._match_value(linked_opt, value)

        orig_value = value

        # convert to proper case if it is a string
        if isinstance(value, str) and not self._case_sensitive[option]:
            value = value.lower()

        if option in self._legal_values:
            if value in self._legal_values[option]:
                if option in self._alias and value in self._alias[option]:
                    return self._alias[option][value]
                return value

            # as it is not a value try and match it with a prefix of a value
            matches = [val for val in self._legal_values[option] if val.startswith(value)]
            if matches and all(m.startswith(matches[0]) for m in matches):
                val = matches[0]
                if option in self._alias and val in self._alias[option]:
                    return self._alias[option][val]
                return val

        if option in self._checker and self._checker[option](value):
            return value

        # if we are still here this is not a good value!

        # replace any value alias with its "real" value
        if option in self._alias and value in self._alias[option]:
            orig_value = self._alias[option][value]
        raise ValueError('%s is not a valid value for %s in the options for %s' % (orig_value, option, self._name))

    def _default_value(self, option):
        r"""
        Return the default value of the option.

        EXAMPLES::

            sage: from sage.structure.global_options import GlobalOptions
            sage: class FoodOptions(GlobalOptions):
            ....:     NAME = 'daily meal'
            ....:     food = dict(default='apple', values=dict(apple='a fruit', pair='of what?'))
            ....:     drink = dict(default='water', values=dict(water='a drink', coffee='a lifestyle'))
            sage: FoodOptions._default_value('food')
            'apple'
        """
        option = self._match_option(option)
        if option in self.__default_value:
            return self.__default_value[option]
        else:
            link, linked_opt = self._linked_value[option]
            return link._default_value(linked_opt)

    def _dispatch(self, obj, dispatch_to, option, *args, **kargs):
        r"""
        The *dispatchable* options are options which dispatch related methods of
        the corresponding class - or user defined methods which are passed to
        :class:`GlobalOptions`. The format for specifying a dispatchable option
        is to include ``dispatch_to = <option name>`` in the specifications for
        the options and then to add the options to the (element) class. Each
        option is then assumed to be a method of the element class with a name
        of the form ``<option name> + '_' + <current vale for this option'``.
        These options are called by the element class via::

            return self.options._dispatch(self, dispatch_to, option, *args, **kargs)

        Note that the argument ``self`` is necessary here because the
        dispatcher is a method of the options class and not of
        ``self``. The value of the option can also be set to a
        user-defined function, with arguments ``self`` and ``option``;
        in this case the user's function is called instead.

        EXAMPLES:

        Here is a contrived example::

            sage: from sage.structure.global_options import GlobalOptions
            sage: class DelimitedListOptions(GlobalOptions):
            ....:           delim=dict(default='b', values={'b':'brackets', 'p':'parentheses'})
            sage: class DelimitedList(SageObject):
            ....:    options = DelimitedListOptions
            ....:    def __init__(self, L):
            ....:        self._list = L
            ....:    def _repr_b(self): return '[%s]' % ','.join('%s'%i for i in self._list)
            ....:    def _repr_p(self): return '(%s)' % ','.join('%s'%i for i in self._list)
            ....:    def _repr_(self): return self.options._dispatch(self, '_repr_','delim')
            sage: dlist = DelimitedList([1,2,3]); dlist
            [1,2,3]
            sage: dlist.options.delim='p'; dlist
            (1,2,3)
            sage: dlist.options(delim=lambda self: '<%s>' % ','.join('%s'%i for i in self._list)); dlist
            <1,2,3>
        """
        if callable(self._value[option]):
            try:
                return self._value[option](obj, *args, **kargs)
            except TypeError:
                raise ValueError('the user defined dispatcher function failed!')
        else:
            if dispatch_to[-1] == '_':
                dispatch_to = dispatch_to[:-1]
            dispatch = getattr(obj, dispatch_to + '_' + self._value[option])
            return dispatch(*args, **kargs)

        raise ValueError('%s is not a dispatchable option!' % option)

    def _reset(self, option=None):
        r"""
        Reset options to their default value.

        INPUT:

        - ``option`` -- (default: ``None``) the name of an option as a string
          or ``None``. If ``option`` is specified only this option is reset to
          its default value; otherwise all options are reset.

        EXAMPLES::

            sage: from sage.structure.global_options import GlobalOptions
            sage: class Meal():
            ....:     class options(GlobalOptions):
            ....:         NAME = 'daily meal'
            ....:         food = dict(default='bread', values=dict(bread='rye bread', salmon='a fish'))
            ....:         drink = dict(default='water',values=dict(water='essential for life', wine='essential'))
            sage: Meal.options.food='salmon'; Meal.options
            Current options for daily meal
              - drink: water
              - food:  salmon
            sage: Meal.options._reset('drink'); Meal.options()
            Current options for daily meal
              - drink: water
              - food:  salmon
            sage: Meal.options._reset(); Meal.options()
            Current options for daily meal
              - drink: water
              - food:  bread
        """
        if option is None:
            for option in self.__default_value:
                self._value[option] = self.__default_value[option]
                if not self._case_sensitive[option] and isinstance(self._value[option], str):
                    self._value[option] = self._value[option].lower()
            for option in self._linked_value:
                link, linked_opt = self._linked_value[option]
                link._reset(linked_opt)
        else:
            option = self._match_option(option)
            if option in self.__default_value:
                self._value[option] = self.__default_value[option]
                if not self._case_sensitive[option] and isinstance(self._value[option], str):
                    self._value[option] = self._value[option].lower()
            elif option in self._linked_value:
                link, linked_opt = self._linked_value[option]
                link._reset(linked_opt)
