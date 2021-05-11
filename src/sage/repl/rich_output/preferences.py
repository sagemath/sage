# -*- encoding: utf-8 -*-
r"""
Display Preferences

This class is used to express display preferences that are not simply
a choice of a particular output format. For example, whether to prefer
vector over raster graphics. By convention, the value ``None`` is
always a valid value for a preference and means no particular
preference.

EXAMPLES::

    sage: from sage.repl.rich_output.preferences import DisplayPreferences
    sage: prefs = DisplayPreferences()
    sage: prefs.available_options()
    (graphics, supplemental_plot, text)
    sage: prefs.text is None
    True
    sage: prefs.text = 'ascii_art'
    sage: prefs.text
    'ascii_art'
    sage: prefs
    Display preferences:
    * graphics is not specified
    * supplemental_plot is not specified
    * text = ascii_art

Properties can be unset by deleting them or by assigning ``None``::

    sage: prefs.text = 'ascii_art'
    sage: del prefs.text
    sage: prefs.text is None
    True

    sage: prefs.text = 'ascii_art'
    sage: prefs.text = None
    sage: prefs.text is None
    True

Properties have documentation attached::

    sage: import pydoc
    sage: doc = pydoc.render_doc(prefs)
    sage: assert ' graphics' in doc
    sage: assert '     Preferred graphics format' in doc
    sage: assert ' text' in doc
    sage: assert '     Which textual representation is preferred' in doc

Values can also be specified as keyword arguments to the constructor::

    sage: DisplayPreferences(text='latex')
    Display preferences:
    * graphics is not specified
    * supplemental_plot is not specified
    * text = latex

.. TODO::

    A value-checking preference system should be used elsewhere in
    Sage, too. The class here is just a simple implementation, a
    proper implementation would use a metaclass to construct the
    preference items.
"""

# ****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from textwrap import dedent

from sage.structure.sage_object import SageObject


class Property(property):

    def __init__(self, name, allowed_values, doc=None):
        r"""
        Preference item

        INPUT:

        - ``name`` -- string. The name of the property.

        - ``allowed_values`` -- list/tuple/iterable of allowed values.

        - ``doc`` -- string (optional). The docstring of the property.

        EXAMPLES::

            sage: from sage.repl.rich_output.preferences import Property
            sage: prop = Property('foo', [0, 1, 2], 'The Foo Property')
            sage: prop.__doc__
            'The Foo Property\n\nAllowed values:\n\n* ``None`` (default): no preference\n\n* 0\n\n* 1\n\n* 2'
            sage: prop.allowed_values
            (0, 1, 2)
        """
        self.name = name
        self.underscore_name = '_{0}'.format(name)
        self.allowed_values = tuple(allowed_values)
        self.__doc__ = doc = self._make_doc(doc)
        super(Property, self).__init__(
            fget=self.getter, fset=self.setter, fdel=self.deleter, doc=doc)

    def _make_doc(self, doc):
        """
        Generate the documentation.

        INPUT:

        - ``doc`` -- the title line of the documentation.

        OUTPUT:

        String. The docstring with auto-generated documentation about
        the allowed values added.

        EXAMPLES::

            sage: from sage.repl.rich_output.preferences import Property
            sage: prop = Property('foo', [0, 1, 2], 'The Foo Property')
            sage: print(prop._make_doc('this is the title'))
            this is the title
            <BLANKLINE>
            Allowed values:
            <BLANKLINE>
            * ``None`` (default): no preference
            <BLANKLINE>
            * 0
            <BLANKLINE>
            * 1
            <BLANKLINE>
            * 2
        """
        doc = dedent(doc)
        doc += '\n\n'
        doc += 'Allowed values:\n\n'
        values_doc = []
        values_doc.append('* ``None`` (default): no preference')
        for value in self.allowed_values:
            values_doc.append('* {0}'.format(repr(value)))
        return doc + '\n\n'.join(values_doc)

    def __repr__(self):
        """
        Return a string representation

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.repl.rich_output.preferences import Property
            sage: prop = Property('foo', [0, 1, 2], 'The Foo Property')
            sage: prop.__repr__()
            'foo'
        """
        return self.name

    def getter(self, prefs):
        """
        Get the current value of the property

        INPUT:

        - ``prefs`` -- the :class:`PreferencesABC` instance that the
          property is bound to.

        OUTPUT:

        One of the allowed values or ``None`` if not set.

        EXAMPLES::

            sage: from sage.repl.rich_output.preferences import Property, PreferencesABC
            sage: prop = Property('foo', [0, 1, 2], 'The Foo Property')
            sage: prefs = PreferencesABC()
            sage: prop.getter(prefs) is None
            True
            sage: prop.setter(prefs, 1)
            sage: prop.getter(prefs)
            1
        """
        try:
            return getattr(prefs, self.underscore_name)
        except AttributeError:
            return None

    def setter(self, prefs, value):
        """
        Get the current value of the property

        INPUT:

        - ``prefs`` -- the :class:`PreferencesABC` instance that the
          property is bound to.

        - ``value`` -- anything. The new value of the
          property. Setting a property to ``None`` is equivalent to
          deleting the value.

        OUTPUT:

        This method does not return anything. A ``ValueError`` is
        raised if the given ``value`` is not one of the allowed
        values.

        EXAMPLES::

            sage: from sage.repl.rich_output.preferences import Property, PreferencesABC
            sage: prop = Property('foo', [0, 1, 2], 'The Foo Property')
            sage: prefs = PreferencesABC()
            sage: prop.getter(prefs) is None
            True
            sage: prop.setter(prefs, 1)
            sage: prop.getter(prefs)
            1

            sage: prop.setter(prefs, None)
            sage: prop.getter(prefs) is None
            True
        """
        if value is None:
            return self.deleter(prefs)
        allowed = self.allowed_values
        if value not in allowed:
            raise ValueError('value must be unset (None) or one of {0}, got {1}'
                             .format(allowed, value))
        setattr(prefs, self.underscore_name, value)

    def deleter(self, prefs):
        """
        Delete the current value of the property

        INPUT:

        - ``prefs`` -- the :class:`PreferencesABC` instance that the
          property is bound to.

        EXAMPLES::

            sage: from sage.repl.rich_output.preferences import Property, PreferencesABC
            sage: prop = Property('foo', [0, 1, 2], 'The Foo Property')
            sage: prefs = PreferencesABC()
            sage: prop.getter(prefs) is None
            True
            sage: prop.setter(prefs, 1)
            sage: prop.deleter(prefs)
            sage: prop.getter(prefs) is None
            True
        """
        underscore_name = self.underscore_name
        try:
            delattr(prefs, underscore_name)
        except AttributeError:
            pass


class PreferencesABC(SageObject):

    def __init__(self, *args, **kwds):
        """
        Preferences for displaying graphics

        These can be preferences expressed by the user or by the
        display backend. They are specified as keyword arguments.

        INPUT:

        - ``*args*`` -- positional arguments are preferences
          instances. The property values will be inherited from left
          to right, that is, later parents override values from
          earlier parents.

        - ``**kwds`` -- keyword arguments. Will be used to initialize
          properties, and override inherited values if necessary.

        EXAMPLES::

            sage: from sage.repl.rich_output.preferences import DisplayPreferences
            sage: p1 = DisplayPreferences(graphics='vector')
            sage: p2 = DisplayPreferences(graphics='raster')
            sage: DisplayPreferences(p1, p2)
            Display preferences:
            * graphics = raster
            * supplemental_plot is not specified
            * text is not specified

        If specified in the opposite order, the setting from ``p1`` is
        inherited::

            sage: DisplayPreferences(p2, p1)
            Display preferences:
            * graphics = vector
            * supplemental_plot is not specified
            * text is not specified

        Further keywords override::

            sage: DisplayPreferences(p2, p1, graphics='disable')
            Display preferences:
            * graphics = disable
            * supplemental_plot is not specified
            * text is not specified
        """
        for parent in args:
            for option in self.available_options():
                value = option.getter(parent)
                if value is not None:
                    option.setter(self, value)
        for key, value in kwds.items():
            setattr(self, key, value)

    @classmethod
    def _add_option(cls, name, values, doc):
        """
        Add an option to the preferences system.

        This method should only be called during the import of
        :mod:`sage.repl.rich_output.preferences`.

        INPUT:

        - ``name`` -- the name of the option.

        - ``values`` -- the allowed values.

        - ``doc`` -- docstring.

        EXAMPLES::

            sage: from sage.repl.rich_output.preferences import PreferencesABC
            sage: class MyPrefs(PreferencesABC):
            ....:     pass
            sage: MyPrefs._add_option('foo', [0, 1, 2], 'The Foo Option')
            sage: prefs = MyPrefs()
            sage: prefs.foo
            sage: prefs.foo = 0
            sage: prefs.foo
            0
        """
        prop = Property(name, values, doc)
        setattr(cls, name, prop)

    def available_options(self):
        """
        Return the available options

        OUTPUT:

        Tuple of the preference items as instances of
        :class:`Property`.

        EXAMPLES::

            sage: from sage.repl.rich_output.preferences import DisplayPreferences
            sage: DisplayPreferences().available_options()
            (graphics, supplemental_plot, text)
        """
        options = []
        for key, value in self.__class__.__dict__.items():
            if isinstance(value, Property):
                options.append(value)
        return tuple(sorted(options, key=str))

    def _repr_(self):
        r"""
        Return a string representation

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.repl.rich_output.preferences import DisplayPreferences
            sage: DisplayPreferences()._repr_()
            'Display preferences:\n* graphics is not specified\n* supplemental_plot is not specified\n* text is not specified'
        """
        s = ['Display preferences:']
        for opt in self.available_options():
            value = opt.getter(self)
            if value is None:
                s += ['* {0} is not specified'.format(opt.name)]
            else:
                s += ['* {0} = {1}'.format(opt.name, value)]
        return '\n'.join(s)


class DisplayPreferences(PreferencesABC):
    pass


DisplayPreferences._add_option(
    'text',
    ('plain', 'ascii_art', 'unicode_art', 'latex'),
    """
    Which textual representation is preferred
    """
)


DisplayPreferences._add_option(
    'graphics',
    ('disable', 'vector', 'raster'),
    """
    Preferred graphics format
    """
)


DisplayPreferences._add_option(
    'supplemental_plot',
    ('always', 'never'),
    """
    Whether to graphically display graphs and other graph-like objects
    that implement rich output. When not specified small objects are
    show graphically and large objects as textual overview.
    """
)
