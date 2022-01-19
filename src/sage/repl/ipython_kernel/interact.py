# -*- coding: utf-8 -*-
r"""
Interacts for the Sage Jupyter notebook

This is mostly the same as the stock ``ipywidgets.interact``, but with
some customizations for Sage.

TESTS:

We need to setup a proper test environment for widgets::

    sage: from ipywidgets.widgets.tests.utils import setup_test_comm
    sage: setup_test_comm()

EXAMPLES::

    sage: from sage.repl.ipython_kernel.interact import interact
    sage: @interact
    ....: def f(x=(0,10)):
    ....:     pass
    Interactive function <function f at ...> with 1 widget
      x: IntSlider(value=5, description='x', max=10)
    sage: f.widget.children
    (IntSlider(value=5, description='x', max=10), Output())
"""

# ****************************************************************************
#       Copyright (C) 2017 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from ipywidgets.widgets import SelectionSlider, ValueWidget, ToggleButtons
from ipywidgets.widgets.interaction import interactive, signature
from collections import OrderedDict
from collections.abc import Iterable, Iterator
from .widgets import EvalText, SageColorPicker
from .widgets_sagenb import input_grid
from sage.structure.element import parent
import sage.rings.abc
from sage.plot.colors import Color
from sage.structure.element import Matrix


class sage_interactive(interactive):
    """
    Wrapper around the ipywidgets interactive which handles some SageNB
    specifics.

    EXAMPLES::

        sage: from sage.repl.ipython_kernel.interact import sage_interactive
        sage: def myfunc(x=10, y="hello", z=None): pass
        sage: sage_interactive(myfunc, x=(0,100), z=["one", "two", "three"])
        Interactive function <function myfunc at ...> with 3 widgets
          x: IntSlider(value=10, description='x')
          y: Text(value='hello', description='y')
          z: Dropdown(description='z', options=('one', 'two', 'three'), value=None)
    """
    def __init__(*args, **kwds):
        """
        See :class:`ipywidgets.widgets.interaction.interactive`

        TESTS::

            sage: from sage.repl.ipython_kernel.interact import sage_interactive
            sage: def myfunc(): pass
            sage: sage_interactive(myfunc, dict(manual=True))
            Manual interactive function <function myfunc ...> with 0 widgets

        ::

            sage: def myfunc(auto_update=False): pass
            sage: sage_interactive(myfunc)
            Manual interactive function <function myfunc ...> with 0 widgets
            sage: def myfunc(auto_update=None): pass
            sage: sage_interactive(myfunc)
            Interactive function <function myfunc ...> with 0 widgets
        """
        # Use *args to avoid name clash with keyword arguments
        if len(args) < 3:
            (self, f) = args
            options = {}
        else:
            (self, f, options) = args
            options = options.copy()

        # Check for auto_update in signature
        sig = signature(f)
        params = OrderedDict(sig.parameters)
        try:
            p_auto_update = params.pop("auto_update")
        except KeyError:
            pass
        else:
            options["manual"] = (p_auto_update.default is False)

        self.__signature = sig.replace(parameters=params.values())
        super(sage_interactive, self).__init__(f, options, **kwds)
        if self.manual:
            # In Sage, manual interacts are always run once
            self.on_displayed(self.update)
        else:
            # In automatic mode, clicking on a ToggleButtons button
            # should also run the interact
            for widget in self.kwargs_widgets:
                if isinstance(widget, ToggleButtons):
                    widget.on_msg(self.update)

    def __repr__(self):
        """
        Textual representation of this interactive function.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.interact import sage_interactive
            sage: def myfunc(): pass
            sage: sage_interactive(myfunc)
            Interactive function <function myfunc ...> with 0 widgets
        """
        s = "Manual interactive" if self.manual else "Interactive"
        widgets = [w for w in self.children if isinstance(w, ValueWidget)]
        n = len(widgets)
        s += " function %r with %s widget%s" % (self.f, n,
                                                "s" if n != 1 else "")
        for w in widgets:
            s += "\n  %s: %s" % (w._kwarg, w)
        return s

    def signature(self):
        """
        Return the fixed signature of the interactive function (after
        a possible ``auto_update`` parameter was removed).

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.interact import sage_interactive
            sage: def myfunc(x=[1,2,3], auto_update=False): pass
            sage: sage_interactive(myfunc).signature().parameters
            mappingproxy({'x': <Parameter "x=[1, 2, 3]">})
        """
        return self.__signature

    @classmethod  # Behaves like a staticmethod, but we need super()
    def widget_from_single_value(cls, abbrev, *args, **kwds):
        """
        Convert a single value (i.e. a non-iterable) to a widget.

        This supports the Sage :class:`Color` and ``Matrix`` classes.
        Any unknown type is changed to a string for evaluating.
        This is meant to support symbolic expressions like ``sin(x)``.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.interact import sage_interactive
            sage: sage_interactive.widget_from_single_value("sin(x)")
            Text(value='sin(x)')
            sage: sage_interactive.widget_from_single_value(sin(x))
            EvalText(value='sin(x)')
            sage: from sage.plot.colors import Color
            sage: sage_interactive.widget_from_single_value(matrix([[1, 2], [3, 4]]))
            Grid(value=[[1, 2], [3, 4]], children=(Label(value=''), VBox(children=(EvalText(value='1', layout=Layout(max_width='5em')), EvalText(value='3', layout=Layout(max_width='5em')))), VBox(children=(EvalText(value='2', layout=Layout(max_width='5em')), EvalText(value='4', layout=Layout(max_width='5em'))))))
            sage: sage_interactive.widget_from_single_value(Color('cornflowerblue'))
            SageColorPicker(value='#6495ed')
        """
        # Support Sage matrices and colors
        if isinstance(abbrev, Matrix):
            return input_grid(abbrev.nrows(), abbrev.ncols(),
                              default=abbrev.list(), to_value=abbrev.parent())
        if isinstance(abbrev, Color):
            return SageColorPicker(value=abbrev.html_color())
        # Get widget from IPython if possible
        widget = super(sage_interactive, cls).widget_from_single_value(abbrev, *args, **kwds)
        if widget is not None or isinstance(abbrev, Iterable):
            return widget
        # If IPython didn't construct a widget and the abbrev is not an
        # iterable, return an EvalText widget
        return EvalText(value=str(abbrev))

    @classmethod  # Behaves like a staticmethod, but we need super()
    def widget_from_tuple(cls, abbrev, *args, **kwds):
        """
        Convert a tuple to a widget.

        This supports two SageNB extensions: ``(description, abbrev)``
        if ``description`` is a string and ``(default, abbrev)`` if
        ``abbrev`` is not a single value.

        Symbolic expressions are changed to a floating-point number.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.interact import sage_interactive
            sage: sage_interactive.widget_from_tuple( (0, 10) )
            IntSlider(value=5, max=10)
            sage: sage_interactive.widget_from_tuple( ("number", (0, 10)) )
            IntSlider(value=5, description='number', max=10)
            sage: sage_interactive.widget_from_tuple( (3, (0, 10)) )
            IntSlider(value=3, max=10)
            sage: sage_interactive.widget_from_tuple((2, dict(one=1, two=2, three=3)))
            Dropdown(index=1, options={'one': 1, 'two': 2, 'three': 3}, value=2)
            sage: sage_interactive.widget_from_tuple( (sqrt(2), pi) )
            FloatSlider(value=2.277903107981444, max=3.141592653589793, min=1.4142135623730951)

        TESTS:

        Symbolic subrings::

            sage: SCR = SR.subring(no_variables=True)
            sage: sage_interactive.widget_from_tuple( (SCR(sqrt(2)), SCR(pi)) )
            FloatSlider(value=2.277903107981444, max=3.141592653589793, min=1.4142135623730951)
        """
        # Support (description, abbrev)
        if len(abbrev) == 2 and isinstance(abbrev[0], str):
            widget = cls.widget_from_abbrev(abbrev[1])
            widget.description = abbrev[0]
            return widget
        # Support (default, abbrev)
        if len(abbrev) == 2 and isinstance(abbrev[1], Iterable):
            widget = cls.widget_from_abbrev(abbrev[1])
            widget.value = abbrev[0]
            return widget
        # Numerically evaluate symbolic expressions

        def n(x):
            if isinstance(parent(x), sage.rings.abc.SymbolicRing):
                return x.numerical_approx()
            else:
                return x
        abbrev = tuple(n(x) for x in abbrev)
        return super(sage_interactive, cls).widget_from_tuple(abbrev, *args, **kwds)

    @classmethod  # Behaves like a staticmethod, but we need super()
    def widget_from_iterable(cls, abbrev, *args, **kwds):
        """
        Convert an unspecified iterable to a widget.

        This behaves like in ipywidgets, except that an iterator (like
        a generator object) becomes a ``SelectionSlider``.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.interact import sage_interactive
            sage: sage_interactive.widget_from_iterable([1..5])
            Dropdown(options=(1, 2, 3, 4, 5), value=1)
            sage: sage_interactive.widget_from_iterable(iter([1..5]))
            SelectionSlider(options=(1, 2, 3, 4, 5), value=1)
            sage: sage_interactive.widget_from_iterable((1..5))
            SelectionSlider(options=(1, 2, 3, 4, 5), value=1)
            sage: sage_interactive.widget_from_iterable(x for x in [1..5])
            SelectionSlider(options=(1, 2, 3, 4, 5), value=1)
            sage: def gen():
            ....:     yield 1; yield 2; yield 3; yield 4; yield 5
            sage: sage_interactive.widget_from_iterable(gen())
            SelectionSlider(options=(1, 2, 3, 4, 5), value=1)
        """
        if isinstance(abbrev, Iterator):
            return SelectionSlider(options=list(abbrev))
        return super(sage_interactive, cls).widget_from_iterable(abbrev, *args, **kwds)


# @interact decorator
interact = sage_interactive.factory()
