# sage_setup: distribution = sagemath-repl
r"""
Display Manager

This is the heart of the rich output system, the display manager
arbitrates between

* Backend capabilities: what can be displayed

* Backend preferences: what gives good quality on the backend

* Sage capabilities: every Sage object can only generate certain
  representations, and

* User preferences: typeset vs. plain text vs. ascii art, etc.

The display manager is a singleton class, Sage always has exactly one
instance of it. Use :func:`get_display_manager` to obtain it.

EXAMPLES::

    sage: from sage.repl.rich_output import get_display_manager
    sage: dm = get_display_manager();  dm
    The Sage display manager using the doctest backend
"""

# ****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


import warnings

from sage.structure.sage_object import SageObject
from sage.repl.rich_output.output_basic import (
    OutputPlainText, OutputAsciiArt, OutputUnicodeArt, OutputLatex,
)
from sage.repl.rich_output.output_browser import (
    OutputHtml,
)
from sage.repl.rich_output.preferences import DisplayPreferences


class DisplayException(Exception):
    """
    Base exception for all rich output-related exceptions.

    EXAMPLES::

        sage: from sage.repl.rich_output.display_manager import DisplayException
        sage: raise DisplayException('foo')
        Traceback (most recent call last):
        ...
        DisplayException: foo
    """
    pass


class OutputTypeException(DisplayException):
    """
    Wrong Output container.

    The output containers are the subclasses of
    :class:`~sage.repl.rich_output.output_basic.OutputBase` that
    contain the entire output. The display backends must create output
    containers of a suitable type depending on the displayed Python
    object. This exception indicates that there is a mistake in the
    backend and it returned the wrong type of output container.

    EXAMPLES::

        sage: from sage.repl.rich_output.display_manager import OutputTypeException
        sage: raise OutputTypeException('foo')
        Traceback (most recent call last):
        ...
        OutputTypeException: foo
    """
    pass


class RichReprWarning(UserWarning):
    """
    Warning that is throws if a call to ``_rich_repr_`` fails.

    If an object implements ``_rich_repr_`` then it must return a
    value, possibly ``None`` to indicate that no rich output can be
    generated. But it may not raise an exception as it is very
    confusing for the user if the displayhook fails.

    EXAMPLES::

        sage: from sage.repl.rich_output.display_manager import RichReprWarning
        sage: raise RichReprWarning('foo')
        Traceback (most recent call last):
        ...
        RichReprWarning: foo
    """
    pass


class restricted_output():

    def __init__(self, display_manager, output_classes):
        """
        Context manager to temporarily restrict the accepted output types.

        In the context, the output is restricted to the output
        container types listed in ``output_classes``. Additionally,
        display preferences are changed not to show graphics.

        INPUT:

        - ``display_manager`` -- the display manager

        - ``output_classes`` -- iterable of output container types

        EXAMPLES::

            sage: from sage.repl.rich_output.display_manager import (
            ....:     get_display_manager, restricted_output)
            sage: dm = get_display_manager()
            sage: restricted_output(dm, [dm.types.OutputPlainText])
            <sage.repl.rich_output.display_manager.restricted_output object at 0x...>
        """
        self._display_manager = display_manager
        self._output_classes = frozenset(output_classes)

    def __enter__(self):
        """
        Enter the restricted output context.

        EXAMPLES::

            sage: from sage.repl.rich_output.display_manager import (
            ....:     get_display_manager, restricted_output)
            sage: dm = get_display_manager()
            sage: len(dm.supported_output()) > 1
            True
            sage: with restricted_output(dm, [dm.types.OutputPlainText]):
            ....:    dm.supported_output()
            frozenset({<class 'sage.repl.rich_output.output_basic.OutputPlainText'>})

            sage: dm.preferences.supplemental_plot
            'never'
            sage: dm.preferences.supplemental_plot = 'always'
            sage: with restricted_output(dm, [dm.types.OutputPlainText]):
            ....:    dm.preferences
            Display preferences:
            * align_latex is not specified
            * graphics = disable
            * supplemental_plot = never
            * text is not specified
            sage: dm.preferences.supplemental_plot = 'never'
        """
        dm = self._display_manager
        self._original = dm._supported_output
        dm._supported_output = self._output_classes
        self._original_prefs = DisplayPreferences(dm.preferences)
        dm.preferences.graphics = 'disable'
        dm.preferences.supplemental_plot = 'never'

    def __exit__(self, exception_type, value, traceback):
        """
        Exit the restricted output context.

        EXAMPLES::

            sage: from sage.repl.rich_output.display_manager import (
            ....:     get_display_manager, restricted_output)
            sage: dm = get_display_manager()
            sage: with restricted_output(dm, [dm.types.OutputPlainText]):
            ....:     assert len(dm.supported_output()) == 1
            sage: assert len(dm.supported_output()) > 1
        """
        dm = self._display_manager
        dm._supported_output = self._original
        dm.preferences.graphics = self._original_prefs.graphics
        dm.preferences.supplemental_plot = self._original_prefs.supplemental_plot


class DisplayManager(SageObject):

    _instance = None

    def __init__(self):
        """
        The Display Manager.

        Used to decide what kind of rich output is best.

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: get_display_manager()
            The Sage display manager using the doctest backend
        """
        assert DisplayManager._instance is None
        DisplayManager._instance = self
        from sage.repl.rich_output.backend_base import BackendSimple
        self.switch_backend(BackendSimple())

    @classmethod
    def get_instance(cls):
        """
        Get the singleton instance.

        This class method is equivalent to
        :func:`get_display_manager`.

        OUTPUT: the display manager singleton

        EXAMPLES::

            sage: from sage.repl.rich_output.display_manager import DisplayManager
            sage: DisplayManager.get_instance()
            The Sage display manager using the doctest backend
        """
        if cls._instance is not None:
            return cls._instance
        else:
            return cls()

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: get_display_manager()
            The Sage display manager using the doctest backend
        """
        s = 'The Sage display manager using the {0} backend'.format(self._backend)
        return s

    @property
    def types(self):
        """
        Catalog of all output container types.

        Note that every output type must be registered in
        :mod:`sage.repl.rich_output.output_catalog`.

        OUTPUT:

        Returns the :mod:`sage.repl.rich_output.output_catalog`
        module.

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: dm.types.OutputPlainText
            <class 'sage.repl.rich_output.output_basic.OutputPlainText'>
        """
        import sage.repl.rich_output.output_catalog
        return sage.repl.rich_output.output_catalog

    def switch_backend(self, backend, **kwds):
        """
        Switch to a new backend.

        INPUT:

        - ``backend`` -- instance of
          :class:`~sage.repl.rich_output.backend_base.BackendBase`

        - ``kwds`` -- optional keyword arguments that are passed on to
          the :meth:`~sage.repl.rich_output.backend_base.BackendBase.install`
          method

        OUTPUT: the previous backend

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_base import BackendSimple
            sage: simple = BackendSimple()
            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager();  dm
            The Sage display manager using the doctest backend

            sage: previous = dm.switch_backend(simple)
            sage: dm
            The Sage display manager using the simple backend

        Restore the doctest backend::

            sage: dm.switch_backend(previous) is simple
            True
        """
        from sage.repl.rich_output.backend_base import BackendBase
        if not isinstance(backend, BackendBase):
            raise ValueError('backend must be instance of BackendBase class')
        supported = backend.supported_output()
        if not any(issubclass(out, OutputPlainText) for out in supported):
            raise ValueError('every backend must support plain text')
        try:
            self._backend.uninstall()
        except AttributeError:
            pass   # first time we switch
        # clear caches
        self._output_promotions = dict()
        self._supported_output = frozenset(
            map(self._demote_output_class, backend.supported_output()))
        # install new backend
        try:
            old_backend = self._backend
        except AttributeError:
            old_backend = None
        self._backend = backend
        self._backend.install(**kwds)
        self._preferences = DisplayPreferences(self._backend.default_preferences())
        return old_backend

    @property
    def preferences(self):
        """
        Return the preferences.

        OUTPUT:

        The display preferences as instance of
        :class:`~sage.repl.rich_output.preferences.DisplayPreferences`.

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: dm.preferences
            Display preferences:
            * align_latex is not specified
            * graphics is not specified
            * supplemental_plot = never
            * text is not specified
        """
        return self._preferences

    def is_in_terminal(self):
        """
        Test whether the UI is meant to run in a terminal.

        When this method returns ``True``, you can assume that it is
        possible to use ``raw_input`` or launch external programs that
        take over the input.

        Otherwise, you should assume that the backend runs remotely or
        in a pty controlled by another program. Then you should not
        launch external programs with a (text or graphical) UI.

        This is used to enable/disable interpreter consoles.

        OUTPUT:

        Boolean.
        """
        return self._backend.is_in_terminal()

    def check_backend_class(self, backend_class):
        """
        Check that the current backend is an instance of
        ``backend_class``.

        This is, for example, used by the Sage IPython display
        formatter to ensure that the IPython backend is in use.

        INPUT:

        - ``backend_class`` -- type of a backend class

        OUTPUT:

        This method returns nothing. A :exc:`RuntimeError` is raised if
        ``backend_class`` is not the type of the current backend.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_base import BackendSimple
            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: dm.check_backend_class(BackendSimple)
            Traceback (most recent call last):
            ...
            RuntimeError: check failed: current backend is invalid
        """
        if not isinstance(self._backend, backend_class):
            raise RuntimeError('check failed: current backend is invalid')

    def _demote_output_class(self, output_class):
        """
        Helper for :meth:`switch_backend`.

        INPUT:

        - ``output_class`` -- a possibly derived class from one of the
          output container classes in :meth:`types`

        OUTPUT: the underlying container class that it was derived from

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: dm._demote_output_class(dm.types.OutputPlainText)
            <class 'sage.repl.rich_output.output_basic.OutputPlainText'>
        """
        from sage.repl.rich_output.output_basic import OutputBase
        if not issubclass(output_class, OutputBase):
            raise OutputTypeException(
                'invalid output container type: {0} is not subclass of OutputBase'
                .format(output_class))
        result = None
        for type_name in dir(self.types):
            if type_name.startswith('_'):
                continue
            tp = getattr(self.types, type_name)
            if not issubclass(tp, OutputBase):
                continue
            if issubclass(output_class, tp):
                if result is not None:
                    raise OutputTypeException(
                        '{0} inherits from multiple output classes'
                        .format(output_class))
                else:
                    self._output_promotions[tp] = output_class
                    result = tp
        if result is None:
            raise OutputTypeException(
                '{0} does not inherit from any known output class'
                .format(output_class))
        return result

    def _promote_output(self, output):
        """
        Promote output container to a backend-specific subclass.

        INPUT:

        - ``output`` -- instance of a subclass of
          :class:`~sage.repl.rich_output.output_basic.OutputBase`. Backend-agnostic
          output container.

        OUTPUT:

        If the backend returned a subclass of the type of ``output``:
        promote the class to the subclass and return it. Otherwise,
        the unchanged output is returned.

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: out = dm._promote_output(dm.types.OutputPlainText('test'))
            sage: type(out)
            <class 'sage.repl.rich_output.output_basic.OutputPlainText'>
        """
        try:
            specialized_class = self._output_promotions[type(output)]
        except KeyError:
            return output
        output.__class__ = specialized_class
        return output

    def _preferred_text_formatter(self, obj, plain_text=None, **kwds):
        r"""
        Return the preferred textual representation.

        INPUT:

        - ``obj`` -- anything; the objects to format

        - ``plain_text`` -- ``None`` (default) or string. The plain
          text representation. If specified, this will be used for
          plain text output.

        OUTPUT:

        One of
        :class:`~sage.repl.rich_output.output_basic.OutputPlainText`,
        :class:`~sage.repl.rich_output.output_basic.OutputAsciiArt`,
        :class:`~sage.repl.rich_output.output_basic.OutputUnicodeArt`,
        or
        :class:`~sage.repl.rich_output.output_basic.OutputLatex`
        containing the preferred textual representation of ``obj``
        supported by the backend.

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: dm.preferences.text is None
            True
            sage: dm._preferred_text_formatter([1/42])
            OutputPlainText container

            sage: dm.preferences.text = 'plain'
            sage: dm._preferred_text_formatter([1/42])
            OutputPlainText container

            sage: dm.preferences.text = 'ascii_art'
            sage: dm._preferred_text_formatter([1/42])
            OutputAsciiArt container

            sage: dm.preferences.text = 'unicode_art'
            sage: dm._preferred_text_formatter([1/42])
            OutputUnicodeArt container

            sage: dm.preferences.text = 'latex'
            sage: dm._preferred_text_formatter([1/42])  # doctest backend does not support html
            OutputPlainText container

            sage: del dm.preferences.text   # reset to default
        """
        want = self.preferences.text
        supported = self._backend.supported_output()
        if want == 'ascii_art' and OutputAsciiArt in supported:
            out = self._backend.ascii_art_formatter(obj, **kwds)
            if type(out) is not OutputAsciiArt:
                raise OutputTypeException('backend returned wrong output type, require AsciiArt')
            return out
        if want == 'unicode_art' and OutputUnicodeArt in supported:
            out = self._backend.unicode_art_formatter(obj, **kwds)
            if type(out) is not OutputUnicodeArt:
                raise OutputTypeException('backend returned wrong output type, require UnicodeArt')
            return out
        if want == 'latex' and OutputHtml in supported:
            out = self._backend.latex_formatter(obj, **kwds)
            if type(out) is not OutputHtml:
                raise OutputTypeException('backend returned wrong output type, require Html')
            return out
        if plain_text is not None:
            if type(plain_text) is not OutputPlainText:
                raise OutputTypeException('backend returned wrong output type, require PlainText')
            return plain_text
        out = self._backend.plain_text_formatter(obj, **kwds)
        if type(out) is not OutputPlainText:
            raise OutputTypeException('backend returned wrong output type, require PlainText')
        return out

    def _call_rich_repr(self, obj, rich_repr_kwds):
        """
        Call the ``_rich_repr_`` method.

        This method calls ``obj._rich_repr_``. If this raises an
        exception, it is caught and a suitable warning is displayed.

        INPUT:

        - ``obj`` -- anything

        - ``rich_repr_kwds`` -- dictionary; optional keyword arguments
          that are passed through to ``obj._rich_repr_``

        OUTPUT:

        Whatever ``_rich_repr_`` returned. If it raises an exception,
        then a :class:`DisplayFormatterWarning` is displayed and
        ``None`` is returned.

        EXAMPLES::

            sage: class Foo(SageObject):
            ....:     def _rich_repr_(self, display_manager, **kwds):
            ....:         raise ValueError('reason')

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: dm._call_rich_repr(Foo(), {})
            doctest:...: ...RichReprWarning: Exception in _rich_repr_ while displaying object: reason
        """
        if rich_repr_kwds:
            # do not ignore errors from invalid options
            return obj._rich_repr_(self, **rich_repr_kwds)
        try:
            return obj._rich_repr_(self)
        except NotImplementedError:
            # don't warn on NotImplementedErrors
            return None
        except Exception as e:
            warnings.warn(
                'Exception in _rich_repr_ while displaying object: {0}'.format(e),
                RichReprWarning,
            )

    def _rich_output_formatter(self, obj, rich_repr_kwds):
        """
        Generate appropriate rich output containers.

        INPUT:

        - ``obj`` -- anything

        - ``rich_repr_kwds`` -- dictionary; optional keyword arguments
          that are passed through to ``obj._rich_repr_``

        OUTPUT:

        A pair of rich output containers. The first is plain text,
        that is, an instance of
        :class:`~sage.repl.rich_output.output_basic.OutputPlainText`. The
        second is the best rich output available, subject to the
        constraints of the backend.

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: dm._rich_output_formatter(1/2, dict())
            (OutputPlainText container, OutputPlainText container)
        """
        rich_output = None
        plain_text = None
        has_rich_repr = isinstance(obj, SageObject) and hasattr(obj, '_rich_repr_')
        if has_rich_repr:
            rich_output = self._call_rich_repr(obj, rich_repr_kwds)
        if isinstance(rich_output, OutputPlainText):
            plain_text = rich_output
        elif has_rich_repr:
            with restricted_output(self, [OutputPlainText]):
                plain_text = self._call_rich_repr(obj, rich_repr_kwds)
        if plain_text is None:
            plain_text = self._backend.plain_text_formatter(obj, **rich_repr_kwds)
        if rich_output is None:
            rich_output = self._preferred_text_formatter(
                obj, plain_text=plain_text, **rich_repr_kwds)
        # promote output container types to backend-specific containers
        plain_text = self._promote_output(plain_text)
        rich_output = self._promote_output(rich_output)
        # check that the output container types are valid for the backend
        supported = self._backend.supported_output()
        if type(plain_text) not in supported:
            raise OutputTypeException(
                'text output container not supported: {0}'.format(type(plain_text)))
        if type(rich_output) not in supported:
            raise OutputTypeException(
                'output container not supported: {0}'.format(type(rich_output)))
        return plain_text, rich_output

    def graphics_from_save(self, save_function, save_kwds,
                           file_extension, output_container,
                           figsize=None, dpi=None):
        r"""
        Helper to construct graphics.

        This method can be used to simplify the implementation of a
        ``_rich_repr_`` method of a graphics object if there is
        already a function to save graphics to a file.

        INPUT:

        - ``save_function`` -- callable that can save graphics to a file
          and accepts options like :meth:`sage.plot.graphics.Graphics.save`

        - ``save_kwds`` -- dictionary; keyword arguments that are
          passed to the save function

        - ``file_extension`` -- string starting with ``'.'``; the file
          extension of the graphics file

        - ``output_container`` -- subclass of
          :class:`sage.repl.rich_output.output_basic.OutputBase`. The
          output container to use. Must be one of the types in
          :meth:`supported_output`.

        - ``figsize`` -- pair of integers (optional). The desired graphics
          size in pixels. Suggested, but need not be respected by the
          output.

        - ``dpi`` -- integer (optional). The desired resolution in dots
          per inch. Suggested, but need not be respected by the output.

        OUTPUT: an instance of ``output_container``

        EXAMPLES::

            sage: # needs sage.plot sage.symbolic
            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: plt = plot(sin)
            sage: out = dm.graphics_from_save(plt.save, dict(), '.png',
            ....:                             dm.types.OutputImagePng)
            sage: out
            OutputImagePng container
            sage: out.png.get().startswith(b'\x89PNG')
            True
            sage: out.png.filename()   # random
            '/home/user/.sage/temp/localhost.localdomain/23903/tmp_pu5woK.png'
        """
        import os
        if not file_extension.startswith(os.path.extsep):
            raise ValueError('file_extension must start with a period')
        if output_container not in self.supported_output():
            raise OutputTypeException('output_container is not supported by backend')
        from sage.misc.temporary_file import tmp_filename
        filename = tmp_filename(ext=file_extension)
        # Call the save_function with the right arguments
        kwds = dict(save_kwds)
        if figsize is not None:
            kwds['figsize'] = figsize
        if dpi is not None:
            kwds['dpi'] = dpi
        save_function(filename, **kwds)
        from sage.repl.rich_output.buffer import OutputBuffer
        buf = OutputBuffer.from_file(filename)
        return output_container(buf)

    def threejs_scripts(self, online):
        """
        Return Three.js script tag for the current backend.

        INPUT:

        - ``online`` -- boolean determining script usage context

        OUTPUT: string containing script tag

        .. NOTE::

            This base method handles ``online=True`` case only, serving CDN
            script tag. Location of script for offline usage is
            backend-specific.

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: get_display_manager().threejs_scripts(online=True)                    # needs sage.plot
            '...<script src="https://cdn.jsdelivr.net/gh/sagemath/threejs-sage@...'
            sage: get_display_manager().threejs_scripts(online=False)                   # needs sage.plot
            Traceback (most recent call last):
            ...
            ValueError: current backend does not support
            offline threejs graphics
        """
        from sage.features.threejs import Threejs
        if online:
            version = Threejs().required_version()
            return """
<script src="https://cdn.jsdelivr.net/gh/sagemath/threejs-sage@{0}/build/three.min.js"></script>
            """.format(version)
        try:
            return self._backend.threejs_offline_scripts()
        except AttributeError:
            raise ValueError(
                'current backend does not support offline threejs graphics')

    def supported_output(self):
        """
        Return the output container classes that can be used.

        OUTPUT:

        Frozen set of subclasses of
        :class:`~sage.repl.rich_output.output_basic.OutputBase`. If
        the backend defines derived container classes, this method
        will always return their base classes.

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: dm.types.OutputPlainText in dm.supported_output()
            True
            sage: type(dm.supported_output())
            <... 'frozenset'>
        """
        return self._supported_output

    def displayhook(self, obj):
        """
        Implementation of the displayhook.

        Every backend must pass the value of the last statement of a
        line / cell to this method. See also
        :meth:`display_immediately` if you want do display rich output
        while a program is running.

        INPUT:

        - ``obj`` -- anything; the object to be shown

        OUTPUT: whatever the backend's displayhook method returned

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: dm.displayhook(1/2)
            1/2
        """
        if obj is None:
            return
        self._backend.set_underscore_variable(obj)
        plain_text, rich_output = self._rich_output_formatter(obj, dict())
        return self._backend.displayhook(plain_text, rich_output)

    def display_immediately(self, obj, **rich_repr_kwds):
        """
        Show output without going back to the command line prompt.

        This method must be called to create rich output from an
        object when we are not returning to the command line prompt,
        for example during program execution. Typically, it is being
        called by :meth:`sage.plot.graphics.Graphics.show`.

        INPUT:

        - ``obj`` -- anything; the object to be shown

        - ``rich_repr_kwds`` -- optional keyword arguments that are
          passed through to ``obj._rich_repr_``

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: dm.display_immediately(1/2)
            1/2
        """
        plain_text, rich_output = self._rich_output_formatter(obj, rich_repr_kwds)
        self._backend.display_immediately(plain_text, rich_output)


get_display_manager = DisplayManager.get_instance
