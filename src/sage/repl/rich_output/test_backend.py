# sage_setup: distribution = sagemath-repl
r"""
Test Backend

This backend is only for doctesting purposes.

EXAMPLES:

We switch to the test backend for the remainder of this file::

    sage: from sage.repl.rich_output import get_display_manager
    sage: dm = get_display_manager()
    sage: from sage.repl.rich_output.test_backend import BackendTest, TestObject
    sage: doctest_backend = dm.switch_backend(BackendTest())
    sage: dm
    The Sage display manager using the test backend

    sage: dm._output_promotions
    {<class 'sage.repl.rich_output.output_basic.OutputPlainText'>:
     <class 'sage.repl.rich_output.test_backend.TestOutputPlainText'>}
    sage: dm.displayhook(1/2)
    1/2 [TestOutputPlainText]
    TestOutputPlainText container

    sage: test = TestObject()
    sage: test
    called the _repr_ method
    sage: dm.displayhook(test)
    called the _rich_repr_ method [TestOutputPlainText]
    TestOutputPlainText container
"""

# ****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.structure.sage_object import SageObject
from sage.repl.rich_output.backend_base import BackendBase
from sage.repl.rich_output.output_catalog import OutputPlainText, OutputImagePng


class TestOutputPlainText(OutputPlainText):

    def __init__(self, *args, **kwds):
        """
        Backend-specific subclass of the plain text output container.

        Backends must not influence how the display system constructs
        output containers, they can only control how the output
        container is displayed. In particular, we cannot override the
        constructor (only the
        :class:`~sage.repl.rich_output.output_basic.OutputPlainText`
        constructor is used).

        EXAMPLES::

            sage: from sage.repl.rich_output.test_backend import TestOutputPlainText
            sage: TestOutputPlainText()
            Traceback (most recent call last):
            ...
            AssertionError: cannot override constructor
        """
        raise AssertionError('cannot override constructor')

    def print_to_stdout(self):
        """
        Write the data to stdout.

        This is just a convenience method to help with debugging.

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: test_output = dm.displayhook(123)
            123 [TestOutputPlainText]
            sage: type(test_output)
            <class 'sage.repl.rich_output.test_backend.TestOutputPlainText'>
            sage: test_output.print_to_stdout()
            123 [TestOutputPlainText]
        """
        print('{0} [{1}]'.format(self.text.get_str(), self.__class__.__name__))


class TestObject(SageObject):
    """
    Test object with both :meth:`_repr_` and :meth:`_rich_repr_`
    """

    def _repr_(self):
        """
        Return string representation.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.repl.rich_output.test_backend import TestObject
            sage: obj = TestObject()
            sage: obj._repr_()
            'called the _repr_ method'
        """
        return 'called the _repr_ method'

    def _rich_repr_(self, display_manager):
        """
        Rich Output Magic Method.

        See :mod:`sage.repl.rich_output` for details.

        EXAMPLES::

            sage: display_manager = sage.repl.rich_output.get_display_manager()
            sage: from sage.repl.rich_output.test_backend import TestObject
            sage: obj = TestObject()
            sage: rich_output = obj._rich_repr_(display_manager);  rich_output
            OutputPlainText container
            sage: rich_output.text.get_str()
            'called the _rich_repr_ method'
        """
        tp = display_manager.types
        if tp.OutputPlainText in display_manager.supported_output():
            return tp.OutputPlainText('called the _rich_repr_ method')


class BackendTest(BackendBase):

    def _repr_(self):
        """
        Return the string representation.

        OUTPUT: string

        EXAMPLES::

            sage: display_manager = sage.repl.rich_output.get_display_manager()
            sage: backend = display_manager._backend
            sage: backend._repr_()
            'test'
        """
        return 'test'

    def supported_output(self):
        """
        Return the outputs that are supported by the backend.

        OUTPUT:

        Iterable of output container classes. Only the
        :class:`~sage.repl.rich_repr.test_backend.TestOutputPlainText`
        output container is supported by the test backend.

        EXAMPLES::

            sage: display_manager = sage.repl.rich_output.get_display_manager()
            sage: backend = display_manager._backend
            sage: list(backend.supported_output())
            [<class 'sage.repl.rich_output.test_backend.TestOutputPlainText'>]

        The output of this method is used by the display manager to
        set up the actual supported outputs. Compare::

            sage: list(display_manager.supported_output())
            [<class 'sage.repl.rich_output.output_basic.OutputPlainText'>]
        """
        return set([TestOutputPlainText])

    def display_immediately(self, plain_text, rich_output):
        """
        Show output without going back to the command line prompt.

        INPUT:

        Same as :meth:`displayhook`.

        OUTPUT:

        This method returns the rich output for doctesting
        convenience. The actual display framework ignores the return
        value.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_basic import OutputPlainText
            sage: plain_text = OutputPlainText.example()
            sage: from sage.repl.rich_output.test_backend import BackendTest
            sage: backend = BackendTest()
            sage: backend.display_immediately(plain_text, plain_text)
            Example plain text output
            OutputPlainText container
        """
        plain_text.print_to_stdout()
        if plain_text is not rich_output:
            print('rich output type: [{0}]'.format(rich_output.__class__.__name__))
        return rich_output
