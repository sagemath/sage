"""
Magics for each of the Sage interfaces

This module defines magic functions for interpreters. As an example,
consider the GAP interpreter which can evaluate a gap command given as
a string::

    sage: gap('SymmetricGroup(4)')        # not tested
    SymmetricGroup( [ 1 .. 4 ] )

Magics are syntactic sugar to avoid writing the Python string. They
are either called as line magics::

    sage: %gap SymmetricGroup(4)          # not tested

or as cell magics, that is, spanning multiple lines::

    sage: %%gap                           # not tested
    ....: G := SymmetricGroup(4);
    ....: Display(G);

Note that the cell magic needs semicolons, this is required by the GAP
language to separate multiple commands.
"""
# Note: no magics in doctests, hence # not tested


# ****************************************************************************
#       Copyright (C) 2016 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.repl.rich_output.display_manager import get_display_manager


LINE_DOCSTRING = """
Interact with {name}

The line magic %{name} sends a single line to the {name} interface.

INPUT:

Single {name} command

OUTPUT:

The result of the evaluated {name} command as an interface wrapper
object.

EXAMPLES::

    sage: %{name} 1 + 2 + 3     # not tested
"""


CELL_DOCSTRING = """
Interact with {name}

The cell magic %%{name} sends multiple lines to the {name} interface.

INPUT:

Multiple lines of valid {name}-commands

OUTPUT:

The result of the evaluated {name}-commands is printed.

EXAMPLES::

    sage: %%{name}                    # not tested
    ....: 1 + 2 + 3;
    ....: some_{name}_command();
"""


class InterfaceMagic(object):

    @classmethod
    def all_iter(cls):
        """
        Iterate over the available interfaces

        EXAMPLES::

            sage: from sage.repl.interface_magic import InterfaceMagic
            sage: next(InterfaceMagic.all_iter())
            <sage.repl.interface_magic.InterfaceMagic object at 0x...>
        """
        import sage.interfaces.all
        for name, obj in sage.interfaces.all.__dict__.items():
            if isinstance(obj, sage.interfaces.interface.Interface):
                yield cls(name, obj)

    @classmethod
    def register_all(cls, shell=None):
        """
        Register all available interfaces

        EXAMPLES::

            sage: class MockShell():
            ....:     magics = set()
            ....:     def register_magic_function(self, fn, magic_name, magic_kind):
            ....:         self.magics.add(magic_name)
            ....:         print(magic_name, magic_kind)
            sage: from sage.repl.interface_magic import InterfaceMagic
            sage: InterfaceMagic.register_all(MockShell())    # random output
            ('gp', 'line')
            ('gp', 'cell')
            ('mwrank', 'line')
            ('mwrank', 'cell')
            ...
            ('maxima', 'line')
            ('maxima', 'cell')
            sage: 'gap' in MockShell.magics
            True
            sage: 'maxima' in MockShell.magics
            True
        """
        if shell is None:
            shell = get_ipython()
        for interface in cls.all_iter():
            shell.register_magic_function(
                interface.line_magic_factory(),
                magic_name=interface._name,
                magic_kind='line'
            )
            shell.register_magic_function(
                interface.cell_magic_factory(),
                magic_name=interface._name,
                magic_kind='cell'
            )

    @classmethod
    def find(cls, name):
        """
        Find a particular magic by name

        This method is for doctesting purposes only.

        INPUT:

        - ``name`` -- string. The name of the interface magic to
          search for.

        OUTPUT:

        The corresponding :class:`InterfaceMagic` instance.

        EXAMPLES::

            sage: from sage.repl.interface_magic import InterfaceMagic
            sage: InterfaceMagic.find('gap')
            <sage.repl.interface_magic.InterfaceMagic object at 0x...>
        """
        for magic in cls.all_iter():
            if magic._name == name:
                return magic

    def __init__(self, name, interface):
        """
        Interface Magic

        This class is a wrapper around interface objects to provide
        them with magics.

        INPUT:

        - ``name`` -- string. The interface name

        - ``interface`` -- :class:`sage.interfaces.expect.Expect`. The
          interface to wrap.

        EXAMPLES::

            sage: from sage.repl.interface_magic import InterfaceMagic
            sage: InterfaceMagic.find('gap')
            <sage.repl.interface_magic.InterfaceMagic object at 0x...>
        """
        self._name = name
        self._interface = interface

    def line_magic_factory(self):
        """
        Factory for line magic

        OUTPUT:

        A function suitable to be used as line magic.

        EXAMPLES::

            sage: from sage.repl.interface_magic import InterfaceMagic
            sage: line_magic = InterfaceMagic.find('gap').line_magic_factory()
            sage: output = line_magic('1+1')
            sage: output
            2
            sage: type(output)
            <class 'sage.interfaces.gap.GapElement'>

        This is how the built line magic is used in practice::

            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.run_cell('%gap 1+1')
            2
            sage: shell.run_cell('%gap?')
            Docstring:
            Interact with gap
            <BLANKLINE>
            The line magic %gap sends a single line to the gap interface.
            ...
        """
        terminal = get_display_manager().is_in_terminal()

        def line_magic(line):
            if line:
                return self._interface(line)
            else:
                if terminal:
                    self._interface.interact()
                else:
                    raise SyntaxError('{0} command required'.format(self._name))
        line_magic.__doc__ = LINE_DOCSTRING.format(name=self._name)
        return line_magic

    def cell_magic_factory(self):
        r"""
        Factory for cell magic

        OUTPUT:

        A function suitable to be used as cell magic.

        EXAMPLES::

            sage: from sage.repl.interface_magic import InterfaceMagic
            sage: cell_magic = InterfaceMagic.find('gap').cell_magic_factory()
            sage: output = cell_magic('', '1+1;')
            2
            sage: output is None
            True
            sage: cell_magic('foo', '1+1;')
            Traceback (most recent call last):
            ...
            SyntaxError: Interface magics have no options, got "foo"

        This is how the built cell magic is used in practice::

            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.run_cell('%%gap\nG:=SymmetricGroup(5);\n1+1;Order(G);')
            Sym( [ 1 .. 5 ] )
            2
            120
            sage: shell.run_cell('%%gap foo\n1+1;\n')
            ...File "<string>", line unknown
            SyntaxError: Interface magics have no options, got "foo"
            <BLANKLINE>
            sage: shell.run_cell('%%gap?')
            Docstring:
            Interact with gap
            <BLANKLINE>
            The cell magic %%gap sends multiple lines to the gap interface.
            ...
        """
        def cell_magic(line, cell):
            """
            Evaluate cell magic

            Docstring is overwritten in the instance

            INPUT:

            - ``line`` -- string. The option part of the cell magic.

            - ``cell`` -- string. The lines of the cell magic.

            OUTPUT:

            Prints the interface output.

            RAISES:

            ``SyntaxError`` if a line is specified; Interfaces have no
            options.
            """
            if line:
                raise SyntaxError('Interface magics have no options, got "{0}"'.format(line))
            output = self._interface.eval(cell)
            print(output)
        cell_magic.__doc__ = CELL_DOCSTRING.format(name=self._name)
        return cell_magic
