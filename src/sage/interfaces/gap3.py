r"""
Interface to GAP3

This module implements an interface to GAP3.

AUTHORS:

-  Franco Saliola (February 2010)
-  Christian Stump (March 2016)

.. WARNING::

    The experimental package for GAP3 is Jean Michel's pre-packaged GAP3,
    which is a minimal GAP3 distribution containing packages that have
    no equivalent in GAP4, see :issue:`20107` and also

        https://webusers.imj-prg.fr/~jean.michel/gap3/

Obtaining GAP3
--------------

Instead of installing the experimental GAP3 package, one can as well install
by hand either of the following two versions of GAP3:

- Frank Luebeck maintains a GAP3 Linux executable, optimized
  for i686 and statically linked for jobs of 2 GByte or more:

    http://www.math.rwth-aachen.de/~Frank.Luebeck/gap/GAP3

- or you can download GAP3 from the GAP website below. Since GAP3
  is no longer supported, it may not be easy to install this version.

    https://www.gap-system.org/Gap3/Download3/download.html

Changing which GAP3 is used
---------------------------

.. WARNING::

    There is a bug in the pexpect module (see :issue:`8471`) that
    prevents the following from working correctly. For now, just make sure
    that ``gap3`` is in your ``PATH``.

Sage assumes that GAP3 can be launched with the command ``gap3``; that is,
Sage assumes that the command ``gap3`` is in your ``PATH``. If this is not
the case, then you can start GAP3 using the following command::

    sage: gap3 = Gap3(command='/usr/local/bin/gap3')            # not tested

Functionality and Examples
--------------------------

The interface to GAP3 offers the following functionality.

#.  ``gap3(expr)`` -- evaluation of arbitrary GAP3 expressions, with the
    result returned as a Sage object wrapping the corresponding GAP3 element::

        sage: # optional - gap3
        sage: a = gap3('3+2')
        sage: a
        5
        sage: type(a)
        <class 'sage.interfaces.gap3.GAP3Element'>

    ::

        sage: # optional - gap3
        sage: S5 = gap3('SymmetricGroup(5)')
        sage: S5
        Group( (1,5), (2,5), (3,5), (4,5) )
        sage: type(S5)
        <class 'sage.interfaces.gap3.GAP3Record'>

    This provides a Pythonic interface to GAP3. If ``gap_function`` is the
    name of a GAP3 function, then the syntax ``gap_element.gap_function()``
    returns the ``gap_element`` obtained by evaluating the command
    ``gap_function(gap_element)`` in GAP3::

        sage: # optional - gap3
        sage: S5.Size()
        120
        sage: S5.CharTable()
        CharTable( Group( (1,5), (2,5), (3,5), (4,5) ) )

    Alternatively, you can instead use the syntax
    ``gap3.gap_function(gap_element)``::

        sage: gap3.DerivedSeries(S5)                       #optional - gap3
        [ Group( (1,5), (2,5), (3,5), (4,5) ),
          Subgroup( Group( (1,5), (2,5), (3,5), (4,5) ),
                    [ (1,2,5), (1,3,5), (1,4,5) ] ) ]

    If ``gap_element`` corresponds to a GAP3 record, then
    ``gap_element.recfield`` provides a means to access the record element
    corresponding to the field ``recfield``::

        sage: # optional - gap3
        sage: S5.IsRec()
        true
        sage: S5.recfields()
        ['isDomain', 'isGroup', 'identity', 'generators', 'operations',
        'isPermGroup', 'isFinite', '1', '2', '3', '4', 'degree']
        sage: S5.identity
        ()
        sage: S5.degree
        5
        sage: S5.1
        (1,5)
        sage: S5.2
        (2,5)

#.  By typing ``%gap3`` or ``gap3.interact()`` at the command-line, you can
    interact directly with the underlying GAP3 session.

    ::

        sage: gap3.interact()                            # not tested

          --> Switching to Gap3 <--

        gap3:

#.  You can start a new GAP3 session as follows::

        sage: gap3.console()                             # not tested

                     ########            Lehrstuhl D fuer Mathematik
                   ###    ####           RWTH Aachen
                  ##         ##
                 ##          #             #######            #########
                ##                        #      ##          ## #     ##
                ##           #           #       ##             #      ##
                ####        ##           ##       #             #      ##
                 #####     ###           ##      ##             ##    ##
                   ######### #            #########             #######
                             #                                  #
                            ##           Version 3              #
                           ###           Release 4.4            #
                          ## #           18 Apr 97              #
                         ##  #
                        ##   #  Alice Niemeyer, Werner Nickel,  Martin Schoenert
                       ##    #  Johannes Meier, Alex Wegner,    Thomas Bischops
                      ##     #  Frank Celler,   Juergen Mnich,  Udo Polis
                      ###   ##  Thomas Breuer,  Goetz Pfeiffer, Hans U. Besche
                       ######   Volkmar Felsch, Heiko Theissen, Alexander Hulpke
                                Ansgar Kaup,    Akos Seress,    Erzsebet Horvath
                                Bettina Eick
                                For help enter: ?<return>
        gap>

#.  The interface also has access to the GAP3 help system::

        sage: gap3.help('help', pager=False)             # not tested
        Help _______________________________________________________...

        This  section describes  together with  the following sections the   GAP
        help system.  The help system lets you read the manual interactively...

Common Pitfalls
---------------

#.  If you want to pass a string to GAP3, then you need to wrap it in
    single quotes as follows::

        sage: gap3('"This is a GAP3 string"')              #optional - gap3
        "This is a GAP3 string"

    This is particularly important when a GAP3 package is loaded via the
    ``RequirePackage`` method (note that one can instead use the
    ``load_package`` method)::

        sage: gap3.RequirePackage('"chevie"')             #optional - gap3

Examples
--------

Load a GAP3 package::

    sage: # optional - gap3
    sage: gap3.load_package("chevie")
    sage: gap3.version() # random  # not tested
    'lib: v3r4p4 1997/04/18, src: v3r4p0 1994/07/10, sys: usg gcc ansi'

Working with GAP3 lists. Note that GAP3 lists are 1-indexed::

    sage: # optional - gap3
    sage: L = gap3([1,2,3])
    sage: L[1]
    1
    sage: L[2]
    2
    sage: 3 in L
    True
    sage: 4 in L
    False
    sage: m = gap3([[1,2],[3,4]])
    sage: m[2,1]
    3
    sage: [1,2] in m
    True
    sage: [3,2] in m
    False
    sage: gap3([1,2]) in m
    True

Controlling variable names used by GAP3::

    sage: # optional - gap3
    sage: gap3('2', name='x')
    2
    sage: gap3('x')
    2
    sage: gap3.unbind('x')
    sage: gap3('x')
    Traceback (most recent call last):
    ...
    TypeError: Gap3 produced error output
    Error, Variable: 'x' must have a value
    ...
"""

# ****************************************************************************
#       Copyright (C) 2010 Franco Saliola <saliola@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import os

from sage.misc.cachefunc import cached_method
from sage.interfaces.expect import Expect
from sage.interfaces.gap import Gap_generic, GapElement_generic
from sage.cpython.string import bytes_to_str


# gap3_cmd should point to the gap3 executable
gap3_cmd = 'gap3'


class Gap3(Gap_generic):
    r"""
    A simple Expect interface to GAP3.

    EXAMPLES::

        sage: from sage.interfaces.gap3 import Gap3
        sage: gap3 = Gap3(command='gap3')

    TESTS::

        sage: gap3(2) == gap3(3)                           #optional - gap3
        False
        sage: gap3(2) == gap3(2)                           #optional - gap3
        True
        sage: gap3._tab_completion()                       #optional - gap3
        []

    We test the interface behaves correctly after a keyboard interrupt::

        sage: gap3(2)                                      #optional - gap3
        2
        sage: try:
        ....:     gap3._keyboard_interrupt()
        ....: except KeyboardInterrupt:
        ....:     pass
        sage: gap3(2)                                      #optional - gap3
        2

    We test that the interface busts out of GAP3's break loop correctly::

        sage: f = gap3('function(L) return L[0]; end;;')   #optional - gap3
        sage: f([1,2,3])                                   #optional - gap3
        Traceback (most recent call last):
        ...
        RuntimeError: Gap3 produced error output
        Error, List Element: <position> must be a positive integer at
        return L[0] ...

    AUTHORS:

    - Franco Saliola (Feb 2010)
    """
    _identical_function = "IsIdentical"

    def __init__(self, command=gap3_cmd):
        r"""
        Initialize the GAP3 interface and start a session.

        INPUT:

        - command -- string (default: ``'gap3'``); points to the gap3
          executable on your system. By default, it is assumed the
          executable is in your path.

        EXAMPLES::

            sage: # optional - gap3
            sage: gap3 = Gap3()
            sage: gap3.is_running()
            False
            sage: gap3._start()
            sage: gap3.is_running()
            True
        """
        self.__gap3_command_string = command
        # Explanation of additional command-line options passed to gap3:
        #
        #     -p invokes the internal programmatic interface, which is how Sage
        #     talks to GAP4. This allows reuse some of the GAP4 interface code.
        #
        #     -y -- sets the number of lines of the terminal; controls how many
        #     lines of text are output by GAP3 before the pager is invoked.
        #     This option is useful in dealing with the GAP3 help system.
        Expect.__init__(self,
                        name='gap3',
                        prompt='gap> ',
                        command=self.__gap3_command_string + " -p -b -y 500",
                        server=None,
                        ulimit=None,
                        script_subdirectory=None,
                        restart_on_ctrlc=True,
                        verbose_start=False,
                        init_code=[],
                        max_startup_time=None,
                        logfile=None,
                        eval_using_file_cutoff=100,
                        do_cleaner=True,
                        remote_cleaner=False,
                        path=None)

    def _start(self):
        r"""
        Initialize the interface and start gap3.

        EXAMPLES::

            sage: # optional - gap3
            sage: gap3 = Gap3()
            sage: gap3.is_running()
            False
            sage: gap3._start()
            sage: gap3.is_running()
            True

        Check that :issue:`23142` is fixed::

            sage: # optional - gap3
            sage: gap3.eval("1+1")
            '2'
            sage: gap3.quit()
        """
        Expect._start(self)
        # The -p command-line option to GAP3 produces the following
        # funny-looking patterns in the interface. We compile the patterns
        # now, and use them later for interpreting interface messages.
        self._compiled_full_pattern = self._expect.compile_pattern_list([
            r'@p\d+\.', '@@', '@[A-Z]', r'@[123456!"#$%&][^+]*\+', '@e', '@c',
            '@f', '@h', '@i', '@m', '@n', '@r', r'@s\d', r'@w.*\+', '@x', '@z'])
        self._compiled_small_pattern = self._expect.compile_pattern_list('@J')
        self._expect.expect("@i")

    def _object_class(self):
        r"""
        Return the class used for constructing GAP3 elements.

        TESTS::

            sage: gap3._object_class()
            <class 'sage.interfaces.gap3.GAP3Element'>
        """
        return GAP3Element

    def _execute_line(self, line, wait_for_prompt=True, expect_eof=False):
        r"""
        Execute a line of code in GAP3 and parse the output.

        TESTS:

        Test that syntax errors are handled correctly::

            sage: gap3('syntax error', name='x')               #optional - gap3
            Traceback (most recent call last):
            ...
            TypeError: Gap3 produced error output
            Syntax error: ; expected
            x:=syntax error;
                          ^
            gap>
            ...

        Test that error messages are detected and reported correctly::

            sage: gap3.SymmetricGroup(5,3)                     #optional - gap3
            Traceback (most recent call last):
            ...
            RuntimeError: Gap3 produced error output
            Error, Record: left operand must be a record at
            return arg[1].operations.SymmetricGroup( arg[1], arg[2] ) ... in
            SymmetricGroup( ..., ... ) called from
            main loop
            brk> quit;
            ...

        Test that break loops are detected and exited properly::

            sage: f = gap3('function(L) return L[0]; end;;')   #optional - gap3
            sage: f([1,2,3])                                   #optional - gap3
            Traceback (most recent call last):
            ...
            RuntimeError: Gap3 produced error output
            Error, List Element: <position> must be a positive integer at
            return L[0] ... in
            ... called from
            main loop
            brk> quit;
            ...
        """
        # It seems that GAP3 does not classify syntax errors as regular error
        # messages, so the generic GAP interface processing code does not
        # detect it. So we test for a syntax error explicitly.
        normal_output, error_output = \
            super()._execute_line(line, wait_for_prompt=True, expect_eof=False)
        normal = bytes_to_str(normal_output)
        if normal.startswith("Syntax error:"):
            normal_output, error_output = "", normal_output
        return (normal_output, error_output)

    def help(self, topic, pager=True):
        r"""
        Print help on the given topic.

        INPUT:

        - ``topic`` -- string

        EXAMPLES::

            sage: gap3.help('help', pager=False)           #optional - gap3
            Help _______________________________________________________...
            <BLANKLINE>
            This  section describes  together with  the following sectio...
            help system.  The help system lets you read the manual inter...

        ::

            sage: gap3.help('SymmetricGroup', pager=False) #optional - gap3
            no section with this name was found

        TESTS::

            sage: # optional - gap3
            sage: m = gap3([[1,2,3],[4,5,6]]); m
            [ [ 1, 2, 3 ], [ 4, 5, 6 ] ]
            sage: gap3.help('help', pager=False)
            Help _______________________________________________________...
            sage: m
            [ [ 1, 2, 3 ], [ 4, 5, 6 ] ]
            sage: m.Print()
            [ [ 1, 2, 3 ], [ 4, 5, 6 ] ]
            sage: gap3.help('Group', pager=False)
            Group ______________________________________________________...
            sage: m
            [ [ 1, 2, 3 ], [ 4, 5, 6 ] ]
            sage: m.Print()
            [ [ 1, 2, 3 ], [ 4, 5, 6 ] ]
        """

        import pexpect
        if self._expect is None:
            self._start()
        E = self._expect
        helptext = []

        # we request the help document
        E.sendline("? %s" % topic)
        # it seems necessary to skip TWO echoes (as is done in the GAP4 interface)
        E.expect("\r\n")
        E.expect("\r\n")

        # next we process the help document; translating special characters, etc.
        while True:
            try:
                x = E.expect_list(self._compiled_full_pattern, timeout=2)
            except pexpect.TIMEOUT:
                break
            if x == 1:
                # matched @@; replace with @
                helptext.append('@')
                helptext.append(E.before)
            elif x == 2:
                # matched a special char; convert and insert
                helptext.append(chr(ord(E.after[1:2]) - ord('A') + 1))
                helptext.append(E.before)
            elif x == 10:
                # matched @n (normal input mode); it seems we're done
                break
            elif x == 11:
                # matched @r (echoing input); skip to end of line
                E.expect_list(self._compiled_small_pattern)

        # merge the help text into one string and print it.
        helptext = "".join(bytes_to_str(line) for line in helptext).strip()
        if pager is True:
            from sage.misc.pager import pager as pag
            pag()(helptext)
        else:
            print(helptext)

    def cputime(self, t=None):
        r"""
        Return the amount of CPU time that the GAP session has used in seconds.

        If ``t`` is not None, then it returns the difference
        between the current CPU time and ``t``.

        EXAMPLES::

            sage: # optional - gap3
            sage: t = gap3.cputime()
            sage: t  # random
            0.02
            sage: gap3.SymmetricGroup(5).Size()
            120
            sage: gap3.cputime()  # random
            0.14999999999999999
            sage: gap3.cputime(t)  # random
            0.13
        """
        if t is not None:
            return self.cputime() - t
        else:
            return eval(self.eval('Runtime();')) / 1000.0

    def console(self):
        r"""
        Spawn a new GAP3 command-line session.

        EXAMPLES::

            sage: gap3.console()                             # not tested

                         ########            Lehrstuhl D fuer Mathematik
                       ###    ####           RWTH Aachen
                      ##         ##
                     ##          #             #######            #########
                    ##                        #      ##          ## #     ##
                    ##           #           #       ##             #      ##
                    ####        ##           ##       #             #      ##
                     #####     ###           ##      ##             ##    ##
                       ######### #            #########             #######
                                 #                                  #
                                ##           Version 3              #
                               ###           Release 4.4            #
                              ## #           18 Apr 97              #
                             ##  #
                            ##   #  Alice Niemeyer, Werner Nickel,  Martin Schoenert
                           ##    #  Johannes Meier, Alex Wegner,    Thomas Bischops
                          ##     #  Frank Celler,   Juergen Mnich,  Udo Polis
                          ###   ##  Thomas Breuer,  Goetz Pfeiffer, Hans U. Besche
                           ######   Volkmar Felsch, Heiko Theissen, Alexander Hulpke
                                    Ansgar Kaup,    Akos Seress,    Erzsebet Horvath
                                    Bettina Eick
                                    For help enter: ?<return>
            gap>
        """
        os.system(self.__gap3_command_string)

    def _install_hints(self):
        r"""

        TESTS::

            sage: gap3 = Gap3(command='/wrongpath/gap3')
            sage: gap3('3+2')
            Traceback (most recent call last):
            ...
            TypeError: unable to start gap3 because the command '/wrongpath/gap3 ...' failed: The command was not found or was not executable: /wrongpath/gap3.
            <BLANKLINE>
                Your attempt to start GAP3 failed, either because you do not have
                have GAP3 installed, or because it is not configured correctly.
            <BLANKLINE>
                - If you do not have GAP3 installed, then you must either...
            sage: print(gap3._install_hints())
            <BLANKLINE>
                Your attempt to start GAP3 failed, either because you do not have
                have GAP3 installed, or because it is not configured correctly.
            <BLANKLINE>
                - If you do not have GAP3 installed, then you must either...
        """
        return r"""
    Your attempt to start GAP3 failed, either because you do not have
    have GAP3 installed, or because it is not configured correctly.

    - If you do not have GAP3 installed, then you must either install
      the optional package, see :issue:`20107`, or you download and
      install it yourself.
      Here are two other ways to obtain GAP3:

        - Frank Luebeck maintains a GAP3 Linux executable, optimized
          for i686 and statically linked for jobs of 2 GByte or more:
            http://www.math.rwth-aachen.de/~Frank.Luebeck/gap/GAP3

        - Finally, you can download GAP3 from the GAP website below. Since
          GAP3 is no longer an officially supported distribution of GAP, it
          may not be easy to install this version.
            https://www.gap-system.org/Gap3/Download3/download.html

    - If you have GAP3 installed, then perhaps it is not configured
      correctly. Sage assumes that you can start GAP3 with the command
      %s. Alternatively, you can use the following command
      to point Sage to the correct command for your system.

          gap3 = Gap3(command='/usr/local/bin/gap3')
        """ % self.__gap3_command_string

    @cached_method
    def _tab_completion(self):
        """
        Return additional tab completion entries.

        Currently this is empty.

        OUTPUT: list of strings

        EXAMPLES::

            sage: gap3._tab_completion()
            []
        """
        return []


gap3 = Gap3()


class GAP3Element(GapElement_generic):
    r"""
    A GAP3 element.

    .. NOTE::

        If the corresponding GAP3 element is a GAP3 record,
        then the class is changed to a ``GAP3Record``.

    INPUT:

    - ``parent`` -- the GAP3 session

    - ``value`` -- the GAP3 command as a string

    - ``is_name`` -- boolean (default: ``False``); if ``True``, then ``value`` is
      the variable name for the object

    - ``name`` -- string (default: ``None``); the variable name to use for the
      object. If ``None``, then a variable name is generated

    .. NOTE::

        If you pass ``E``, ``X`` or ``Z`` for ``name``, then an error is
        raised because these are sacred variable names in GAP3 that should
        never be redefined. Sage raises an error because GAP3 does not!

    EXAMPLES::

        sage: # optional - gap3
        sage: from sage.interfaces.gap3 import GAP3Element
        sage: gap3 = Gap3()
        sage: GAP3Element(gap3, value='3+2')
        5
        sage: GAP3Element(gap3, value='sage0', is_name=True)
        5

    TESTS::

        sage: GAP3Element(gap3, value='3+2', is_name=False, name='X') #optional - gap3
        Traceback (most recent call last):
        ...
        ValueError: you are attempting to redefine X; but you should never redefine E, X or Z in gap3 (because things will break!)

    AUTHORS:

    - Franco Saliola (Feb 2010)
    """
    def __init__(self, parent, value, is_name=False, name=None):
        r"""
        See ``GAP3Element`` for full documentation.

        EXAMPLES::

            sage: # optional - gap3
            sage: from sage.interfaces.gap3 import GAP3Element
            sage: gap3 = Gap3()
            sage: GAP3Element(gap3, value='3+2')
            5
            sage: GAP3Element(gap3, value='sage0', is_name=True)
            5

        TESTS::

            sage: GAP3Element(gap3, value='3+2', is_name=False, name='X') #optional - gap3
            Traceback (most recent call last):
            ...
            ValueError: you are attempting to redefine X; but you should never redefine E, X or Z in gap3 (because things will break!)
        """
        # Warning: One should not redefine E, X or Z in gap3, because
        # things will break, but gap3 raises no errors if one does this!
        if name in ["E", "X", "Z"]:
            raise ValueError("you are attempting to redefine %s; but you should never redefine E, X or Z in gap3 (because things will break!)" % name)

        # initialize the superclass
        super().__init__(parent, value, is_name, name)

        # check for a GAP record; if so then change the class
        parent._synchronize()
        if parent.eval("IsRec(%s)" % self._name) == "true":
            self.__class__ = GAP3Record

    def __getitem__(self, n):
        r"""
        EXAMPLES::

            sage: # optional - gap3
            sage: l = gap3('[1,2,3]')
            sage: l[1]
            1
            sage: a = gap3([1,2,3])
            sage: a[1]
            1
            sage: m = gap3([[1,2,3],[4,5,6],[7,8,9]])
            sage: m[1,3]
            3
            sage: m[2][1]
            4
        """
        gap3_session = self._check_valid()
        if not isinstance(n, tuple):
            return gap3_session.new('%s[%s]' % (self.name(), n))
        return gap3_session.new('%s%s' % (self.name(),
                                          ''.join('[%s]' % x for x in n)))

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: # optional - gap3
            sage: s = gap3("[[1,2], [3/4, 5/6]]")
            sage: s._latex_()
            '\\left(\\begin{array}{rr} 1&2\\\\ 3/4&\\frac{5}{6}\\\\ \\end{array}\\right)'
            sage: latex(s)
            \left(\begin{array}{rr} 1&2\\ 3/4&\frac{5}{6}\\ \end{array}\right)
        """
        gap3_session = self._check_valid()
        try:
            s = gap3_session.eval('FormatLaTeX(%s)' % self.name())
            s = s.replace('\\\\', '\\').replace('"', '')
            s = s.replace('%\\n', ' ')
            return s
        except RuntimeError:
            return str(self)


class GAP3Record(GAP3Element):
    r"""
    A GAP3 record.

    .. NOTE::

        This class should not be called directly, use GAP3Element instead.
        If the corresponding GAP3 element is a GAP3 record, then the class
        is changed to a ``GAP3Record``.

    AUTHORS:

    - Franco Saliola (Feb 2010)
    """
    def recfields(self):
        r"""
        Return a list of the fields for the record. (Record fields are akin
        to object attributes in Sage.)

        OUTPUT: list of strings -- the field records

        EXAMPLES::

            sage: S5 = gap3.SymmetricGroup(5)              #optional - gap3
            sage: S5.recfields()                           #optional - gap3
            ['isDomain', 'isGroup', 'identity', 'generators',
             'operations', 'isPermGroup', 'isFinite', '1', '2',
             '3', '4', 'degree']
            sage: S5.degree                                      #optional - gap3
            5
        """
        gap3_session = self._check_valid()
        if not hasattr(self, "_gap_recfields"):
            s = str(gap3_session.eval("RecFields(%s)" % self._name))
            s = s.strip('[] ').replace('\n', '')
            self._gap_recfields = [ss.strip('" ') for ss in s.split(',')]
        return getattr(self, "_gap_recfields")

    def operations(self):
        r"""
        Return a list of the GAP3 operations for the record.

        OUTPUT: list of strings -- operations of the record

        EXAMPLES::

            sage: S5 = gap3.SymmetricGroup(5)              #optional - gap3
            sage: S5.operations()                          #optional - gap3
            [..., 'NormalClosure', 'NormalIntersection', 'Normalizer',
            'NumberConjugacyClasses', 'PCore', 'Radical', 'SylowSubgroup',
            'TrivialSubgroup', 'FusionConjugacyClasses', 'DerivedSeries', ...]
            sage: S5.DerivedSeries()                       #optional - gap3
            [ Group( (1,5), (2,5), (3,5), (4,5) ),
              Subgroup( Group( (1,5), (2,5), (3,5), (4,5) ),
                        [ (1,2,5), (1,3,5), (1,4,5) ] ) ]
        """
        gap3_session = self._check_valid()
        if not hasattr(self, "_gap_operations"):
            s = str(gap3_session.eval("RecFields(%s.operations)" % self._name))
            s = s.strip('[] ').replace('\n', '')
            self._gap_operations = [ss.strip('" ') for ss in s.split(',')]
        return getattr(self, "_gap_operations")

    def __getattr__(self, attrname):
        r"""
        OUTPUT:

        - ``GAP3Record`` -- if ``attrname`` is a field of the GAP record

        - ``ExpectFunction`` -- if ``attrname`` is the name of a GAP3 function

        EXAMPLES::

            sage: # optional - gap3
            sage: S5 = gap3.SymmetricGroup(5)
            sage: S5.__getattr__('Size')
            Size
            sage: gap3.IsFunc(S5.__getattr__('Size'))
            true
            sage: S5.__getattr__('generators')
            [ (1,5), (2,5), (3,5), (4,5) ]
        """
        gap3_session = self._check_valid()
        if attrname[:1] == "_":
            raise AttributeError
        if attrname in self.recfields():
            return gap3_session.new('%s.%s' % (self.name(), attrname))
        return gap3_session._function_element_class()(self, attrname)

    def _tab_completion(self):
        r"""
        Define the list of methods and attributes that will appear for tab
        completion.

        OUTPUT:

        - list of strings -- the available fields and operations of the
          record

        EXAMPLES::

            sage: S5 = gap3.SymmetricGroup(5)              #optional - gap3
            sage: S5._tab_completion()                     #optional - gap3
            [..., 'ConjugacyClassesTry', 'ConjugateSubgroup', 'ConjugateSubgroups',
            'Core', 'DegreeOperation', 'DerivedSeries', 'DerivedSubgroup',
            'Difference', 'DimensionsLoewyFactors', 'DirectProduct', ...]
        """
        names = self.recfields() + self.operations()
        names.sort()
        return names


def gap3_console():
    r"""
    Spawn a new GAP3 command-line session.

    EXAMPLES::

        sage: gap3.console()                             # not tested

                     ########            Lehrstuhl D fuer Mathematik
                   ###    ####           RWTH Aachen
                  ##         ##
                 ##          #             #######            #########
                ##                        #      ##          ## #     ##
                ##           #           #       ##             #      ##
                ####        ##           ##       #             #      ##
                 #####     ###           ##      ##             ##    ##
                   ######### #            #########             #######
                             #                                  #
                            ##           Version 3              #
                           ###           Release 4.4            #
                          ## #           18 Apr 97              #
                         ##  #
                        ##   #  Alice Niemeyer, Werner Nickel,  Martin Schoenert
                       ##    #  Johannes Meier, Alex Wegner,    Thomas Bischops
                      ##     #  Frank Celler,   Juergen Mnich,  Udo Polis
                      ###   ##  Thomas Breuer,  Goetz Pfeiffer, Hans U. Besche
                       ######   Volkmar Felsch, Heiko Theissen, Alexander Hulpke
                                Ansgar Kaup,    Akos Seress,    Erzsebet Horvath
                                Bettina Eick
                                For help enter: ?<return>
        gap>
    """
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%gap3 magics instead.')
    os.system(gap3_cmd)


def gap3_version():
    r"""
    Return the version of GAP3 that you have in your PATH on your computer.

    EXAMPLES::

        sage: gap3_version()                                           # random, optional - gap3
        'lib: v3r4p4 1997/04/18, src: v3r4p0 1994/07/10, sys: usg gcc ansi'
    """
    return gap3.eval('VERSION')[1:-1]
