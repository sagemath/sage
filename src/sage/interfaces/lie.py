r"""
Interface to LiE

LiE is a software package under development at CWI since
January 1988.  Its purpose is to enable mathematicians and
physicists to obtain on-line information as well as to
interactively perform computations of a Lie group theoretic
nature.  It focuses on the representation theory of complex
semisimple (reductive) Lie groups and algebras, and on the
structure of their Weyl groups and root systems.

Type ``lie.[tab]`` for a list of all the functions available
from your LiE install.  Type ``lie.[tab]?`` for LiE's
help about a given function.  Type ``lie(...)`` to create
a new LiE object, and ``lie.eval(...)`` to run a string
using LiE (and get the result back as a string).

To access the LiE interpreter directly, run lie_console().


EXAMPLES::

    sage: a4 = lie('A4')             # optional - lie
    sage: lie.diagram('A4')          # optional - lie
    O---O---O---O
    1   2   3   4
    A4

    sage: lie.diagram(a4)            # optional - lie
    O---O---O---O
    1   2   3   4
    A4

    sage: a4.diagram()               # optional - lie
    O---O---O---O
    1   2   3   4
    A4

    sage: a4.Cartan()                # optional - lie
         [[ 2,-1, 0, 0]
         ,[-1, 2,-1, 0]
         ,[ 0,-1, 2,-1]
         ,[ 0, 0,-1, 2]
         ]
    sage: lie.LR_tensor([3,1],[2,2]) # optional - lie
    1X[5,3]


Tutorial
--------

The following examples are taken from Section 2.1 of the LiE manual.

You can perform basic arithmetic operations in LiE. ::

    sage: lie.eval('19+68') # optional - lie
    '87'
    sage: a = lie('1111111111*1111111111') # optional - lie
    sage: a # optional - lie
    1234567900987654321
    sage: a/1111111111 # optional - lie
    1111111111
    sage: a = lie('345') # optional - lie
    sage: a^2+3*a-5 # optional - lie
    120055
    sage: _ / 7*a # optional - lie
    5916750

Vectors in LiE are created using square brackets.  Notice that
the indexing in LiE is 1-based, unlike Python/Sage which is 0-based. ::

    sage: v = lie('[3,2,6873,-38]') # optional - lie
    sage: v # optional - lie
    [3,2,6873,-38]
    sage: v[3] # optional - lie
    6873
    sage: v+v # optional - lie
    [6,4,13746,-76]
    sage: v*v # optional - lie
    47239586
    sage: v+234786 # optional - lie
    [3,2,6873,-38,234786]
    sage: v-3 # optional - lie
    [3,2,-38]
    sage: v^v # optional - lie
    [3,2,6873,-38,3,2,6873,-38]

You can also work with matrices in LiE. ::

    sage: m = lie('[[1,0,3,3],[12,4,-4,7],[-1,9,8,0],[3,-5,-2,9]]') # optional - lie
    sage: m # optional - lie
         [[ 1, 0, 3,3]
         ,[12, 4,-4,7]
         ,[-1, 9, 8,0]
         ,[ 3,-5,-2,9]
         ]
    sage: print(lie.eval('*'+m._name))  # optional - lie
         [[1,12,-1, 3]
         ,[0, 4, 9,-5]
         ,[3,-4, 8,-2]
         ,[3, 7, 0, 9]
         ]

    sage: m^3 # optional - lie
         [[ 220,   87, 81, 375]
         ,[-168,-1089, 13,1013]
         ,[1550,  357,-55,1593]
         ,[-854, -652, 98,-170]
         ]
    sage: v*m # optional - lie
    [-6960,62055,55061,-319]
    sage: m*v # optional - lie
    [20508,-27714,54999,-14089]
    sage: v*m*v # optional - lie
    378549605
    sage: m+v # optional - lie
         [[ 1, 0,   3,  3]
         ,[12, 4,  -4,  7]
         ,[-1, 9,   8,  0]
         ,[ 3,-5,  -2,  9]
         ,[ 3, 2,6873,-38]
         ]

    sage: m-2 # optional - lie
         [[ 1, 0, 3,3]
         ,[-1, 9, 8,0]
         ,[ 3,-5,-2,9]
         ]


LiE handles multivariate (Laurent) polynomials. ::

    sage: lie('X[1,2]') # optional - lie
    1X[1,2]
    sage: -3*_ # optional - lie
    -3X[1,2]
    sage: _ + lie('4X[-1,4]') # optional - lie
    4X[-1,4] - 3X[ 1,2]
    sage: _^2 # optional - lie
    16X[-2,8] - 24X[ 0,6] +  9X[ 2,4]
    sage: lie('(4X[-1,4]-3X[1,2])*(X[2,0]-X[0,-4])') # optional - lie
    -4X[-1, 0] + 3X[ 1,-2] + 4X[ 1, 4] - 3X[ 3, 2]
    sage: _ - _ # optional - lie
    0X[0,0]


You can call LiE's built-in functions using ``lie.functionname``. ::

    sage: lie.partitions(6) # optional - lie
         [[6,0,0,0,0,0]
         ,[5,1,0,0,0,0]
         ,[4,2,0,0,0,0]
         ,[4,1,1,0,0,0]
         ,[3,3,0,0,0,0]
         ,[3,2,1,0,0,0]
         ,[3,1,1,1,0,0]
         ,[2,2,2,0,0,0]
         ,[2,2,1,1,0,0]
         ,[2,1,1,1,1,0]
         ,[1,1,1,1,1,1]
         ]
    sage: lie.diagram('E8') # optional - lie
            O 2
            |
            |
    O---O---O---O---O---O---O
    1   3   4   5   6   7   8
    E8


You can define your own functions in LiE using lie.eval .  Once you've defined
a function (say f), you can call it using lie.f ; however, user-defined functions
do not show up when using tab-completion. ::

    sage: lie.eval('f(int x) = 2*x') # optional - lie
    ''
    sage: lie.f(984) # optional - lie
    1968
    sage: lie.eval('f(int n) = a=3*n-7; if a < 0 then a = -a fi; 7^a+a^3-4*a-57') # optional - lie
    ''
    sage: lie.f(2) # optional - lie
    -53
    sage: lie.f(5) # optional - lie
    5765224



LiE's help can be accessed through lie.help('functionname') where
functionname is the function you want to receive help for. ::

   sage: print(lie.help('diagram'))  # optional - lie
   diagram(g).   Prints the Dynkin diagram of g, also indicating
      the type of each simple component printed, and labeling the nodes as
      done by Bourbaki (for the second and further simple components the
      labels are given an offset so as to make them disjoint from earlier
      labels). The labeling of the vertices of the Dynkin diagram prescribes
      the order of the coordinates of root- and weight vectors used in LiE.

This can also be accessed with lie.functionname? .



With the exception of groups, all LiE data types can be converted into
native Sage data types by calling the .sage() method.

Integers::

    sage: a = lie('1234') # optional - lie
    sage: b = a.sage(); b # optional - lie
    1234
    sage: type(b) # optional - lie
    <class 'sage.rings.integer.Integer'>

Vectors::

    sage: a = lie('[1,2,3]') # optional - lie
    sage: b = a.sage(); b # optional - lie
    [1, 2, 3]
    sage: type(b) # optional - lie
    <... 'list'>

Matrices::

    sage: a = lie('[[1,2],[3,4]]') # optional - lie
    sage: b = a.sage(); b # optional - lie
    [1 2]
    [3 4]
    sage: type(b) # optional - lie
    <class 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>

Polynomials::

    sage: a = lie('X[1,2] - 2*X[2,1]') # optional - lie
    sage: b = a.sage(); b              # optional - lie
    -2*x0^2*x1 + x0*x1^2
    sage: type(b)                      # optional - lie
    <class 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>

Text::

    sage: a = lie('"text"') # optional - lie
    sage: b = a.sage(); b # optional - lie
    'text'
    sage: type(b) # optional - lie
    <... 'str'>


LiE can be programmed using the Sage interface as well. Section 5.1.5
of the manual gives an example of a function written in LiE's language
which evaluates a polynomial at a point.  Below is a (roughly) direct
translation of that program into Python / Sage. ::

    sage: def eval_pol(p, pt): # optional - lie
    ....:     s = 0
    ....:     for i in range(1,p.length().sage()+1):
    ....:         m = 1
    ....:         for j in range(1,pt.size().sage()+1):
    ....:             m *= pt[j]^p.expon(i)[j]
    ....:         s += p.coef(i)*m
    ....:     return s
    sage: a = lie('X[1,2]') # optional - lie
    sage: b1 = lie('[1,2]') # optional - lie
    sage: b2 = lie('[2,3]') # optional - lie
    sage: eval_pol(a, b1) # optional - lie
    4
    sage: eval_pol(a, b2) # optional - lie
    18



AUTHORS:

- Mike Hansen 2007-08-27
- William Stein (template)
"""

##########################################################################
#
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
#
##########################################################################

from .expect import Expect, ExpectElement, ExpectFunction, FunctionElement
from sage.interfaces.interface import AsciiArtString
from sage.misc.misc_c import prod
from sage.env import DOT_SAGE, LIE_INFO_DIR
from sage.misc.sage_eval import sage_eval
from sage.interfaces.tab_completion import ExtraTabCompletion
from sage.docs.instancedoc import instancedoc
import os


COMMANDS_CACHE = '%s/lie_commandlist_cache.sobj' % DOT_SAGE
HELP_CACHE = '%s/lie_helpdict_cache.sobj' % DOT_SAGE


class LiE(ExtraTabCompletion, Expect):
    r"""
    Interface to the LiE interpreter.

    Type ``lie.[tab]`` for a list of all the functions available
    from your LiE install.  Type ``lie.[tab]?`` for LiE's
    help about a given function.  Type ``lie(...)`` to create
    a new LiE object, and ``lie.eval(...)`` to run a string
    using LiE (and get the result back as a string).
    """

    def __init__(self,
                 maxread=None, script_subdirectory=None,
                 logfile=None,
                 server=None):
        """
        EXAMPLES::

            sage: lie == loads(dumps(lie))
            True
        """
        Expect.__init__(self,

                        # The capitalized version of this is used for printing.
                        name='LiE',

                        # This is regexp of the input prompt.  If you can change
                        # it to be very obfuscated that would be better.   Even
                        # better is to use sequence numbers.
                        prompt='> ',

                        # This is the command that starts up your program
                        command="bash lie",

                        server=server,
                        script_subdirectory=script_subdirectory,

                        # If this is true, then whenever the user presses Control-C to
                        # interrupt a calculation, the whole interface is restarted.
                        restart_on_ctrlc=False,

                        # If true, print out a message when starting
                        # up the command when you first send a command
                        # to this interface.
                        verbose_start=False,

                        logfile=logfile,

                        # If an input is longer than this number of characters, then
                        # try to switch to outputting to a file.
                        eval_using_file_cutoff=1024)

        self._seq = 0

        self._tab_completion_dict = None
        self._tab_completion_list = None
        self._help_dict = None

    def _read_info_files(self, use_disk_cache=True):
        """
        EXAMPLES::

            sage: from sage.interfaces.lie import LiE
            sage: lie = LiE()
            sage: lie._tab_completion_list is None
            True
            sage: lie._read_info_files(use_disk_cache=False) #optional - lie
            sage: lie._tab_completion_list # optional - lie
            ['Adams',
             ...
             'history',
             ...
             'sort',
             ...
             'version',
             'void',
             'write']
        """
        import sage.misc.persist
        if use_disk_cache:
            try:
                trait_dict = sage.misc.persist.load(COMMANDS_CACHE)
                help_dict = sage.misc.persist.load(HELP_CACHE)
                v = []
                for key in trait_dict:
                    v += trait_dict[key]
                self._tab_completion_list = sorted(v)
                self._tab_completion_dict = trait_dict
                self._help_dict = help_dict
                return
            except IOError:
                pass

        # Go through INFO.3 and get the necessary information
        filenames = ['INFO.3', 'INFO.0']
        commands = {}
        commands['vid'] = []
        help = {}

        for f in filenames:
            info = open(os.path.join(LIE_INFO_DIR, f))
            prev_command = ""
            help_text = ""
            for line in info:
                # If the line doesn't start with an "@", then
                # it is part of the help text for the previous
                # command
                if len(line) == 0 or line[0] != "@":
                    if prev_command != "":
                        help_text += line
                    continue

                # Do not add not completions that do not start with an
                # alphabetical character or that contain 'silence'
                if len(line) > 1 and (not line[1].isalpha() or line.find('silence') != -1):
                    help[prev_command] = help.get(prev_command, "") + help_text
                    help_text = ""
                    prev_command = ""
                    continue

                # At this point we should be at the start of a new
                # command definition

                # Get the type of the first argument of the command
                i = line.find('(')
                if line[i + 1] == ")":
                    t = 'vid'
                else:
                    t = line[i + 1:i + 4]

                # Save the help text for the command
                help[prev_command] = help.get(prev_command, "") + help_text
                help_text = ""
                prev_command = line[1:i]

                # Add the command
                if t in commands:
                    commands[t].append(line[1:i])
                else:
                    commands[t] = [line[1:i]]

            # Take care of the last help text which doesn't get processed
            # since there's no following @ symbol
            help[prev_command] = help.get(prev_command, "") + help_text

        info.close()

        # Build the list of all possible command completions
        l = []
        for key in commands:
            l += commands[key]

        # Save the data
        self._tab_completion_dict = commands
        self._tab_completion_list = sorted(l)
        self._help_dict = help

        # Write them to file
        if use_disk_cache:
            sage.misc.persist.save(commands, COMMANDS_CACHE)
            sage.misc.persist.save(help, HELP_CACHE)

    def _repr_(self) -> str:
        """
        EXAMPLES::

            sage: lie
            LiE Interpreter
        """
        return 'LiE Interpreter'

    def __reduce__(self):
        """
        EXAMPLES::

            sage: lie.__reduce__()
            (<function reduce_load_lie at 0x...>, ())
        """
        return reduce_load_lie, tuple([])

    def _function_class(self):
        """
        EXAMPLES::

            sage: lie._function_class()
            <class 'sage.interfaces.lie.LiEFunction'>
        """
        return LiEFunction

    def _quit_string(self) -> str:
        """
        EXAMPLES::

            sage: lie._quit_string()
            'quit'
        """
        return 'quit'

    def _read_in_file_command(self, filename):
        """
        EXAMPLES::

            sage: lie._read_in_file_command('testfile')
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _tab_completion(self, type=None, verbose=False, use_disk_cache=True):
        """
        EXAMPLES::

            sage: lie._tab_completion() # optional - lie
            ['Adams',
             ...
             'Cartan_type',
             ...
             'cent_roots',
             ...
             'n_comp',
             ...
             'write']
        """
        if self._tab_completion_dict is None:
            self._read_info_files()
        if type:
            return sorted(self._tab_completion_dict[type])
        else:
            return self._tab_completion_list

    def _an_element_impl(self):
        """
        EXAMPLES::

            sage: lie._an_element_impl() # optional - lie
            0
        """
        return self(0)

    def read(self, filename):
        r"""
        EXAMPLES::

            sage: filename = tmp_filename()
            sage: with open(filename, 'w') as f:
            ....:     _ = f.write('x = 2\n')
            sage: lie.read(filename)  # optional - lie
            sage: lie.get('x')        # optional - lie
            '2'
            sage: import os
            sage: os.unlink(filename)
        """
        self.eval('read %s' % filename)

    def console(self):
        """
        Spawn a new LiE command-line session.

        EXAMPLES::

            sage: lie.console()                    # not tested
            LiE version 2.2.2 created on Sep 26 2007 at 18:13:19
            Authors: Arjeh M. Cohen, Marc van Leeuwen, Bert Lisser.
            Free source code distribution
            ...
        """
        lie_console()

    def version(self) -> str:
        """
        EXAMPLES::

            sage: lie.version() # optional - lie
            '2.2'
        """
        return lie_version()

    def _object_class(self):
        """
        EXAMPLES::

            sage: lie._object_class()
            <class 'sage.interfaces.lie.LiEElement'>
        """
        return LiEElement

    def _true_symbol(self) -> str:
        """
        EXAMPLES::

            sage: lie._true_symbol()
            '1'
        """
        return '1'

    def _false_symbol(self) -> str:
        """
        EXAMPLES::

            sage: lie._false_symbol()
            '0'
        """
        return '0'

    def _equality_symbol(self) -> str:
        """
        EXAMPLES::

            sage: lie._equality_symbol()
            '=='
        """
        return '=='

    def help(self, command):
        """
        Return a string of the LiE help for command.

        EXAMPLES::

            sage: lie.help('diagram') # optional - lie
            'diagram(g)...'
        """
        # return help on a given command.
        if self._help_dict is None:
            self._read_info_files()
        try:
            return self._help_dict[command]
        except KeyError:
            return "Could not find help for " + command

    def _eval_line(self, line, allow_use_file=True, wait_for_prompt=True, restart_if_needed=False):
        """
        EXAMPLES::

            sage: lie._eval_line('2+2') # optional - lie
            '     4'
            sage: lie._eval_line('diagram(2)') # optional - lie
            Traceback (most recent call last):
            ...
            RuntimeError: An error occurred running a LiE command:
            Argument types do not match in call. Types are: diagram(bin).
            Valid argument types are for instance: diagram(grp).
        """
        out = Expect._eval_line(self, line, allow_use_file=allow_use_file, wait_for_prompt=wait_for_prompt)
        # Check to see if an error has occurred
        err = max(out.find("\n(in"), out.find('not defined'), out.find('Argument types'))
        if err != -1:
            raise RuntimeError("An error occurred running a LiE command:\n%s" % (out.replace('\r\n', '\n')))
        return out

    def eval(self, code, strip=True, **kwds):
        """
        EXAMPLES::

            sage: lie.eval('2+2')  # optional - lie
            '4'
        """
        s = Expect.eval(self, code, strip=True, **kwds)
        # return s.strip()
        if len(s) > 0 and s.find("\n") != -1:
            return s
        else:
            return s.strip()

    def set(self, var, value):
        """
        Set the variable var to the given value.

        EXAMPLES::

            sage: lie.set('x', '2')  # optional - lie
            sage: lie.get('x')       # optional - lie
            '2'
        """
        cmd = '%s=%s' % (var, value)
        out = self.eval(cmd)
        i = min(out.find('not defined'), out.find(r'\(in'), out.find('Argument types'))
        if i != -1:
            raise RuntimeError(out)

    def get(self, var):
        """
        Get the value of the variable var.

        EXAMPLES::

            sage: lie.set('x', '2')  # optional - lie
            sage: lie.get('x')       # optional - lie
            '2'
        """
        return self.eval('%s' % var)

    def get_using_file(self, var):
        """
        EXAMPLES::

            sage: lie.get_using_file('x')
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def function_call(self, function, args=None, kwds=None):
        """
        EXAMPLES::

            sage: lie.function_call("diagram", args=['A4']) # optional - lie
            O---O---O---O
            1   2   3   4
            A4
        """
        # If function just prints something on the screen rather than
        # returning an object, then we return an AsciiArtString rather
        # than a LiEElement
        if function in ['diagram', 'setdefault', 'print_tab', 'type', 'factor', 'void', 'gcol']:
            args, kwds = self._convert_args_kwds(args, kwds)
            cmd = "%s(%s)" % (function, ",".join(s.name() for s in args))
            return AsciiArtString(self.eval(cmd))

        return Expect.function_call(self, function, args, kwds)

    def _function_element_class(self):
        """
        EXAMPLES::

            sage: lie._function_element_class()
            <class 'sage.interfaces.lie.LiEFunctionElement'>
        """
        return LiEFunctionElement


@instancedoc
class LiEElement(ExtraTabCompletion, ExpectElement):
    def _tab_completion(self):
        """
        Return the possible tab completions for ``self``.

        EXAMPLES::

            sage: a4 = lie('A4')   # optional - lie
            sage: a4._tab_completion() # optional - lie
            ['Cartan',
             ...
             'center',
             'det_Cartan',
             'diagram',
             ...
             'n_comp',
             ...
             'res_mat']
        """
        return self.parent()._tab_completion(type=self.type())

    def type(self):
        """
        EXAMPLES::

            sage: m = lie('[[1,0,3,3],[12,4,-4,7],[-1,9,8,0],[3,-5,-2,9]]') # optional - lie
            sage: m.type() # optional - lie
            'mat'
        """
        t = self.parent().eval('type(%s)' % self._name)
        i = t.find(':')
        return t[i + 1:].strip()

    def _matrix_(self, R=None):
        """
        EXAMPLES::

            sage: m = lie('[[1,0,3,3],[12,4,-4,7],[-1,9,8,0],[3,-5,-2,9]]') # optional - lie
            sage: matrix(m)  # optional - lie
            [ 1  0  3  3]
            [12  4 -4  7]
            [-1  9  8  0]
            [ 3 -5 -2  9]
            sage: matrix(RDF, m) # optional - lie
            [ 1.0  0.0  3.0  3.0]
            [12.0  4.0 -4.0  7.0]
            [-1.0  9.0  8.0  0.0]
            [ 3.0 -5.0 -2.0  9.0]
        """
        self._check_valid()
        if self.type() == 'mat':
            m = self.sage()
            if R is not None:
                m = m.change_ring(R)
            return m
        else:
            raise ValueError("not a matrix")

    def _sage_(self):
        """
        EXAMPLES::

            sage: m = lie('[[1,0,3,3],[12,4,-4,7],[-1,9,8,0],[3,-5,-2,9]]') # optional - lie
            sage: m.sage()  # optional - lie
            [ 1  0  3  3]
            [12  4 -4  7]
            [-1  9  8  0]
            [ 3 -5 -2  9]

        """
        t = self.type()
        if t == 'grp':
            raise ValueError("cannot convert Lie groups to native Sage objects")
        elif t == 'mat':
            import sage.matrix.constructor
            data = sage_eval(str(self).replace('\n', '').strip())
            return sage.matrix.constructor.matrix(data)
        elif t == 'pol':
            from sage.rings.all import PolynomialRing, QQ

            # Figure out the number of variables
            s = str(self)
            open_bracket = s.find('[')
            close_bracket = s.find(']')
            nvars = len(s[open_bracket:close_bracket].split(','))

            # create the polynomial ring
            R = PolynomialRing(QQ, nvars, 'x')
            x = R.gens()
            pol = R.zero()

            # Split up the polynomials into terms
            terms = []
            for termgrp in s.split(' - '):
                # The first entry in termgrp has
                # a negative coefficient
                termgrp = "-" + termgrp.strip()
                terms += termgrp.split('+')
            # Make sure we don't accidentally add a negative
            # sign to the first monomial
            if s[0] != "-":
                terms[0] = terms[0][1:]

            # go through all the terms in s
            for term in terms:
                xpos = term.find('X')
                coef = eval(term[:xpos].strip())
                exps = eval(term[xpos + 1:].strip())
                monomial = prod([x[i]**exps[i] for i in range(nvars)])
                pol += coef * monomial

            return pol
        elif t == 'tex':
            return repr(self)
        elif t == 'vid':
            return None
        else:
            return ExpectElement._sage_(self)


@instancedoc
class LiEFunctionElement(FunctionElement):
    def _instancedoc_(self):
        """
        EXAMPLES::

            sage: a4 = lie('A4')  # optional - lie
            sage: a4.diagram.__doc__  # optional - lie
            'diagram(g)...'
        """
        M = self._obj.parent()
        return M.help(self._name)


@instancedoc
class LiEFunction(ExpectFunction):
    def _instancedoc_(self):
        """
        Return the help for ``self``.

        EXAMPLES::

            sage: lie.diagram.__doc__  # optional - lie
            'diagram(g)...'
        """
        M = self._parent
        return M.help(self._name)


def is_LiEElement(x) -> bool:
    """
    EXAMPLES::

        sage: from sage.interfaces.lie import is_LiEElement
        sage: l = lie(2) # optional - lie
        sage: is_LiEElement(l) # optional - lie
        True
        sage: is_LiEElement(2)
        False
    """
    return isinstance(x, LiEElement)


# An instance
lie = LiE()


def reduce_load_lie():
    """
    EXAMPLES::

        sage: from sage.interfaces.lie import reduce_load_lie
        sage: reduce_load_lie()
        LiE Interpreter
    """
    return lie


def lie_console():
    """
    Spawn a new LiE command-line session.

    EXAMPLES::

        sage: from sage.interfaces.lie import lie_console
        sage: lie_console()                    # not tested
        LiE version 2.2.2 created on Sep 26 2007 at 18:13:19
        Authors: Arjeh M. Cohen, Marc van Leeuwen, Bert Lisser.
        Free source code distribution
        ...
    """
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%lie magics instead.')
    os.system('bash `which lie`')


def lie_version():
    """
    EXAMPLES::

        sage: from sage.interfaces.lie import lie_version
        sage: lie_version() # optional - lie
        '2.2'
    """
    with open(os.path.join(LIE_INFO_DIR, 'INFO.0')) as f:
        lines = f.readlines()
    i = lines.index('@version()\n')
    return lines[i + 1].split()[1]
