r"""
Interface to GNU Octave

GNU Octave is a free software (GPL) MATLAB-like program with numerical
routines for integrating, solving systems of equations, special
functions, and solving (numerically) differential equations. Please see
http://octave.org/ for more details.

The commands in this section only work if you have the optional
"octave" interpreter installed and available in your PATH. It's not
necessary to install any special Sage packages.

EXAMPLES::

    sage: octave.eval('2+2')    # optional - octave
    'ans = 4'

    sage: a = octave(10)        # optional - octave
    sage: a**10                 # optional - octave
    1e+10

LOG: - creation (William Stein) - ? (David Joyner, 2005-12-18) -
Examples (David Joyner, 2005-01-03)

Computation of Special Functions
--------------------------------

Octave implements computation of the following special functions
(see the maxima and gp interfaces for even more special
functions)::

    airy
        Airy functions of the first and second kind, and their derivatives.
        airy(0,x) = Ai(x), airy(1,x) = Ai'(x), airy(2,x) = Bi(x), airy(3,x) = Bi'(x)
    besselj
        Bessel functions of the first kind.
    bessely
        Bessel functions of the second kind.
    besseli
        Modified Bessel functions of the first kind.
    besselk
        Modified Bessel functions of the second kind.
    besselh
        Compute Hankel functions of the first (k = 1) or second (k = 2) kind.
    beta
        The Beta function,
              beta (a, b) = gamma (a) * gamma (b) / gamma (a + b).
    betainc
        The incomplete Beta function,
    erf
        The error function,
    erfinv
        The inverse of the error function.
    gamma
        The Gamma function,
    gammainc
        The incomplete gamma function,

For example,

::

    sage: # optional - octave
    sage: octave("airy(3,2)")
    4.10068
    sage: octave("beta(2,2)")
    0.166667
    sage: octave("betainc(0.2,2,2)")
    0.104
    sage: octave("besselh(0,2)")
    (0.223891,0.510376)
    sage: octave("besselh(0,1)")
    (0.765198,0.088257)
    sage: octave("besseli(1,2)")
    1.59064
    sage: octave("besselj(1,2)")
    0.576725
    sage: octave("besselk(1,2)")
    0.139866
    sage: octave("erf(0)")
    0
    sage: octave("erf(1)")
    0.842701
    sage: octave("erfinv(0.842)")
    0.998315
    sage: octave("gamma(1.5)")
    0.886227
    sage: octave("gammainc(1.5,1)")
    0.77687

Tutorial
--------

EXAMPLES::

    sage: # optional - octave
    sage: octave('4+10')
    14
    sage: octave('date')
    18-Oct-2007
    sage: octave('5*10 + 6')
    56
    sage: octave('(6+6)/3')
    4
    sage: octave('9')^2
    81
    sage: a = octave(10); b = octave(20); c = octave(30)
    sage: avg = (a+b+c)/3
    sage: avg
    20
    sage: parent(avg)
    Octave

::

    sage: # optional - octave
    sage: my_scalar = octave('3.1415')
    sage: my_scalar
    3.1415
    sage: my_vector1 = octave('[1,5,7]')
    sage: my_vector1
    1     5     7
    sage: my_vector2 = octave('[1;5;7]')
    sage: my_vector2
    1
    5
    7
    sage: my_vector1 * my_vector2
    75
"""

# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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
from .expect import Expect, ExpectElement
import pexpect
from sage.misc.verbose import verbose
from sage.misc.instancedoc import instancedoc
from sage.misc.temporary_file import tmp_filename
from sage.cpython.string import bytes_to_str


class Octave(Expect):
    r"""
    Interface to the Octave interpreter.

    EXAMPLES::

        sage: octave.eval("a = [ 1, 1, 2; 3, 5, 8; 13, 21, 33 ]").strip()    # optional - octave
        'a =\n\n 1 1 2\n 3 5 8\n 13 21 33'
        sage: octave.eval("b = [ 1; 3; 13]").strip()                         # optional - octave
        'b =\n\n 1\n 3\n 13'

    The following solves the linear equation: a*c = b::

        sage: octave.eval(r"c=a \ b").strip()          # optional - octave  # abs tol 0.01
        'c =\n\n 1\n -0\n 0'
        sage: octave.eval("c").strip()                 # optional - octave  # abs tol 0.01
        'c =\n\n 1\n -0\n 0'

    TESTS:

    We check that the interface can handle large inputs (see :issue:`940`)::

        sage: t = '"{}"'.format(10^10000)
        sage: a = octave(t)                     # optional - octave
        sage: str(a) == ' {}'.format(10^10000)  # optional - octave
        True
    """

    def __init__(self, maxread=None, script_subdirectory=None, logfile=None,
            server=None, server_tmpdir=None, seed=None, command=None):
        """
        EXAMPLES::

            sage: octave == loads(dumps(octave))
            True
        """
        if command is None:
            command = os.getenv('SAGE_OCTAVE_COMMAND') or 'octave-cli'
        if server is None:
            server = os.getenv('SAGE_OCTAVE_SERVER') or None
        # Use a temporary workspace file.
        workspace_file = tmp_filename()
        Expect.__init__(self,
                        name='octave',
                        # We want the prompt sequence to be unique to avoid confusion with syntax error messages containing >>>
                        prompt=r'octave\:\d+> ',
                        # We don't want any pagination of output
                        command=command + f" --no-line-editing --silent --eval 'PS2(PS1());more off; octave_core_file_name (\"{workspace_file}\")' --persist",
                        maxread=maxread,
                        server=server,
                        server_tmpdir=server_tmpdir,
                        script_subdirectory=script_subdirectory,
                        restart_on_ctrlc=False,
                        verbose_start=False,
                        logfile=logfile,
                        eval_using_file_cutoff=100)
        self._seed = seed

    def set_seed(self, seed=None):
        """
        Set the seed for the random number generator
        for this octave interpreter.

        EXAMPLES::

            sage: o = Octave() # optional - octave
            sage: o.set_seed(1) # optional - octave
            1
            sage: [o.rand() for i in range(5)] # optional - octave
            [ 0.134364,  0.847434,  0.763775,  0.255069,  0.495435]
        """
        if seed is None:
            seed = self.rand_seed()
        self.eval("rand('state',%d)" % seed)
        self._seed = seed
        return seed

    def __reduce__(self):
        """
        EXAMPLES::

            sage: Octave().__reduce__()
            (<function reduce_load_Octave at 0x...>, ())
        """
        return reduce_load_Octave, tuple([])

    def _read_in_file_command(self, filename):
        """
        EXAMPLES::

            sage: filename = tmp_filename()
            sage: octave._read_in_file_command(filename)
            'source("...");'
        """
        return 'source("%s");' % filename

    def _quit_string(self):
        """
        EXAMPLES::

            sage: octave._quit_string()
            'quit;'
        """
        return 'quit;'

    def _install_hints(self):
        """
        Return hints on how to install Octave.

        EXAMPLES::

            sage: print(octave._install_hints())
            You must get ...
        """
        return """
        You must get the program "octave" in order to use Octave
        from Sage.   You can read all about Octave at
                https://www.gnu.org/software/octave/

        LINUX:
           Do apt-get install octave as root on your machine (Ubuntu/Debian).
           Other Linux systems have octave too.

        OS X:
           * This website has details on OS X builds of Octave:
                    http://wiki.octave.org/Octave_for_MacOS_X
           * Darwin ports and fink have Octave as well.
        """

    def _eval_line(self, line, reformat=True, allow_use_file=False,
                   wait_for_prompt=True, restart_if_needed=False):
        """
        EXAMPLES::

            sage: print(octave._eval_line('2+2'))  #optional - octave
              ans =  4
        """
        from pexpect.exceptions import EOF
        if not wait_for_prompt:
            return Expect._eval_line(self, line)
        if line == '':
            return ''
        if self._expect is None:
            self._start()
        if allow_use_file and len(line) > 3000:
            return self._eval_line_using_file(line)
        try:
            E = self._expect
            # debug
            # self._synchronize(cmd='1+%s\n')
            verbose("in = '%s'" % line, level=3)
            E.sendline(line)
            E.expect(self._prompt)
            out = bytes_to_str(E.before)
            # debug
            verbose("out = '%s'" % out, level=3)
        except EOF:
            if self._quit_string() in line:
                return ''
        except KeyboardInterrupt:
            self._keyboard_interrupt()
        try:
            if reformat:
                if 'syntax error' in out:
                    raise SyntaxError(out)
            out = "\n".join(out.splitlines()[1:])
            return out
        except NameError:
            return ''

    def _keyboard_interrupt(self):
        print("CntrlC: Interrupting %s..." % self)
        if self._restart_on_ctrlc:
            try:
                self._expect.close(force=1)
            except pexpect.ExceptionPexpect as msg:
                raise RuntimeError("THIS IS A BUG -- PLEASE REPORT. This should never happen.\n" + msg)
            self._start()
            raise KeyboardInterrupt("Restarting %s (WARNING: all variables defined in previous session are now invalid)" % self)
        else:
            self._expect.send('\003')  # control-c
            raise KeyboardInterrupt("Ctrl-c pressed while running %s" % self)

    def quit(self, verbose=False):
        """
        EXAMPLES::

            sage: o = Octave()
            sage: o._start()    # optional - octave
            sage: o.quit(True)  # optional - octave
            Exiting spawned Octave process.
        """
        # Don't bother, since it just hangs in some cases, and it
        # isn't necessary, since octave behaves well with respect
        # to signals.
        if self._expect is not None:
            if verbose:
                print("Exiting spawned %s process." % self)
        return

    def _start(self):
        """
        Start the Octave process.

        EXAMPLES::

            sage: # optional - octave
            sage: o = Octave()
            sage: o.is_running()
            False
            sage: o._start()
            sage: o.is_running()
            True
        """
        Expect._start(self)
        self.eval("page_screen_output=0;")
        self.eval("format none;")
        # set random seed
        self.set_seed(self._seed)

    def _equality_symbol(self):
        """
        EXAMPLES::

            sage: octave('0 == 1')  # optional - octave
             0
            sage: octave('1 == 1')  # optional - octave
             1
        """
        return '=='

    def _true_symbol(self):
        """
        EXAMPLES::

            sage: octave('1 == 1')  # optional - octave
             1
        """
        return '1'

    def _false_symbol(self):
        """
        EXAMPLES::

            sage: octave('0 == 1')  # optional - octave
             0
        """
        return '0'

    def set(self, var, value):
        """
        Set the variable ``var`` to the given ``value``.

        EXAMPLES::

            sage: octave.set('x', '2') # optional - octave
            sage: octave.get('x') # optional - octave
            ' 2'
        """
        cmd = '%s=%s;' % (var, value)
        out = self.eval(cmd)
        if out.find("error") != -1 or out.find("Error") != -1:
            raise TypeError("Error executing code in Octave\nCODE:\n\t%s\nOctave ERROR:\n\t%s" % (cmd, out))

    def get(self, var):
        """
        Get the value of the variable ``var``.

        EXAMPLES::

            sage: octave.set('x', '2') # optional - octave
            sage: octave.get('x') # optional - octave
            ' 2'
        """
        s = self.eval('%s' % var)
        i = s.find('=')
        return s[i+1:]

    def clear(self, var):
        """
        Clear the variable named var.

        EXAMPLES::

            sage: octave.set('x', '2') # optional - octave
            sage: octave.clear('x')    # optional - octave
            sage: octave.get('x')      # optional - octave
            "error: 'x' undefined near line ... column 1"
        """
        self.eval('clear %s' % var)

    def console(self):
        """
        Spawn a new Octave command-line session.

        This requires that the optional octave program be installed and in
        your PATH, but no optional Sage packages need be installed.

        EXAMPLES::

            sage: octave_console()         # not tested
            GNU Octave, version 2.1.73 (i386-apple-darwin8.5.3).
            Copyright (C) 2006 John W. Eaton.
            ...
            octave:1> 2+3
            ans = 5
            octave:2> [ctl-d]

        Pressing ctrl-d exits the octave console and returns you to Sage.
        octave, like Sage, remembers its history from one session to
        another.
        """
        octave_console()

    def version(self):
        r"""
        Return the version of Octave.

        OUTPUT: string

        EXAMPLES::

            sage: v = octave.version()   # optional - octave
            sage: v                      # optional - octave; random
            '2.13.7'

            sage: import re
            sage: assert re.match(r"\d+\.\d+\.\d+", v)  is not None # optional - octave
        """
        return str(self("version")).strip()

    def solve_linear_system(self, A, b):
        r"""
        Use Octave to compute a solution x to ``A*x = b``, as a list.

        INPUT:

        - ``A`` -- mxn matrix A with entries in `\QQ` or `\RR`

        - ``b`` -- m-vector b entries in `\QQ` or `\RR` (resp)

        OUTPUT: list x (if it exists) which solves ``M*x = b``

        EXAMPLES::

            sage: M33 = MatrixSpace(QQ,3,3)
            sage: A   = M33([1,2,3,4,5,6,7,8,0])
            sage: V3  = VectorSpace(QQ,3)
            sage: b   = V3([1,2,3])
            sage: octave.solve_linear_system(A,b)    # optional - octave (and output is slightly random in low order bits)
            [-0.33333299999999999, 0.66666700000000001, -3.5236600000000002e-18]

        AUTHORS:

        - David Joyner and William Stein
        """
        m = A.nrows()
        if m != len(b):
            raise ValueError("dimensions of A and b must be compatible")
        from sage.matrix.matrix_space import MatrixSpace
        from sage.rings.rational_field import QQ
        MS = MatrixSpace(QQ, m, 1)
        b = MS(list(b))  # converted b to a "column vector"
        sA = self.sage2octave_matrix_string(A)
        sb = self.sage2octave_matrix_string(b)
        self.eval("a = " + sA)
        self.eval("b = " + sb)
        soln = octave.eval(r"c = a \ b")
        soln = soln.replace("\n\n ", "[")
        soln = soln.rstrip() + "]"
        soln = soln.replace("\n", ",")
        sol = soln[3:]
        return eval(sol)

    def sage2octave_matrix_string(self, A):
        """
        Return an octave matrix from a Sage matrix.

        INPUT:

        - ``A`` - Sage matrix with entries in the rationals or reals

        OUTPUT: string that evaluates to an Octave matrix

        EXAMPLES::

            sage: M33 = MatrixSpace(QQ,3,3)
            sage: A = M33([1,2,3,4,5,6,7,8,0])
            sage: octave.sage2octave_matrix_string(A)   # optional - octave
            '[1, 2, 3; 4, 5, 6; 7, 8, 0]'

        AUTHORS:

        - David Joyner and William Stein
        """
        return str(A.rows()).replace('), (', '; ').replace('(', '').replace(')', '')

    def de_system_plot(self, f, ics, trange):
        r"""
        Plot (using octave's interface to gnuplot) the solution to a
        `2\times 2` system of differential equations.

        INPUT:

        - ``f`` -- a pair of strings representing the
          differential equations; the independent variable must be called x
          and the dependent variable must be called y

        - ``ics`` -- a pair [x0,y0] such that x(t0) = x0, y(t0) = y0

        - ``trange`` -- a pair [t0,t1]

        OUTPUT: a gnuplot window appears

        EXAMPLES::

            sage: octave.de_system_plot(['x+y','x-y'], [1,-1], [0,2])  # not tested -- does this actually work (on OS X it fails for me -- William Stein, 2007-10)

        This should yield the two plots `(t,x(t)), (t,y(t))` on the
        same graph (the `t`-axis is the horizontal axis) of the
        system of ODEs

        .. MATH::

                       x' = x+y, x(0) = 1;\qquad y' = x-y, y(0) = -1,                     \quad\text{for}\quad 0 < t < 2.
        """
        eqn1 = f[0].replace('x', 'x(1)').replace('y', 'x(2)')
        eqn2 = f[1].replace('x', 'x(1)').replace('y', 'x(2)')
        fcn = "function xdot = f(x,t) xdot(1) = %s; xdot(2) = %s; endfunction" % (eqn1, eqn2)
        self.eval(fcn)
        x0_eqn = "x0 = [%s; %s]" % (ics[0], ics[1])
        self.eval(x0_eqn)
        t_eqn = "t = linspace(%s, %s, 200)'" % (trange[0], trange[1])
        self.eval(t_eqn)
        x_eqn = 'x = lsode("f",x0,t);'
        self.eval(x_eqn)
        self.eval("plot(t,x)")

    def _object_class(self):
        """
        EXAMPLES::

            sage: octave._object_class()
            <class 'sage.interfaces.octave.OctaveElement'>
        """
        return OctaveElement


octave_functions = set()


def to_complex(octave_string, R):
    r"""
    Helper function to convert octave complex number.

    TESTS::

        sage: from sage.interfaces.octave import to_complex
        sage: to_complex('(0,1)', CDF)
        1.0*I
        sage: to_complex('(1.3231,-0.2)', CDF)
        1.3231 - 0.2*I
    """
    real, imag = octave_string.strip('() ').split(',')
    return R(float(real), float(imag))


@instancedoc
class OctaveElement(ExpectElement):
    def _get_sage_ring(self):
        r"""
        TESTS::

            sage: octave('1')._get_sage_ring()  # optional - octave
            Real Double Field
            sage: octave('I')._get_sage_ring()  # optional - octave
            Complex Double Field
            sage: octave('[]')._get_sage_ring() # optional - octave
            Real Double Field
        """
        if self.isinteger():
            import sage.rings.integer_ring
            return sage.rings.integer_ring.ZZ
        elif self.isreal():
            import sage.rings.real_double
            return sage.rings.real_double.RDF
        elif self.iscomplex():
            import sage.rings.complex_double
            return sage.rings.complex_double.CDF
        else:
            raise TypeError("no Sage ring associated to this element.")

    def __bool__(self):
        r"""
        Test whether this element is nonzero.

        EXAMPLES::

            sage: # optional - octave
            sage: bool(octave('0'))
            False
            sage: bool(octave('[]'))
            False
            sage: bool(octave('[0,0]'))
            False
            sage: bool(octave('[0,0,0;0,0,0]'))
            False

            sage: bool(octave('0.1'))               # optional - octave
            True
            sage: bool(octave('[0,1,0]'))           # optional - octave
            True
            sage: bool(octave('[0,0,-0.1;0,0,0]'))  # optional - octave
            True
        """
        return str(self) != ' [](0x0)' and any(x != '0' for x in str(self).split())

    def _matrix_(self, R=None):
        r"""
        Return Sage matrix from this octave element.

        EXAMPLES::

            sage: A = octave('[1,2;3,4.5]')     # optional - octave
            sage: matrix(A)                     # optional - octave
            [1.0 2.0]
            [3.0 4.5]
            sage: _.base_ring()                 # optional - octave
            Real Double Field

            sage: A = octave('[I,1;-1,0]')      # optional - octave
            sage: matrix(A)                     # optional - octave
            [1.0*I   1.0]
            [ -1.0   0.0]
            sage: _.base_ring()                 # optional - octave
            Complex Double Field

            sage: # optional - octave
            sage: A = octave('[1,2;3,4]')
            sage: matrix(ZZ, A)
            [1 2]
            [3 4]
            sage: A = octave('[1,2;3,4.5]')
            sage: matrix(RR, A)
            [1.00000000000000 2.00000000000000]
            [3.00000000000000 4.50000000000000]
        """
        if not self.ismatrix():
            raise TypeError('not an octave matrix')
        if R is None:
            R = self._get_sage_ring()

        s = str(self).strip('\n ')
        w = [u.strip().split(' ') for u in s.split('\n')]
        nrows = len(w)
        ncols = len(w[0])

        if self.iscomplex():
            w = [[to_complex(x, R) for x in row] for row in w]

        from sage.matrix.matrix_space import MatrixSpace
        return MatrixSpace(R, nrows, ncols)(w)

    def _vector_(self, R=None):
        r"""
        Return Sage vector from this octave element.

        EXAMPLES::

            sage: # optional - octave
            sage: A = octave('[1,2,3,4]')
            sage: vector(ZZ, A)
            (1, 2, 3, 4)
            sage: A = octave('[1,2.3,4.5]')
            sage: vector(A)
            (1.0, 2.3, 4.5)
            sage: A = octave('[1,I]')
            sage: vector(A)
            (1.0, 1.0*I)
        """
        from sage.modules.free_module import FreeModule
        if not self.isvector():
            raise TypeError('not an octave vector')
        if R is None:
            R = self._get_sage_ring()

        s = str(self).strip('\n ')
        w = s.strip().split(' ')
        nrows = len(w)

        if self.iscomplex():
            w = [to_complex(x, R) for x in w]

        return FreeModule(R, nrows)(w)

    def _scalar_(self):
        """
        Return Sage scalar from this octave element.

        EXAMPLES::

            sage: A = octave('2833')      # optional - octave
            sage: As = A.sage(); As       # optional - octave
            2833.0
            sage: As.parent()             # optional - octave
            Real Double Field

            sage: B = sqrt(A)             # optional - octave
            sage: Bs = B.sage(); Bs       # optional - octave
            53.2259
            sage: Bs.parent()             # optional - octave
            Real Double Field

            sage: C = sqrt(-A)            # optional - octave
            sage: Cs = C.sage(); Cs       # optional - octave
            53.2259*I
            sage: Cs.parent()             # optional - octave
            Complex Double Field
        """
        if not self.isscalar():
            raise TypeError("not an octave scalar")

        R = self._get_sage_ring()
        if self.iscomplex():
            return to_complex(str(self), R)
        else:
            return R(str(self))

    def _sage_(self):
        """
        Try to parse the octave object and return a sage object.

        EXAMPLES::

            sage: # optional - octave
            sage: A = octave('2833')
            sage: A.sage()
            2833.0
            sage: B = sqrt(A)
            sage: B.sage()
            53.2259
            sage: C = sqrt(-A)
            sage: C.sage()
            53.2259*I
            sage: A = octave('[1,2,3,4]')
            sage: A.sage()
            (1.0, 2.0, 3.0, 4.0)
            sage: A = octave('[1,2.3,4.5]')
            sage: A.sage()
            (1.0, 2.3, 4.5)
            sage: A = octave('[1,2.3+I,4.5]')
            sage: A.sage()
            (1.0, 2.3 + 1.0*I, 4.5)
        """
        if self.isscalar():
            return self._scalar_()
        elif self.isvector():
            return self._vector_()
        elif self.ismatrix():
            return self._matrix_()
        else:
            raise NotImplementedError('octave type is not recognized')


# An instance
octave = Octave()


def reduce_load_Octave():
    """
    EXAMPLES::

        sage: from sage.interfaces.octave import reduce_load_Octave
        sage: reduce_load_Octave()
        Octave
    """
    return octave


def octave_console():
    """
    Spawn a new Octave command-line session.

    This requires that the optional octave program be installed and in
    your PATH, but no optional Sage packages need be installed.

    EXAMPLES::

        sage: octave_console()         # not tested
        GNU Octave, version 2.1.73 (i386-apple-darwin8.5.3).
        Copyright (C) 2006 John W. Eaton.
        ...
        octave:1> 2+3
        ans = 5
        octave:2> [ctl-d]

    Pressing ctrl-d exits the octave console and returns you to Sage.
    octave, like Sage, remembers its history from one session to
    another.
    """
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%octave magics instead.')
    os.system('octave-cli')
