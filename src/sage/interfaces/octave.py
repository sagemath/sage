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

    sage: octave("airy(3,2)")         # optional - octave
    4.10068
    sage: octave("beta(2,2)")         # optional - octave
    0.166667
    sage: octave("betainc(0.2,2,2)")  # optional - octave
    0.104
    sage: octave("besselh(0,2)")      # optional - octave
    (0.223891,0.510376)
    sage: octave("besselh(0,1)")      # optional - octave
    (0.765198,0.088257)
    sage: octave("besseli(1,2)")      # optional - octave
    1.59064
    sage: octave("besselj(1,2)")      # optional - octave
    0.576725
    sage: octave("besselk(1,2)")      # optional - octave
    0.139866
    sage: octave("erf(0)")            # optional - octave
    0
    sage: octave("erf(1)")            # optional - octave
    0.842701
    sage: octave("erfinv(0.842)")     # optional - octave
    0.998315
    sage: octave("gamma(1.5)")        # optional - octave
    0.886227
    sage: octave("gammainc(1.5,1)")   # optional - octave
    0.77687

The Octave interface reads in even very long input (using files) in
a robust manner::

    sage: t = '"%s"'%10^10000   # ten thousand character string.
    sage: a = octave.eval(t + ';')    # optional - octave, < 1/100th of a second
    sage: a = octave(t)               # optional - octave

Note that actually reading a back out takes forever. This *must*
be fixed ASAP - see
http://trac.sagemath.org/sage_trac/ticket/940/.

Tutorial
--------

EXAMPLES::

    sage: octave('4+10')              # optional - octave
    14
    sage: octave('date')              # optional - octave; random output
    18-Oct-2007
    sage: octave('5*10 + 6')          # optional - octave
    56
    sage: octave('(6+6)/3')           # optional - octave
    4
    sage: octave('9')^2               # optional - octave
    81
    sage: a = octave(10); b = octave(20); c = octave(30)    # optional - octave
    sage: avg = (a+b+c)/3             # optional - octave
    sage: avg                         # optional - octave
    20
    sage: parent(avg)                 # optional - octave
    Octave

::

    sage: my_scalar = octave('3.1415')       # optional - octave
    sage: my_scalar                          # optional - octave
    3.1415
    sage: my_vector1 = octave('[1,5,7]')     # optional - octave
    sage: my_vector1                         # optional - octave
    1     5     7
    sage: my_vector2 = octave('[1;5;7]')     # optional - octave
    sage: my_vector2                         # optional - octave
    1
    5
    7
    sage: my_vector1 * my_vector2            # optional - octave
    75
"""

#*****************************************************************************
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os
from expect import Expect, ExpectElement
from sage.misc.misc import verbose


class Octave(Expect):
    r"""
    Interface to the Octave interpreter.

    EXAMPLES::

        sage: octave.eval("a = [ 1, 1, 2; 3, 5, 8; 13, 21, 33 ]")    # optional - octave
        'a =\n\n 1 1 2\n 3 5 8\n 13 21 33\n\n'
        sage: octave.eval("b = [ 1; 3; 13]")                         # optional - octave
        'b =\n\n 1\n 3\n 13\n\n'
        sage: octave.eval("c=a \\ b") # solves linear equation: a*c = b  # optional - octave; random output
        'c =\n\n 1\n 7.21645e-16\n -7.21645e-16\n\n'
        sage: octave.eval("c")                                 # optional - octave; random output
        'c =\n\n 1\n 7.21645e-16\n -7.21645e-16\n\n'
    """

    def __init__(self, maxread=100, script_subdirectory=None, logfile=None, server=None, server_tmpdir=None):
        """
        EXAMPLES::

            sage: octave == loads(dumps(octave))
            True
        """
        Expect.__init__(self,
                        name = 'octave',
                        prompt = '>',
                        command = "sage-native-execute octave --no-line-editing --silent",
                        maxread = maxread,
                        server = server,
                        server_tmpdir = server_tmpdir,
                        script_subdirectory = script_subdirectory,
                        restart_on_ctrlc = False,
                        verbose_start = False,
                        logfile = logfile,
                        eval_using_file_cutoff=100)

    def __reduce__(self):
        """
        EXAMPLES::

            sage: octave.__reduce__()
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
        return 'source("%s");'%filename

    def _quit_string(self):
        """
        EXAMPLES::

            sage: octave._quit_string()
            'quit;'
        """
        return 'quit;'

    def _install_hints(self):
        """
        Returns hints on how to install Octave.

        EXAMPLES::

            sage: print octave._install_hints()
            You must get ...
        """
        return """
        You must get the program "octave" in order to use Octave
        from Sage.   You can read all about Octave at
                http://www.gnu.org/software/octave/

        LINUX / WINDOWS (colinux):
           Do apt-get install octave as root on your machine
           (or, in Windows, in the colinux console).

        OS X:
           * This website has links to binaries for OS X PowerPC
             and OS X Intel builds of the latest version of Octave:
                     http://hpc.sourceforge.net/
             Once you get the tarball from there, go to the / directory
             and type
                     tar zxvf octave-intel-bin.tar.gz
             to extract it to usr/local/...   Make sure /usr/local/bin
             is in your PATH.  Then type "octave" and verify that
             octave starts up.
           * Darwin ports and fink have Octave as well.
        """

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
        if not self._expect is None:
            if verbose:
                print "Exiting spawned %s process."%self
        return

    def _start(self):
        """
        Starts the Octave process.

        EXAMPLES::

            sage: o = Octave()    # optional - octave
            sage: o.is_running()  # optional - octave
            False
            sage: o._start()      # optional - octave
            sage: o.is_running()  # optional - octave
            True
        """
        Expect._start(self)
        self.eval("page_screen_output=0;")
        self.eval("format none;")

    def set(self, var, value):
        """
        Set the variable var to the given value.

        EXAMPLES::

            sage: octave.set('x', '2') # optional - octave
            sage: octave.get('x') # optional - octave
            ' 2'
        """
        cmd = '%s=%s;'%(var,value)
        out = self.eval(cmd)
        if out.find("error") != -1:
            raise TypeError("Error executing code in Octave\nCODE:\n\t%s\nOctave ERROR:\n\t%s"%(cmd, out))

    def get(self, var):
        """
        Get the value of the variable var.

        EXAMPLES::

            sage: octave.set('x', '2') # optional - octave
            sage: octave.get('x') # optional - octave
            ' 2'
        """
        s = self.eval('%s'%var)
        i = s.find('=')
        return s[i+1:]

    def clear(self, var):
        """
        Clear the variable named var.

        EXAMPLES::

            sage: octave.set('x', '2') # optional - octave
            sage: octave.clear('x') # optional - octave
            sage: octave.get('x') # optional - octave
            "error: `x' undefined near line ... column 1"
        """
        self.eval('clear %s'%var)

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
        """
        Return the version of Octave.

        OUTPUT: string

        EXAMPLES::

            sage: octave.version()   # optional - octave; random output depending on version
            '2.1.73'
        """
        return octave_version()

    def solve_linear_system(self, A, b):
        """
        Use octave to compute a solution x to A\*x = b, as a list.

        INPUT:


        -  ``A`` - mxn matrix A with entries in QQ or RR

        -  ``b`` - m-vector b entries in QQ or RR (resp)


        OUTPUT: An list x (if it exists) which solves M\*x = b

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
        n = A.ncols()
        if m != len(b):
            raise ValueError("dimensions of A and b must be compatible")
        from sage.matrix.all import MatrixSpace
        from sage.rings.all import QQ
        MS = MatrixSpace(QQ,m,1)
        b  = MS(list(b)) # converted b to a "column vector"
        sA = self.sage2octave_matrix_string(A)
        sb = self.sage2octave_matrix_string(b)
        self.eval("a = " + sA )
        self.eval("b = " + sb )
        soln = octave.eval("c = a \\ b")
        soln = soln.replace("\n\n ","[")
        soln = soln.replace("\n\n","]")
        soln = soln.replace("\n",",")
        sol  = soln[3:]
        return eval(sol)


    def sage2octave_matrix_string(self, A):
        """
        Return an octave matrix from a Sage matrix.

        INPUT: A Sage matrix with entries in the rationals or reals.

        OUTPUT: A string that evaluates to an Octave matrix.

        EXAMPLES::

            sage: M33 = MatrixSpace(QQ,3,3)
            sage: A = M33([1,2,3,4,5,6,7,8,0])
            sage: octave.sage2octave_matrix_string(A)   # optional - octave
            '[1, 2, 3; 4, 5, 6; 7, 8, 0]'

        AUTHORS:

        - David Joyner and William Stein
        """
        return str(A.rows()).replace('), (', '; ').replace('(', '').replace(')','')

    def de_system_plot(self, f, ics, trange):
        r"""
        Plots (using octave's interface to gnuplot) the solution to a
        `2\times 2` system of differential equations.

        INPUT:


        -  ``f`` - a pair of strings representing the
           differential equations; The independent variable must be called x
           and the dependent variable must be called y.

        -  ``ics`` - a pair [x0,y0] such that x(t0) = x0, y(t0)
           = y0

        -  ``trange`` - a pair [t0,t1]


        OUTPUT: a gnuplot window appears

        EXAMPLES::

            sage: octave.de_system_plot(['x+y','x-y'], [1,-1], [0,2])  # not tested -- does this actually work (on OS X it fails for me -- William Stein, 2007-10)

        This should yield the two plots `(t,x(t)), (t,y(t))` on the
        same graph (the `t`-axis is the horizontal axis) of the
        system of ODEs

        .. math::

                       x' = x+y, x(0) = 1;\qquad y' = x-y, y(0) = -1,                     \quad\text{for}\quad 0 < t < 2.
        """
        eqn1 = f[0].replace('x','x(1)').replace('y','x(2)')
        eqn2 = f[1].replace('x','x(1)').replace('y','x(2)')
        fcn = "function xdot = f(x,t) xdot(1) = %s; xdot(2) = %s; endfunction"%(eqn1, eqn2)
        self.eval(fcn)
        x0_eqn = "x0 = [%s; %s]"%(ics[0], ics[1])
        self.eval(x0_eqn)
        t_eqn = "t = linspace(%s, %s, 200)'"%(trange[0], trange[1])
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


class OctaveElement(ExpectElement):
    def _matrix_(self, R):
        r"""
        Return Sage matrix from this octave element.

        EXAMPLES::

            sage: A = octave('[1,2;3,4]')       # optional - octave
            sage: matrix(ZZ, A)                 # optional - octave
            [1 2]
            [3 4]
            sage: A = octave('[1,2;3,4.5]')     # optional - octave
            sage: matrix(RR, A)                 # optional - octave
            [1.00000000000000 2.00000000000000]
            [3.00000000000000 4.50000000000000]
        """
        from sage.matrix.all import MatrixSpace
        s = str(self).strip()
        v = s.split('\n ')
        nrows = len(v)
        if nrows == 0:
            return MatrixSpace(R,0,0)(0)
        ncols = len(v[0].split())
        M = MatrixSpace(R, nrows, ncols)
        v = sum([[x for x in w.split()] for w in v], [])
        return M(v)


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


import os
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
    os.system('octave')


def octave_version():
    """
    Return the version of Octave installed.

    EXAMPLES::

        sage: octave_version()    # optional - octave; and output is random
        '2.9.12'
    """
    return str(octave('version')).strip()
