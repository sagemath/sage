r"""
Abstract interface to Maxima

Maxima is a free GPL'd general purpose computer algebra system
whose development started in 1968 at MIT. It contains symbolic
manipulation algorithms, as well as implementations of special
functions, including elliptic functions and generalized
hypergeometric functions. Moreover, Maxima has implementations of
many functions relating to the invariant theory of the symmetric
group `S_n`. (However, the commands for group invariants,
and the corresponding Maxima documentation, are in French.) For many
links to Maxima documentation see
http://maxima.sourceforge.net/docs.shtml/.

AUTHORS:

- William Stein (2005-12): Initial version

- David Joyner: Improved documentation

- William Stein (2006-01-08): Fixed bug in parsing

- William Stein (2006-02-22): comparisons (following suggestion of
  David Joyner)

- William Stein (2006-02-24): *greatly* improved robustness by adding
  sequence numbers to IO bracketing in _eval_line

- Robert Bradshaw, Nils Bruin, Jean-Pierre Flori (2010,2011): Binary library
  interface

This is an abstract class implementing the functions shared between the Pexpect
and library interfaces to Maxima.
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
import re
import sys
import subprocess

from sage.env import DOT_SAGE, MAXIMA
COMMANDS_CACHE = '%s/maxima_commandlist_cache.sobj' % DOT_SAGE

from sage.cpython.string import bytes_to_str

from sage.misc.multireplace import multiple_replace
from sage.structure.richcmp import richcmp, rich_to_bool

from .interface import (Interface, InterfaceElement, InterfaceFunctionElement,
                        InterfaceFunction, AsciiArtString)
from sage.interfaces.tab_completion import ExtraTabCompletion
from sage.misc.instancedoc import instancedoc

# The Maxima "apropos" command, e.g., apropos(det) gives a list
# of all identifiers that begin in a certain way.  This could
# maybe be useful somehow... (?)  Also maxima has a lot for getting
# documentation from the system -- this could also be useful.


class MaximaAbstract(ExtraTabCompletion, Interface):
    r"""
    Abstract interface to Maxima.

    INPUT:

    - ``name`` -- string

    OUTPUT: the interface

    EXAMPLES:

    This class should not be instantiated directly,
    but through its subclasses Maxima (Pexpect interface)
    or MaximaLib (library interface)::

        sage: m = Maxima()
        sage: from sage.interfaces.maxima_abstract import MaximaAbstract
        sage: isinstance(m,MaximaAbstract)
        True
    """

    def __init__(self, name='maxima_abstract'):
        r"""
        Create an instance of an abstract interface to Maxima.
        See ``MaximaAbstract`` for full documentation.

        EXAMPLES::

            sage: from sage.interfaces.maxima_abstract import MaximaAbstract
            sage: isinstance(maxima,MaximaAbstract)
            True

        TESTS::

            sage: from sage.interfaces.maxima_abstract import MaximaAbstract
            sage: loads(dumps(MaximaAbstract)) == MaximaAbstract
            True
        """
        Interface.__init__(self, name)

    ###########################################
    # System -- change directory, etc.
    ###########################################
    def chdir(self, dir):
        r"""
        Change Maxima's current working directory.

        INPUT:

        - ``dir`` -- string

        OUTPUT: none

        EXAMPLES::

            sage: maxima.chdir('/')
        """
        self.lisp('(ext::cd "%s")' % dir)

    ###########################################
    # Interactive help
    ###########################################
    def _command_runner(self, command, s, redirect=True):
        r"""
        Run ``command`` in a new Maxima session and return its
        output as an ``AsciiArtString``.

        INPUT:

        - ``command`` -- string; function to call

        - ``s`` -- string; argument to the function

        - ``redirect`` -- boolean (default: ``True``); if redirect is set to False,
          then the output of the command is not returned as a string.
          Instead, it behaves like os.system. This is used for interactive
          things like Maxima's demos. See maxima.demo?

        OUTPUT:

        Output of ``command(s)`` as an ``AsciiArtString`` if ``redirect`` is set
        to False; None otherwise.

        EXAMPLES::

            sage: maxima._command_runner('describe', 'gcd')
            -- Function: gcd (<p_1>, <p_2>, <x_1>, ...)
            ...
        """
        cmd = '{} --very-quiet --batch-string="{}({});" '.format(MAXIMA, command, s)
        env = os.environ.copy()

        if redirect:
            res = bytes_to_str(subprocess.check_output(cmd, shell=True,
                                                       env=env))
            # We get a few lines of commented verbosity every time Maxima starts
            while res.startswith(';;;'):
                newline = res.find('\n')
                if newline == -1:
                    break
                res = res[newline + 1:]
            # The input is echoed, so we need to get rid of it
            res = res[res.find('\n')+1:]

            return AsciiArtString(res)
        else:
            subprocess.check_call(cmd, shell=True, env=env)

    def help(self, s):
        r"""
        Return Maxima's help for ``s``.

        INPUT:

        - ``s`` -- string

        OUTPUT: Maxima's help for ``s``

        EXAMPLES::

            sage: maxima.help('gcd')
            -- Function: gcd (<p_1>, <p_2>, <x_1>, ...)
            ...
        """
        # Should this be implemented without launching a new Maxima session
        # i.e. using eval_line ?
        return self._command_runner("describe", s)

    def example(self, s):
        r"""
        Return Maxima's examples for ``s``.

        INPUT:

        - ``s`` -- string

        OUTPUT: Maxima's examples for ``s``

        EXAMPLES::

            sage: maxima.example('arrays')
            a[n]:=n*a[n-1]
                                            a  := n a
                                             n       n - 1
            a[0]:1
            a[5]
                                                  120
            a[n]:=n
            a[6]
                                                   6
            a[4]
                                                  24
                                                 done
        """
        # Should this be implemented without launching a new Maxima session
        # i.e. using eval_line ?
        return self._command_runner("example", s)

    describe = help

    def demo(self, s):
        r"""
        Run Maxima's demo for ``s``.

        INPUT:

        - ``s`` -- string

        OUTPUT: none

        EXAMPLES::

            sage: maxima.demo('cf') # not tested
            read and interpret file: .../share/maxima/5.34.1/demo/cf.dem

            At the '_' prompt, type ';' and <enter> to get next demonstration.
            frac1:cf([1,2,3,4])
            ...
        """
        # Should this be implemented without launching a new Maxima session
        # i.e. using eval_line ?
        return self._command_runner("demo", s, redirect=False)

    def completions(self, s, verbose=True):
        r"""
        Return all commands that complete the command starting with the
        string ``s``. This is like typing s[tab] in the Maxima interpreter.

        INPUT:

        - ``s`` -- string

        - ``verbose`` -- boolean (default: ``True``)

        OUTPUT: array of strings

        EXAMPLES::

            sage: sorted(maxima.completions('gc', verbose=False))
            ['gcd', 'gcdex', 'gcfactor', 'gctime']
        """
        if verbose:
            print(s, end="")
            sys.stdout.flush()
        # in Maxima 5.19.1, apropos returns all commands that contain
        # the given string, instead of all commands that start with
        # the given string
        #
        # Maxima 5.44 changed DEFMFUN so that it creates both $NAME
        # and $NAME-IMPL (although the documentation suggests it would
        # create NAME-IMPL, without the leading $).  This causes
        # name-impl to show up in $APROPOS.  We remove it.
        # https://sourceforge.net/p/maxima/bugs/3643/
        cmd_list = self._eval_line('apropos("%s")' % s, error_check=False)
        cmd_list = cmd_list.replace(' ', '').replace('\n', '').replace('\\ - ', '-').replace('\\-', '-')
        cmd_list = [x for x in cmd_list[1:-1].split(',') if x[0] != '?' and not x.endswith('-impl')]
        return [x for x in cmd_list if x.find(s) == 0]

    def _commands(self, verbose=True):
        """
        Return list of all commands defined in Maxima.

        INPUT:

        - ``verbose`` -- boolean (default: ``True``)

        OUTPUT: array of strings

        EXAMPLES::

            # The output is kind of random
            sage: sorted(maxima._commands(verbose=False))
            [...
             'display',
             ...
             'gcd',
             ...
             'verbose',
             ...]
        """
        try:
            return self.__commands
        except AttributeError:
            self.__commands = sum(
                [self.completions(chr(65+n), verbose=verbose) +
                 self.completions(chr(97+n), verbose=verbose)
                 for n in range(26)], [])
        return self.__commands

    def _tab_completion(self, verbose=True, use_disk_cache=True):
        r"""
        Return all Maxima commands, which is useful for tab completion.

        INPUT:

        - ``verbose`` -- boolean (default: ``True``)

        - ``use_disk_cache`` -- boolean (default: ``True``); if set to True,
          try to read cached result from disk

        OUTPUT: array of strings

        EXAMPLES::

            sage: t = maxima._tab_completion(verbose=False)
            sage: 'gcd' in t
            True
            sage: len(t)    # random output
            1840
        """
        try:
            return self.__tab_completion
        except AttributeError:
            import sage.misc.persist
            if use_disk_cache:
                try:
                    self.__tab_completion = sage.misc.persist.load(COMMANDS_CACHE)
                    return self.__tab_completion
                except OSError:
                    pass
            if verbose:
                print("\nBuilding Maxima command completion list (this takes")
                print("a few seconds only the first time you do it).")
                print("To force rebuild later, delete %s." % COMMANDS_CACHE)
            v = self._commands(verbose=verbose)
            if verbose:
                print("\nDone!")
            self.__tab_completion = v
            sage.misc.persist.save(v, COMMANDS_CACHE)
            return v

    def console(self):
        r"""
        Start the interactive Maxima console. This is a completely separate
        maxima session from this interface. To interact with this session,
        you should instead use ``maxima.interact()``.

        OUTPUT: none

        EXAMPLES::

            sage: maxima.console()             # not tested (since we can't)
            Maxima 5.46.0 https://maxima.sourceforge.io
            using Lisp ECL 21.2.1
            Distributed under the GNU Public License. See the file COPYING.
            Dedicated to the memory of William Schelter.
            This is a development version of Maxima. The function bug_report()
            provides bug reporting information.
            (%i1)

        ::

            sage: maxima.interact()     # not tested
              --> Switching to Maxima <--
            maxima: 2+2
            4
            maxima:
              --> Exiting back to Sage <--
        """
        maxima_console()

    def cputime(self, t=None):
        r"""
        Return the amount of CPU time that this Maxima session has used.

        INPUT:

        - ``t`` -- float (default: ``None``); if \var{t} is not None, then
          it returns the difference between the current CPU time and \var{t}

        OUTPUT: float

        EXAMPLES::

            sage: t = maxima.cputime()
            sage: _ = maxima.de_solve('diff(y,x,2) + 3*x = y', ['x','y'], [1,1,1])
            sage: maxima.cputime(t) # output random
            0.568913
        """
        if t:
            return float(self.eval('elapsed_run_time()')) - t
        else:
            return float(self.eval('elapsed_run_time()'))

    def version(self):
        r"""
        Return the version of Maxima that Sage includes.

        OUTPUT: none

        EXAMPLES::

            sage: maxima.version()  # random
            '5.41.0'
        """
        return maxima_version()

    ####
    # Overriding default values
    ###

    def _assign_symbol(self):
        r"""
        Return the assign symbol in Maxima.

        OUTPUT: string

        EXAMPLES::

            sage: maxima._assign_symbol()
            ':'
            sage: maxima.eval('t : 8')
            '8'
            sage: maxima.eval('t')
            '8'
        """
        return ":"

    def _true_symbol(self):
        """
        Return the true symbol in Maxima.

        OUTPUT: string

        EXAMPLES::

            sage: maxima._true_symbol()
            'true'
            sage: maxima.eval('is(2 = 2)')
            'true'
        """
        return 'true'

    def _false_symbol(self):
        """
        Return the false symbol in Maxima.

        OUTPUT: string

        EXAMPLES::

            sage: maxima._false_symbol()
            'false'
            sage: maxima.eval('is(2 = 4)')
            'false'
        """
        return 'false'

    def _equality_symbol(self):
        """
        Return the equality symbol in Maxima.

        OUTPUT: string

        EXAMPLES::

             sage: maxima._equality_symbol()
             '='
             sage: var('x y')
             (x, y)
             sage: maxima(x == y)
             _SAGE_VAR_x = _SAGE_VAR_y
        """
        return '='

    def _inequality_symbol(self):
        """
        Return the inequality symbol in Maxima.

        OUTPUT: string

        EXAMPLES::

             sage: maxima._inequality_symbol()
             '#'
             sage: maxima((x != 1))
             _SAGE_VAR_x # 1
        """
        return '#'

    def _function_class(self):
        """
        Return the Python class of Maxima functions.

        OUTPUT: type

        EXAMPLES::

            sage: maxima._function_class()
            <class 'sage.interfaces.interface.InterfaceFunction'>
        """
        return MaximaAbstractFunction

    def _object_class(self):
        """
        Return the Python class of Maxima elements.

        OUTPUT: type

        EXAMPLES::

            sage: maxima._object_class()
            <class 'sage.interfaces.maxima.MaximaElement'>
        """
        return MaximaAbstractElement

    def _function_element_class(self):
        """
        Return the Python class of Maxima functions of elements.

        OUTPUT: type

        EXAMPLES::

            sage: maxima._function_element_class()
            <class 'sage.interfaces.interface.InterfaceFunctionElement'>
        """
        return MaximaAbstractFunctionElement

    def _object_function_class(self):
        """
        Return the Python class of Maxima user-defined functions.

        OUTPUT: type

        EXAMPLES::

            sage: maxima._object_function_class()
            <class 'sage.interfaces.maxima.MaximaElementFunction'>
        """
        return MaximaAbstractElementFunction

    ####################
    # Maxima functions #
    ####################

    def function(self, args, defn, rep=None, latex=None):
        """
        Return the Maxima function with given arguments and definition.

        INPUT:

        - ``args`` -- string with variable names separated by
          commas

        - ``defn`` -- string (or Maxima expression) that
          defines a function of the arguments in Maxima

        - ``rep`` -- an optional string; if given, this is how
          the function will print

        OUTPUT: Maxima function

        EXAMPLES::

            sage: f = maxima.function('x', 'sin(x)')
            sage: f(3.2)  # abs tol 2e-16
            -0.058374143427579909
            sage: f = maxima.function('x,y', 'sin(x)+cos(y)')
            sage: f(2, 3.5)  # abs tol 2e-16
            sin(2)-0.9364566872907963
            sage: f
            sin(x)+cos(y)

        ::

            sage: g = f.integrate('z')
            sage: g
            (cos(y)+sin(x))*z
            sage: g(1,2,3)
            3*(cos(2)+sin(1))

        The function definition can be a Maxima object::

            sage: an_expr = maxima('sin(x)*gamma(x)')
            sage: t = maxima.function('x', an_expr)
            sage: t
            gamma(x)*sin(x)
            sage: t(2)
             sin(2)
            sage: float(t(2))
            0.9092974268256817
            sage: loads(t.dumps())
            gamma(x)*sin(x)
        """
        name = self._next_var_name()
        if isinstance(defn, MaximaAbstractElement):
            defn = defn.str()
        elif not isinstance(defn, str):
            defn = str(defn)
        if isinstance(args, MaximaAbstractElement):
            args = args.str()
        elif not isinstance(args, str):
            args = str(args)
        cmd = '%s(%s) := %s' % (name, args, defn)
        self._eval_line(cmd)
        if rep is None:
            rep = defn
        f = self._object_function_class()(self, name, rep, args, latex)
        return f

#     def display2d(self, flag=True):
#         """
#         Set the flag that determines whether Maxima objects are
#         printed using their 2-d ASCII art representation.  When the
#         maxima interface starts the default is that objects are not
#         represented in 2-d.

#         INPUT:

#         flag -- boolean (default: ``True``)

#         EXAMPLES::

#             sage: maxima('1/2')
#             1/2
#             sage: maxima.display2d(True)
#             sage: maxima('1/2')
#                                            1
#                                            -
#                                            2
#             sage: maxima.display2d(False)
#         """
#         self._display2d = bool(flag)

    def plot2d(self, *args):
        r"""
        Plot a 2d graph using Maxima / gnuplot.

        maxima.plot2d(f, '[var, min, max]', options)

        INPUT:

        - ``f`` -- string representing a function (such as
          f="sin(x)") [var, xmin, xmax]

        - ``options`` -- an optional string representing plot2d
          options in gnuplot format

        EXAMPLES::

            sage: maxima.plot2d('sin(x)','[x,-5,5]')   # not tested
            sage: opts = '[gnuplot_term, ps], [gnuplot_out_file, "sin-plot.eps"]'
            sage: maxima.plot2d('sin(x)','[x,-5,5]',opts)    # not tested

        The eps file is saved in the current directory.
        """
        self('plot2d(%s)' % (','.join(str(x) for x in args)))

    def plot2d_parametric(self, r, var, trange, nticks=50, options=None):
        r"""
        Plot r = [x(t), y(t)] for t = tmin...tmax using gnuplot with
        options.

        INPUT:

        - ``r`` -- string representing a function (such as
          r="[x(t),y(t)]")

        - ``var`` -- string representing the variable (such
          as var = "t")

        - ``trange`` -- [tmin, tmax] are numbers with tmintmax

        - ``nticks`` -- integer (default: 50)

        - ``options`` -- an optional string representing plot2d
          options in gnuplot format

        EXAMPLES::

            sage: maxima.plot2d_parametric(["sin(t)","cos(t)"], "t",[-3.1,3.1])   # not tested

        ::

            sage: opts = '[gnuplot_preamble, "set nokey"], [gnuplot_term, ps], [gnuplot_out_file, "circle-plot.eps"]'
            sage: maxima.plot2d_parametric(["sin(t)","cos(t)"], "t", [-3.1,3.1], options=opts)   # not tested

        The eps file is saved to the current working directory.

        Here is another fun plot::

            sage: maxima.plot2d_parametric(["sin(5*t)","cos(11*t)"], "t", [0,2*pi()], nticks=400)    # not tested
        """
        tmin = trange[0]
        tmax = trange[1]
        cmd = "plot2d([parametric, %s, %s, [%s, %s, %s], [nticks, %s]]" % (
            r[0], r[1], var, tmin, tmax, nticks)
        if options is None:
            cmd += ")"
        else:
            cmd += ", %s)" % options
        self(cmd)

    def plot3d(self, *args):
        r"""
        Plot a 3d graph using Maxima / gnuplot.

        maxima.plot3d(f, '[x, xmin, xmax]', '[y, ymin, ymax]', '[grid, nx,
        ny]', options)

        INPUT:

        - ``f`` -- string representing a function (such as
          f="sin(x)") [var, min, max]

        - ``args`` should be of the form '[x, xmin, xmax]', '[y, ymin, ymax]',
          '[grid, nx, ny]', options

        EXAMPLES::

            sage: maxima.plot3d('1 + x^3 - y^2', '[x,-2,2]', '[y,-2,2]', '[grid,12,12]')    # not tested
            sage: maxima.plot3d('sin(x)*cos(y)', '[x,-2,2]', '[y,-2,2]', '[grid,30,30]')   # not tested
            sage: opts = '[gnuplot_term, ps], [gnuplot_out_file, "sin-plot.eps"]'
            sage: maxima.plot3d('sin(x+y)', '[x,-5,5]', '[y,-1,1]', opts)    # not tested

        The eps file is saved in the current working directory.
        """
        self('plot3d(%s)' % (','.join(str(x) for x in args)))

    def plot3d_parametric(self, r, vars, urange, vrange, options=None):
        r"""
        Plot a 3d parametric graph with r=(x,y,z), x = x(u,v), y = y(u,v),
        z = z(u,v), for u = umin...umax, v = vmin...vmax using gnuplot with
        options.

        INPUT:

        - ``x``, ``y``, ``z`` -- string representing a function (such
          as ``x="u2+v2"``, ...) vars is a list or two strings
          representing variables (such as vars = ["u","v"])

        - ``urange`` -- [umin, umax]

        - ``vrange`` -- [vmin, vmax] are lists of numbers with
          umin umax, vmin vmax

        - ``options`` -- (optional) string representing plot2d
          options in gnuplot format

        OUTPUT: displays a plot on screen or saves to a file

        EXAMPLES::

            sage: maxima.plot3d_parametric(["v*sin(u)","v*cos(u)","v"], ["u","v"],[-3.2,3.2],[0,3])     # not tested
            sage: opts = '[gnuplot_term, ps], [gnuplot_out_file, "sin-cos-plot.eps"]'
            sage: maxima.plot3d_parametric(["v*sin(u)","v*cos(u)","v"], ["u","v"],[-3.2,3.2],[0,3],opts)      # not tested

        The eps file is saved in the current working directory.

        Here is a torus::

            sage: _ = maxima.eval("expr_1: cos(y)*(10.0+6*cos(x)); expr_2: sin(y)*(10.0+6*cos(x)); expr_3: -6*sin(x);")
            sage: maxima.plot3d_parametric(["expr_1","expr_2","expr_3"], ["x","y"],[0,6],[0,6])  # not tested

        Here is a Möbius strip::

            sage: x = "cos(u)*(3 + v*cos(u/2))"
            sage: y = "sin(u)*(3 + v*cos(u/2))"
            sage: z = "v*sin(u/2)"
            sage: maxima.plot3d_parametric([x,y,z],["u","v"],[-3.1,3.2],[-1/10,1/10])   # not tested
        """
        umin = urange[0]
        umax = urange[1]
        vmin = vrange[0]
        vmax = vrange[1]
        cmd = 'plot3d([%s, %s, %s], [%s, %s, %s], [%s, %s, %s]' % (
            r[0], r[1], r[2], vars[0], umin, umax, vars[1], vmin, vmax)
        if options is None:
            cmd += ')'
        else:
            cmd += ', %s)' % options
        self(cmd)

    def de_solve(self, de, vars, ics=None):
        """
        Solve a 1st or 2nd order ordinary differential equation (ODE) in
        two variables, possibly with initial conditions.

        INPUT:

        - ``de`` -- string representing the ODE

        - ``vars`` -- list of strings representing the two
          variables

        - ``ics`` -- a triple of numbers [a,b1,b2] representing
          y(a)=b1, y'(a)=b2

        EXAMPLES::

            sage: maxima.de_solve('diff(y,x,2) + 3*x = y', ['x','y'], [1,1,1])
            y = 3*x-2*%e^(x-1)
            sage: maxima.de_solve('diff(y,x,2) + 3*x = y', ['x','y'])
            y = %k1*%e^x+%k2*%e^-x+3*x
            sage: maxima.de_solve('diff(y,x) + 3*x = y', ['x','y'])
            y = (%c-3*(...-x...-1)*%e^-x)*%e^x
            sage: maxima.de_solve('diff(y,x) + 3*x = y', ['x','y'],[1,1])
            y = -...%e^-1*(5*%e^x-3*%e*x-3*%e)...
        """
        if not isinstance(vars, str):
            str_vars = '%s, %s' % (vars[1], vars[0])
        else:
            str_vars = vars
        self.eval('depends(%s)' % str_vars)
        m = self(de)
        a = 'ode2(%s, %s)' % (m.name(), str_vars)
        if ics is not None:
            if len(ics) == 3:
                cmd = "ic2("+a+",%s=%s,%s=%s,diff(%s,%s)=%s);" % (vars[0], ics[0], vars[1], ics[1], vars[1], vars[0], ics[2])
                return self(cmd)
            if len(ics) == 2:
                return self("ic1("+a+",%s=%s,%s=%s);" % (vars[0], ics[0],
                                                         vars[1], ics[1]))
        return self(a + ";")

    def de_solve_laplace(self, de, vars, ics=None):
        """
        Solve an ordinary differential equation (ODE) using Laplace
        transforms.

        INPUT:

        - ``de`` -- string representing the ODE (e.g., de =
          "diff(f(x),x,2)=diff(f(x),x)+sin(x)")

        - ``vars`` -- list of strings representing the
          variables (e.g., ``vars = ["x","f"]``)

        - ``ics`` -- list of numbers representing initial
          conditions, with symbols allowed which are represented by strings
          (eg, f(0)=1, f'(0)=2 is ics = [0,1,2])

        EXAMPLES::

            sage: maxima.clear('x'); maxima.clear('f')
            sage: maxima.de_solve_laplace("diff(f(x),x,2) = 2*diff(f(x),x)-f(x)", ["x","f"], [0,1,2])
            f(x) = x*%e^x+%e^x

        ::

            sage: maxima.clear('x'); maxima.clear('f')
            sage: f = maxima.de_solve_laplace("diff(f(x),x,2) = 2*diff(f(x),x)-f(x)", ["x","f"])
            sage: f
            f(x) = x*%e^x*('at('diff(f(x),x,1),x = 0))-f(0)*x*%e^x+f(0)*%e^x
            sage: print(f)
                                               !
                                   x  d        !                  x          x
                        f(x) = x %e  (-- (f(x))!     ) - f(0) x %e  + f(0) %e
                                      dx       !
                                               !x = 0

        .. NOTE::

           The second equation sets the values of `f(0)` and
           `f'(0)` in Maxima, so subsequent ODEs involving these
           variables will have these initial conditions automatically
           imposed.
        """
        if ics is not None:
            d = len(ics)
            for i in range(d - 1):
                ic = 'atvalue(diff(%s(%s), %s, %s), %s = %s, %s)' % (
                    vars[1], vars[0], vars[0], i, vars[0], ics[0], ics[1+i])
                self.eval(ic)
        return self('desolve(%s, %s(%s))' % (de, vars[1], vars[0]))

    def solve_linear(self, eqns, vars):
        """
        Wraps maxima's linsolve.

        INPUT:

        - ``eqns`` -- list of m strings; each representing a linear
          question in m = n variables

        - ``vars`` -- list of n strings; each
          representing a variable

        EXAMPLES::

            sage: eqns = ["x + z = y","2*a*x - y = 2*a^2","y - 2*z = 2"]
            sage: vars = ["x","y","z"]
            sage: maxima.solve_linear(eqns, vars)
            [x = a+1,y = 2*a,z = a-1]
        """
        eqs = "["
        for i in range(len(eqns)):
            if i < len(eqns) - 1:
                eqs = eqs + eqns[i] + ","
            if i == len(eqns) - 1:
                eqs = eqs + eqns[i] + "]"
        vrs = "["
        for i in range(len(vars)):
            if i < len(vars) - 1:
                vrs = vrs + vars[i] + ","
            if i == len(vars) - 1:
                vrs = vrs + vars[i] + "]"
        return self('linsolve(%s, %s)' % (eqs, vrs))

    def unit_quadratic_integer(self, n):
        r"""
        Finds a unit of the ring of integers of the quadratic number field
        `\QQ(\sqrt{n})`, `n>1`, using the qunit maxima command.

        INPUT:

        - ``n`` -- integer

        EXAMPLES::

            sage: u = maxima.unit_quadratic_integer(101); u
            a + 10
            sage: u.parent()
            Number Field in a with defining polynomial x^2 - 101 with a = 10.04987562112089?
            sage: u = maxima.unit_quadratic_integer(13)
            sage: u
            5*a + 18
            sage: u.parent()
            Number Field in a with defining polynomial x^2 - 13 with a = 3.605551275463990?
        """
        from sage.rings.integer import Integer
        from sage.rings.number_field.number_field import QuadraticField
        # Take square-free part so sqrt(n) doesn't get simplified
        # further by maxima
        # (The original version of this function would yield wrong answers if
        # n is not squarefree.)
        n = Integer(n).squarefree_part()
        if n < 1:
            raise ValueError("n (=%s) must be >= 1" % n)
        s = repr(self('qunit(%s)' % n)).lower()
        r = re.compile(r'sqrt\(.*\)')
        Qa = QuadraticField(n, 'a')
        s = r.sub('a', s)
        return Qa(s)

    def plot_list(self, ptsx, ptsy, options=None):
        r"""
        Plots a curve determined by a sequence of points.

        INPUT:

        - ``ptsx`` -- [x1,...,xn], where the xi and yi are
          real,

        - ``ptsy`` -- [y1,...,yn]

        - ``options`` -- string representing maxima plot2d
          options

        The points are (x1,y1), (x2,y2), etc.

        This function requires maxima 5.9.2 or newer.

        .. NOTE::

           More that 150 points can sometimes lead to the program
           hanging. Why?

        EXAMPLES::

            sage: zeta_ptsx = [(pari(1/2 + i*I/10).zeta().real()).precision(1)          # needs sage.libs.pari
            ....:              for i in range(70,150)]
            sage: zeta_ptsy = [(pari(1/2 + i*I/10).zeta().imag()).precision(1)          # needs sage.libs.pari
            ....:              for i in range(70,150)]
            sage: maxima.plot_list(zeta_ptsx, zeta_ptsy)        # not tested            # needs sage.libs.pari
            sage: opts='[gnuplot_preamble, "set nokey"], [gnuplot_term, ps], [gnuplot_out_file, "zeta.eps"]'
            sage: maxima.plot_list(zeta_ptsx, zeta_ptsy, opts)  # not tested            # needs sage.libs.pari
        """
        cmd = 'plot2d([discrete,%s, %s]' % (ptsx, ptsy)
        if options is None:
            cmd += ')'
        else:
            cmd += ', %s)' % options
        self(cmd)

    def plot_multilist(self, pts_list, options=None):
        r"""
        Plots a list of list of points pts_list=[pts1,pts2,...,ptsn],
        where each ptsi is of the form [[x1,y1],...,[xn,yn]] x's must be
        integers and y's reals options is a string representing maxima
        plot2d options.

        INPUT:

        - ``pts_lst`` -- list of points; each point must be of the form [x,y]
          where ``x`` is an integer and ``y`` is a real

        - ``var`` -- string; representing Maxima's plot2d options

        Requires maxima 5.9.2 at least.

        .. NOTE::

           More that 150 points can sometimes lead to the program
           hanging.

        EXAMPLES::

            sage: xx = [i/10.0 for i in range(-10,10)]
            sage: yy = [i/10.0 for i in range(-10,10)]
            sage: x0 = [0 for i in range(-10,10)]
            sage: y0 = [0 for i in range(-10,10)]
            sage: zeta_ptsx1 = [(pari(1/2+i*I/10).zeta().real()).precision(1)           # needs sage.libs.pari
            ....:               for i in range(10)]
            sage: zeta_ptsy1 = [(pari(1/2+i*I/10).zeta().imag()).precision(1)           # needs sage.libs.pari
            ....:               for i in range(10)]
            sage: maxima.plot_multilist([[zeta_ptsx1,zeta_ptsy1], [xx,y0], [x0,yy]])    # not tested
            sage: zeta_ptsx1 = [(pari(1/2+i*I/10).zeta().real()).precision(1)           # needs sage.libs.pari
            ....:               for i in range(10,150)]
            sage: zeta_ptsy1 = [(pari(1/2+i*I/10).zeta().imag()).precision(1)           # needs sage.libs.pari
            ....:               for i in range(10,150)]
            sage: maxima.plot_multilist([[zeta_ptsx1,zeta_ptsy1], [xx,y0], [x0,yy]])    # not tested
            sage: opts='[gnuplot_preamble, "set nokey"]'
            sage: maxima.plot_multilist([[zeta_ptsx1,zeta_ptsy1], [xx,y0], [x0,yy]],    # not tested
            ....:                       opts)
        """
        n = len(pts_list)
        cmd = '['
        for i in range(n):
            if i < n-1:
                cmd = cmd+'[discrete,'+str(pts_list[i][0])+','+str(pts_list[i][1])+'],'
            if i == n-1:
                cmd = cmd+'[discrete,'+str(pts_list[i][0])+','+str(pts_list[i][1])+']]'
        if options is None:
            self('plot2d('+cmd+')')
        else:
            self('plot2d('+cmd+','+options+')')


@instancedoc
class MaximaAbstractElement(ExtraTabCompletion, InterfaceElement):
    r"""
    Element of Maxima through an abstract interface.

    EXAMPLES:

    Elements of this class should not be created directly.
    The targeted parent of a concrete inherited class should be used instead::

        sage: from sage.interfaces.maxima_lib import maxima_lib
        sage: xp = maxima(x)
        sage: type(xp)
        <class 'sage.interfaces.maxima.MaximaElement'>
        sage: xl = maxima_lib(x)
        sage: type(xl)
        <class 'sage.interfaces.maxima_lib.MaximaLibElement'>
    """
    _cached_repr = True

    def __str__(self):
        """
        Printing an object explicitly gives ASCII art.

        OUTPUT: string

        EXAMPLES::

            sage: f = maxima('1/(x-1)^3'); f
            1/(x-1)^3
            sage: print(f)
                                                  1
                                               --------
                                                      3
                                               (x - 1)
        """
        return self.display2d(onscreen=False)

    def __bool__(self):
        """
        Convert ``self`` into a boolean.

        OUTPUT: boolean

        EXAMPLES::

            sage: bool(maxima(0))
            False
            sage: bool(maxima(1))
            True
            sage: bool(maxima('false'))
            False
            sage: bool(maxima('true'))
            True
        """
        P = self._check_valid()
        return (P.eval('is({0} = 0 or {0} = false);'.format(self.name()))
                != P._true_symbol())
        # but be careful, since for relations things like is(equal(a,b)) are
        # what Maxima needs

    def _richcmp_(self, other, op):
        """
        Compare this Maxima object with ``other``.

        INPUT:

        - ``other`` -- an object to compare to

        OUTPUT: boolean

        EXAMPLES::

            sage: a = maxima(1); b = maxima(2)
            sage: a == b
            False
            sage: a < b
            True
            sage: a > b
            False
            sage: b < a
            False
            sage: b > a
            True

        We can also compare more complicated object such as functions::

            sage: f = maxima('sin(x)'); g = maxima('cos(x)')
            sage: -f == g.diff('x')
            True
        """
        # thanks to David Joyner for telling me about using "is".
        # but be careful, since for relations things like is(equal(a,b))
        # are what Maxima needs
        P = self.parent()
        try:
            if P.eval("is (%s < %s)" % (self.name(), other.name())) == P._true_symbol():
                return rich_to_bool(op, -1)
            elif P.eval("is (%s > %s)" % (self.name(), other.name())) == P._true_symbol():
                return rich_to_bool(op, 1)
            elif P.eval("is (%s = %s)" % (self.name(), other.name())) == P._true_symbol():
                return rich_to_bool(op, 0)
        except TypeError:
            pass
        return richcmp(repr(self), repr(other), op)
        # everything is supposed to be comparable in Python2,
        # so we define the comparison thus when no comparable
        # in interfaced system.

    def _sage_(self):
        """
        Attempt to make a native Sage object out of this Maxima object.
        This is useful for automatic coercions in addition to other
        things.

        OUTPUT: Sage object

        EXAMPLES::

            sage: a = maxima('sqrt(2) + 2.5'); a
            sqrt(2)+2.5
            sage: b = a._sage_(); b
            sqrt(2) + 2.5
            sage: type(b)
            <class 'sage.symbolic.expression.Expression'>

        We illustrate an automatic coercion::

            sage: c = b + sqrt(3); c
            sqrt(3) + sqrt(2) + 2.5
            sage: type(c)
            <class 'sage.symbolic.expression.Expression'>
            sage: d = sqrt(3) + b; d
            sqrt(3) + sqrt(2) + 2.5
            sage: type(d)
            <class 'sage.symbolic.expression.Expression'>

            sage: a = sage.calculus.calculus.maxima('x^(sqrt(y)+%pi) + sin(%e + %pi)')
            sage: a._sage_()
            x^(pi + sqrt(y)) - sin(e)
            sage: var('x, y')
            (x, y)
            sage: v = sage.calculus.calculus.maxima.vandermonde_matrix([x, y, 1/2])
            sage: v._sage_()
            [  1   x x^2]
            [  1   y y^2]
            [  1 1/2 1/4]

        TESTS:

        Check if :issue:`7661` is fixed::

            sage: var('delta')
            delta
            sage: (2*delta).simplify()
            2*delta

        Check conversion of Booleans (:issue:`28705`)::

            sage: maxima('true')._sage_(), maxima('false')._sage_()
            (True, False)
        """
        import sage.calculus.calculus as calculus
        return calculus.symbolic_expression_from_maxima_string(self.name(),
                maxima=self.parent())

    def _symbolic_(self, R):
        """
        Return a symbolic expression equivalent to this Maxima object.

        INPUT:

        - ``R`` -- symbolic ring to convert into

        OUTPUT: symbolic expression

        EXAMPLES::

            sage: t = sqrt(2)._maxima_()
            sage: u = t._symbolic_(SR); u
            sqrt(2)
            sage: u.parent()
            Symbolic Ring

        This is used when converting Maxima objects to the Symbolic Ring::

            sage: SR(t)
            sqrt(2)
        """
        return R(self._sage_())

    def __complex__(self):
        """
        Return a complex number equivalent to this Maxima object.

        OUTPUT: complex

        EXAMPLES::

            sage: complex(maxima('sqrt(-2)+1'))
            (1+1.4142135623730951j)
        """
        return complex(self._sage_())

    def _complex_mpfr_field_(self, C):
        """
        Return a mpfr complex number equivalent to this Maxima object.

        INPUT:

        - ``C`` -- complex numbers field to convert into

        OUTPUT: complex

        EXAMPLES::

            sage: CC(maxima('1+%i'))
             1.00000000000000 + 1.00000000000000*I
            sage: CC(maxima('2342.23482943872+234*%i'))
             2342.23482943872 + 234.000000000000*I
            sage: ComplexField(10)(maxima('2342.23482943872+234*%i'))
             2300. + 230.*I
            sage: ComplexField(200)(maxima('1+%i'))
            1.0000000000000000000000000000000000000000000000000000000000 + 1.0000000000000000000000000000000000000000000000000000000000*I
            sage: ComplexField(200)(maxima('sqrt(-2)'))
            1.4142135623730950488016887242096980785696718753769480731767*I
            sage: N(sqrt(-2), 200)
            8.0751148893563733350506651837615871941533119425962889089783e-62 + 1.4142135623730950488016887242096980785696718753769480731767*I
        """
        return C(self._sage_())

    def _mpfr_(self, R):
        """
        Return a mpfr real number equivalent to this Maxima object.

        INPUT:

        - ``R`` -- real numbers field to convert into

        OUTPUT: real

        EXAMPLES::

            sage: RealField(100)(maxima('sqrt(2)+1'))
            2.4142135623730950488016887242
        """
        return R(self._sage_())

    def _complex_double_(self, C):
        """
        Return a double precision complex number equivalent to this Maxima object.

        INPUT:

        - ``C`` -- double precision complex numbers field to convert into

        OUTPUT: complex

        EXAMPLES::

            sage: CDF(maxima('sqrt(2)+1'))
            2.414213562373095
        """
        return C(self._sage_())

    def _real_double_(self, R):
        """
        Return a double precision real number equivalent to this Maxima object.

        INPUT:

        - ``R`` -- double precision real numbers field to convert into

        OUTPUT: real

        EXAMPLES::

            sage: RDF(maxima('sqrt(2)+1'))
            2.414213562373095
        """
        return R(self._sage_())

    def real(self):
        """
        Return the real part of this Maxima element.

        OUTPUT: Maxima real

        EXAMPLES::

            sage: maxima('2 + (2/3)*%i').real()
            2
        """
        return self.realpart()

    def imag(self):
        """
        Return the imaginary part of this Maxima element.

        OUTPUT: Maxima real

        EXAMPLES::

            sage: maxima('2 + (2/3)*%i').imag()
            2/3
        """
        return self.imagpart()

    def numer(self):
        """
        Return numerical approximation to ``self`` as a Maxima object.

        OUTPUT: Maxima object

        EXAMPLES::

            sage: a = maxima('sqrt(2)').numer(); a
            1.41421356237309...
            sage: type(a)
            <class 'sage.interfaces.maxima.MaximaElement'>
        """
        return self.comma('numer')

    def str(self):
        """
        Return string representation of this Maxima object.

        OUTPUT: string

        EXAMPLES::

            sage: maxima('sqrt(2) + 1/3').str()
            'sqrt(2)+1/3'
        """
        P = self._check_valid()
        return P.get(self._name)

    def diff(self, var='x', n=1):
        """
        Return the `n`-th derivative of ``self``.

        INPUT:

        - ``var`` -- variable (default: ``'x'``)

        - ``n`` -- integer (default: 1)

        OUTPUT: `n`-th derivative of ``self`` with respect to the variable var

        EXAMPLES::

            sage: f = maxima('x^2')
            sage: f.diff()
            2*x
            sage: f.diff('x')
            2*x
            sage: f.diff('x', 2)
            2
            sage: maxima('sin(x^2)').diff('x',4)
            16*x^4*sin(x^2)-12*sin(x^2)-48*x^2*cos(x^2)

        ::

            sage: f = maxima('x^2 + 17*y^2')
            sage: f.diff('x')
            34*y*'diff(y,x,1)+2*x
            sage: f.diff('y')
            34*y
        """
        return InterfaceElement.__getattr__(self, 'diff')(var, n)

    derivative = diff

    def nintegral(self, var='x', a=0, b=1,
                  desired_relative_error='1e-8',
                  maximum_num_subintervals=200):
        r"""
        Return a numerical approximation to the integral of ``self`` from `a`
        to `b`.

        INPUT:

        - ``var`` -- variable to integrate with respect to

        - ``a`` -- lower endpoint of integration

        - ``b`` -- upper endpoint of integration

        - ``desired_relative_error`` -- (default: ``'1e-8'``) the
          desired relative error

        - ``maximum_num_subintervals`` -- (default: 200)
          maxima number of subintervals

        OUTPUT: approximation to the integral

        - estimated absolute error of the
          approximation

        - the number of integrand evaluations

        - an error code:

            - ``0`` -- no problems were encountered

            - ``1`` -- too many subintervals were done

            - ``2`` -- excessive roundoff error

            - ``3`` -- extremely bad integrand behavior

            - ``4`` -- failed to converge

            - ``5`` -- integral is probably divergent or slowly convergent

            - ``6`` -- the input is invalid

        EXAMPLES::

            sage: maxima('exp(-sqrt(x))').nintegral('x',0,1)
            (0.5284822353142306, 4.163...e-11, 231, 0)

        Note that GP also does numerical integration, and can do so to very
        high precision very quickly::

            sage: gp('intnum(x=0,1,exp(-sqrt(x)))')
            0.5284822353142307136179049194             # 32-bit
            0.52848223531423071361790491935415653022   # 64-bit
            sage: _ = gp.set_precision(80)
            sage: gp('intnum(x=0,1,exp(-sqrt(x)))')
            0.52848223531423071361790491935415653021675547587292866196865279321015401702040079
        """
        from sage.rings.integer import Integer
        v = self.quad_qags(var, a, b, epsrel=desired_relative_error,
                           limit=maximum_num_subintervals)
        return v[0], v[1], Integer(v[2]), Integer(v[3])

    def integral(self, var='x', min=None, max=None):
        r"""
        Return the integral of ``self`` with respect to the variable `x`.

        INPUT:

        - ``var`` -- variable

        - ``min`` -- (default: ``None``)

        - ``max`` -- (default: ``None``)

        OUTPUT: the definite integral if xmin is not ``None``

        - an indefinite integral otherwise

        EXAMPLES::

            sage: maxima('x^2+1').integral()
            x^3/3+x
            sage: maxima('x^2+ 1 + y^2').integral('y')
            y^3/3+x^2*y+y
            sage: maxima('x / (x^2+1)').integral()
            log(x^2+1)/2
            sage: maxima('1/(x^2+1)').integral()
            atan(x)
            sage: maxima('1/(x^2+1)').integral('x', 0, infinity)
            %pi/2
            sage: maxima('x/(x^2+1)').integral('x', -1, 1)
            0

        ::

            sage: f = maxima('exp(x^2)').integral('x',0,1)
            sage: f.sage()
            -1/2*I*sqrt(pi)*erf(I)
            sage: f.numer()
            1.46265174590718...
        """
        I = InterfaceElement.__getattr__(self, 'integrate')
        if min is None:
            return I(var)
        else:
            if max is None:
                raise ValueError("neither or both of min/max must be specified.")
            return I(var, min, max)

    integrate = integral

    def __float__(self):
        """
        Return floating point version of this Maxima element.

        OUTPUT: real

        EXAMPLES::

            sage: float(maxima("3.14"))
            3.14
            sage: float(maxima("1.7e+17"))
            1.7e+17
            sage: float(maxima("1.7e-17"))
            1.7e-17
        """
        try:
            return float(repr(self.numer()))
        except ValueError:
            raise TypeError("unable to coerce '%s' to float" % repr(self))

    def __len__(self):
        """
        Return the length of a list.

        OUTPUT: integer

        EXAMPLES::

            sage: v = maxima('create_list(x^i,i,0,5)')
            sage: len(v)
            6
        """
        P = self._check_valid()
        return int(P.eval('length(%s)' % self.name()))

    def dot(self, other):
        """
        Implement the notation ``self . other``.

        INPUT:

        - ``other`` -- matrix; argument to dot

        OUTPUT: Maxima matrix

        EXAMPLES::

            sage: A = maxima('matrix ([a1],[a2])')
            sage: B = maxima('matrix ([b1, b2])')
            sage: A.dot(B)
            matrix([a1*b1,a1*b2],[a2*b1,a2*b2])
        """
        P = self._check_valid()
        Q = P(other)
        return P('%s . %s' % (self.name(), Q.name()))

    def __getitem__(self, n):
        r"""
        Return the `n`-th element of this list.

        INPUT:

        - ``n`` -- integer

        OUTPUT: Maxima object

        .. NOTE::

           Lists are 0-based when accessed via the Sage interface, not
           1-based as they are in the Maxima interpreter.

        EXAMPLES::

            sage: v = maxima('create_list(i*x^i,i,0,5)'); v
            [0,x,2*x^2,3*x^3,4*x^4,5*x^5]
            sage: v[3]
            3*x^3
            sage: v[0]
            0
            sage: v[10]
            Traceback (most recent call last):
            ...
            IndexError: n = (10) must be between 0 and 5
        """
        n = int(n)
        if n < 0 or n >= len(self):
            raise IndexError("n = (%s) must be between %s and %s" % (n, 0, len(self)-1))
        # If you change the n+1 to n below, better change __iter__ as well.
        return InterfaceElement.__getitem__(self, n+1)

    def __iter__(self):
        """
        Return an iterator for ``self``.

        OUTPUT: iterator

        EXAMPLES::

            sage: v = maxima('create_list(i*x^i,i,0,5)')
            sage: L = list(v)
            sage: [e._sage_() for e in L]
            [0, x, 2*x^2, 3*x^3, 4*x^4, 5*x^5]
        """
        for i in range(len(self)):
            yield self[i]

    def subst(self, val):
        """
        Substitute a value or several values into this Maxima object.

        INPUT:

        - ``val`` -- string representing substitution(s) to perform

        OUTPUT: Maxima object

        EXAMPLES::

            sage: maxima('a^2 + 3*a + b').subst('b=2')
            a^2+3*a+2
            sage: maxima('a^2 + 3*a + b').subst('a=17')
            b+340
            sage: maxima('a^2 + 3*a + b').subst('a=17, b=2')
            342
        """
        return self.comma(val)

    def comma(self, args):
        """
        Form the expression that would be written 'self, args' in Maxima.

        INPUT:

        - ``args`` -- string

        OUTPUT: Maxima object

        EXAMPLES::

            sage: maxima('sqrt(2) + I').comma('numer')
            I+1.41421356237309...
            sage: maxima('sqrt(2) + I*a').comma('a=5')
            5*I+sqrt(2)
        """
        self._check_valid()
        P = self.parent()
        return P('%s, %s' % (self.name(), args))

    def _latex_(self):
        r"""
        Return Latex representation of this Maxima object.

        OUTPUT: string

        This calls the tex command in Maxima, then does a little
        post-processing to fix bugs in the resulting Maxima output.

        EXAMPLES::

            sage: maxima('sqrt(2) + 1/3 + asin(5)')._latex_()
            '\\sin^{-1}\\cdot5+\\sqrt{2}+{{1}\\over{3}}'

            sage: y,d = var('y,d')
            sage: f = function('f')
            sage: latex(maxima(derivative(f(x*y), x)))
            \left(\left.{{{\it \partial}}\over{{\it \partial}\,  {\it \_symbol}_{0}}}\,f\left(  {\it \_symbol}_{0}\right)\right|_{  {\it \_symbol}_{0}={\it x}\,  {\it y}}\right)\,{\it y}
            sage: latex(maxima(derivative(f(x,y,d), d,x,x,y)))
            {{{\it \partial}^4}\over{{\it \partial}\,{\it d}\,  {\it \partial}\,{\it x}^2\,{\it \partial}\,  {\it y}}}\,f\left({\it x} ,  {\it y} , {\it d}\right)
            sage: latex(maxima(d/(d-2)))
            {{{\it d}}\over{{\it d}-2}}
        """
        self._check_valid()
        P = self.parent()
        s = P._eval_line('tex(%s);' % self.name(), reformat=False)
        if '$$' not in s:
            raise RuntimeError("Error texing Maxima object.")
        i = s.find('$$')
        j = s.rfind('$$')
        s = s[i+2:j]
        s = multiple_replace({'\r\n': ' ',
                              '\\%': '',
                              '\\arcsin ': '\\sin^{-1} ',
                              '\\arccos ': '\\cos^{-1} ',
                              '\\arctan ': '\\tan^{-1} ',
                              '\\_SAGE\\_VAR\\_': ''}, s)

        # Fix a maxima bug, which gives a latex representation of multiplying
        # two numbers as a single space. This was really bad when 2*17^(1/3)
        # gets TeXed as '2 17^{\frac{1}{3}}'
        #
        # This regex matches a string of spaces preceded by either a '}', a
        # decimal digit, or a ')', and followed by a decimal digit. The spaces
        # get replaced by a '\cdot'.
        return re.sub(r'(?<=[})\d]) +(?=\d)', r'\\cdot', s)

    def _tab_completion(self, verbose=False):
        """
        Return all Maxima commands, which is useful for tab completion.

        INPUT:

        - ``verbose`` -- boolean

        OUTPUT: list of strings

        EXAMPLES::

            sage: m = maxima(2)
            sage: 'gcd' in m._tab_completion()
            True
        """
        return self.parent()._tab_completion(verbose=False)

    def _matrix_(self, R):
        r"""
        If ``self`` is a Maxima matrix, return the corresponding Sage matrix
        over the Sage ring `R`.

        INPUT:

        - ``R`` -- ring to coerce into

        OUTPUT: matrix

        This may or may not work depending in how complicated the entries
        of ``self`` are! It only works if the entries of ``self`` can be coerced as
        strings to produce meaningful elements of `R`.

        EXAMPLES::

            sage: _ = maxima.eval("f[i,j] := i/j")
            sage: A = maxima('genmatrix(f,4,4)'); A
            matrix([1,1/2,1/3,1/4],[2,1,2/3,1/2],[3,3/2,1,3/4],[4,2,4/3,1])
            sage: A._matrix_(QQ)
            [  1 1/2 1/3 1/4]
            [  2   1 2/3 1/2]
            [  3 3/2   1 3/4]
            [  4   2 4/3   1]

        You can also use the ``matrix`` command (which is
        defined in ``sage.misc.functional``)::

            sage: matrix(QQ, A)
            [  1 1/2 1/3 1/4]
            [  2   1 2/3 1/2]
            [  3 3/2   1 3/4]
            [  4   2 4/3   1]
        """
        from sage.matrix.matrix_space import MatrixSpace
        self._check_valid()
        P = self.parent()
        nrows = int(P.eval('length(%s)' % self.name()))
        if nrows == 0:
            return MatrixSpace(R, 0, 0)(0)
        ncols = int(P.eval('length(%s[1])' % self.name()))
        M = MatrixSpace(R, nrows, ncols)
        s = self.str().replace('matrix', '').replace(',', "','").\
            replace("]','[", "','").replace('([', "['").replace('])', "']")
        s = eval(s)
        return M([R(x) for x in s])

    def partial_fraction_decomposition(self, var='x'):
        """
        Return the partial fraction decomposition of ``self`` with respect to
        the variable var.

        INPUT:

        - ``var`` -- string

        OUTPUT: Maxima object

        EXAMPLES::

            sage: f = maxima('1/((1+x)*(x-1))')
            sage: f.partial_fraction_decomposition('x')
            1/(2*(x-1))-1/(2*(x+1))
            sage: print(f.partial_fraction_decomposition('x'))
                                 1           1
                             --------- - ---------
                             2 (x - 1)   2 (x + 1)
        """
        return self.partfrac(var)

    def _operation(self, operation, other=None):
        r"""
        Return the result of applying the binary operation
        ``operation`` on the arguments ``self`` and ``other``, or the
        unary operation on ``self`` if ``other`` is not given.

        This is a utility function which factors out much of the
        commonality used in the arithmetic operations for interface
        elements.

        INPUT:

        - ``operation`` -- string representing the operation
          being performed. For example, '*', or '1/'

        - ``other`` -- the other operand. If ``other`` is ``None``,
          then the operation is assumed to be unary rather than binary

        OUTPUT: Maxima object

        Note that other's parent should already be Maxima since this should
        be called after coercion has been performed.

        If other is a ``MaximaFunction``, then we convert
        ``self`` to a ``MaximaFunction`` that takes
        no arguments, and let the
        ``MaximaFunction._operation`` code handle everything
        from there.

        EXAMPLES::

            sage: f = maxima.cos(x)
            sage: f._operation("+", f)
            2*cos(_SAGE_VAR_x)
        """
        P = self._check_valid()

        if other is None:
            cmd = '%s %s' % (operation, self._name)
        elif isinstance(other, P._object_function_class()):
            fself = P.function('', repr(self))
            return fself._operation(operation, other)
        else:
            cmd = '%s %s %s' % (self._name, operation, other._name)
        try:
            return P.new(cmd)
        except Exception as msg:
            raise TypeError(msg)


MaximaAbstractFunctionElement = InterfaceFunctionElement
MaximaAbstractFunction = InterfaceFunction


class MaximaAbstractElementFunction(MaximaAbstractElement):
    r"""
    Create a Maxima function with the parent ``parent``,
    name ``name``, definition ``defn``, arguments ``args``
    and latex representation ``latex``.

    INPUT:

    - ``parent`` -- an instance of a concrete Maxima interface

    - ``name`` -- string

    - ``defn`` -- string

    - ``args`` -- string; comma separated names of arguments

    - ``latex`` -- string

    OUTPUT: Maxima function

    EXAMPLES::

        sage: maxima.function('x,y','e^cos(x)')
        e^cos(x)
    """

    def __init__(self, parent, name, defn, args, latex):
        """
        Create a Maxima function.
        See ``MaximaAbstractElementFunction`` for full documentation.

        TESTS::

            sage: from sage.interfaces.maxima_abstract import MaximaAbstractElementFunction
            sage: MaximaAbstractElementFunction == loads(dumps(MaximaAbstractElementFunction))
            True
            sage: f = maxima.function('x,y','sin(x+y)')
            sage: f == loads(dumps(f))
            True
        """
        MaximaAbstractElement.__init__(self, parent, name, is_name=True)
        self.__defn = defn
        self.__args = args
        self.__latex = latex

    def __reduce__(self):
        """
        Implement __reduce__ for ``MaximaAbstractElementFunction``.

        OUTPUT: a couple consisting of:

        - the function to call for unpickling

        - a tuple of arguments for the function

        EXAMPLES::

            sage: f = maxima.function('x,y','sin(x+y)')
            sage: f.__reduce__()
            (<function reduce_load_MaximaAbstract_function at 0x...>,
             (Maxima, 'sin(x+y)', 'x,y', None))
        """
        return reduce_load_MaximaAbstract_function, (self.parent(),
                            self.__defn, self.__args, self.__latex)

    def __call__(self, *args):
        """
        Return the result of calling this Maxima function
        with the given arguments.

        INPUT:

        - ``args`` -- a variable number of arguments

        OUTPUT: Maxima object

        EXAMPLES::

            sage: f = maxima.function('x,y','sin(x+y)')
            sage: f(1,2)
            sin(3)
            sage: f(x,x)
            sin(2*x)
        """
        P = self._check_valid()
        if len(args) == 1:
            args = '(%s)' % args
        return P('%s%s' % (self.name(), args))

    def _repr_(self):
        """
        Return print representation of this Maxima function.

        OUTPUT: string

        EXAMPLES::

            sage: f = maxima.function('x,y','sin(x+y)')
            sage: repr(f)    # indirect doctest
            'sin(x+y)'
        """
        return self.definition()

    def _latex_(self):
        r"""
        Return latex representation of this Maxima function.

        OUTPUT: string

        EXAMPLES::

            sage: f = maxima.function('x,y','sin(x+y)')
            sage: latex(f)
            \mathrm{sin(x+y)}
        """
        if self.__latex is None:
            return r'\mathrm{%s}' % self.__defn
        else:
            return self.__latex

    def arguments(self, split=True):
        r"""
        Return the arguments of this Maxima function.

        INPUT:

        - ``split`` -- boolean; if ``True`` return a tuple of strings,
          otherwise return a string of comma-separated arguments

        OUTPUT: string if ``split`` is False

        - a list of strings if ``split`` is True

        EXAMPLES::

            sage: f = maxima.function('x,y','sin(x+y)')
            sage: f.arguments()
            ['x', 'y']
            sage: f.arguments(split=False)
            'x,y'
            sage: f = maxima.function('', 'sin(x)')
            sage: f.arguments()
            []
        """
        if split:
            return self.__args.split(',') if self.__args != '' else []
        else:
            return self.__args

    def definition(self):
        """
        Return the definition of this Maxima function as a string.

        EXAMPLES::

            sage: f = maxima.function('x,y','sin(x+y)')
            sage: f.definition()
            'sin(x+y)'
        """
        return self.__defn

    def integral(self, var):
        """
        Return the integral of ``self`` with respect to the variable var.

        INPUT:

        - ``var`` -- a variable

        OUTPUT: Maxima function

        Note that integrate is an alias of integral.

        EXAMPLES::

            sage: x,y = var('x,y')
            sage: f = maxima.function('x','sin(x)')
            sage: f.integral(x)
            -cos(x)
            sage: f.integral(y)
            sin(x)*y
        """
        var = str(var)
        P = self._check_valid()
        f = P('integrate(%s(%s), %s)' % (self.name(),
                        self.arguments(split=False), var))

        args = self.arguments()
        if var not in args:
            args.append(var)
        return P.function(",".join(args), repr(f))

    integrate = integral

    def _operation(self, operation, f=None):
        r"""
        This is a utility function which factors out much of the
        commonality used in the arithmetic operations for
        ``MaximaAbstractElementFunction``.

        INPUT:

        - ``operation`` -- string representing the operation
          being performed. For example, '\*', or '1/'

        - ``f`` -- the other operand; if ``f`` is
          ``None``, then the operation is assumed to be unary
          rather than binary

        EXAMPLES::

            sage: f = maxima.function('x,y','sin(x+y)')
            sage: f._operation("+", f)
            2*sin(y+x)
            sage: f._operation("+", 2)
            sin(y+x)+2
            sage: f._operation('-')
            -sin(y+x)
            sage: f._operation('1/')
            1/sin(y+x)
        """
        P = self._check_valid()
        if isinstance(f, P._object_function_class()):
            tmp = sorted(set(self.arguments() + f.arguments()))
            args = ','.join(tmp)
            defn = "(%s)%s(%s)" % (self.definition(), operation, f.definition())
        elif f is None:
            args = self.arguments(split=False)
            defn = "%s(%s)" % (operation, self.definition())
        else:
            args = self.arguments(split=False)
            defn = "(%s)%s(%s)" % (self.definition(), operation, repr(f))

        return P.function(args, P.eval(defn))


def reduce_load_MaximaAbstract_function(parent, defn, args, latex):
    r"""
    Unpickle a Maxima function.

    EXAMPLES::

        sage: from sage.interfaces.maxima_abstract import reduce_load_MaximaAbstract_function
        sage: f = maxima.function('x,y','sin(x+y)')
        sage: _,args = f.__reduce__()
        sage: g = reduce_load_MaximaAbstract_function(*args)
        sage: g == f
        True
    """
    return parent.function(args, defn, defn, latex)


def maxima_version():
    """
    Return Maxima version.

    Currently this calls a new copy of Maxima.

    EXAMPLES::

        sage: from sage.interfaces.maxima_abstract import maxima_version
        sage: maxima_version()  # random
        '5.41.0'
    """
    with os.popen('{} --version'.format(MAXIMA)) as p:
        return p.read().split()[-1]


def maxima_console():
    """
    Spawn a new Maxima command-line session.

    EXAMPLES::

        sage: from sage.interfaces.maxima_abstract import maxima_console
        sage: maxima_console()                    # not tested
        Maxima 5.46.0 https://maxima.sourceforge.io
        ...
    """
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%maxima magics instead.')
    os.system('{}'.format(MAXIMA))
