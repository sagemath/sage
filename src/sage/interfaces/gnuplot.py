# sage_setup: distribution = sagemath-plot
r"""
Interface to the Gnuplot interpreter
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
import time
from sage.structure.sage_object import SageObject


class Gnuplot(SageObject):
    """
    Interface to the Gnuplot interpreter.
    """
    def _quit_string(self):
        return 'quit'

    def gnuplot(self):
        try:
            return self._gnuplot
        except AttributeError:
            try:
                import Gnuplot as GP
                self._gnuplot = GP.Gnuplot()
                return self._gnuplot
            except ImportError:
                raise RuntimeError("Install the gnuplotpy Python module.")

    def __call__(self, line):
        return self.gnuplot()(line)

    def _eval_line(self, line, *args, **kwds):
        self(line)
        return ''

    def _repr_(self):
        """
        TESTS::

            sage: gnuplot               # indirect doctests
            Interface to Gnuplot
        """
        return "Interface to Gnuplot"

    def plot(self, cmd, file=None, verbose=True, reset=True):
        r"""
        Draw the plot described by cmd, and possibly also save to an eps or
        png file.

        INPUT:

        - ``cmd`` -- string

        - ``file`` -- string (default: ``None``); if specified save
          plot to given file, which may be either an eps (default) or png file

        - ``verbose`` -- print some info

        - ``reset`` -- true; reset gnuplot before making graph

        OUTPUT: displays graph

        .. NOTE::

           Note that ``^`` s  are replaced by ``**`` s before being passed to gnuplot.
        """
        if reset:
            self('reset')
        self('set terminal x11')
        cmd = cmd.replace('^', '**')
        self(cmd)
        if file is not None:
            if file[-4:] == '.png':
                self('set terminal png medium')
            else:
                if file[-4:] != '.eps':
                    file += '.eps'
                self('set terminal postscript eps enhanced')
            # self("set output '%s'"%file)
            tmp = 'gnuplot_tmp%s' % file[-4:]
            self("set output '%s'" % tmp)
            print("Saving plot to %s" % file)
            self(cmd)
            time.sleep(0.1)
            os.system('mv %s %s 2>/dev/null' % (tmp, file))
            time.sleep(0.1)
        self('set terminal x11')

    def plot3d(self, f, xmin=-1, xmax=1, ymin=-1, ymax=1, zmin=-1, zmax=1,
               title=None,
               samples=25, isosamples=20, xlabel='x', ylabel='y',
               interact=True):
        if title is None:
            title = str(f)
        f = f.replace('^', '**')
        cmd = """
        set xlabel "%s"
        set ylabel "%s"
        set key top
        set border 4095
        set xrange [%s:%s]
        set yrange [%s:%s]
        set samples %s
        set isosamples %s

        set title "%s"
        set pm3d; set palette
        #show pm3d
        #show palette
        splot %s
        """ % (xlabel, ylabel,
             xmin, xmax, ymin, ymax,  # zmin, zmax,
             samples, isosamples,
             title, f)
        if interact:
            self.interact(cmd)
        else:
            self(cmd)

    def plot3d_parametric(self, f='cos(u)*(3 + v*cos(u/2)), sin(u)*(3 + v*cos(u/2)), v*sin(u/2)',
                          range1='[u=-pi:pi]',
                          range2='[v=-0.2:0.2]', samples=50, title=None,
                          interact=True):
        r"""
        Draw a parametric 3d surface and rotate it interactively.

        INPUT:

        - ``f`` -- string; a function of two variables, e.g.,
          ``'cos(u)\*(3 + v\*cos(u/2)), sin(u)\*(3 + v\*cos(u/2)), v\*sin(u/2)'``

        - ``range1`` -- string; range of values for one
          variable, e.g., ``'[u=-pi:pi]'``

        - ``range2`` -- string; range of values for another
          variable, e.g., ``'[v=-0.2:0.2]'``

        - ``samples`` -- integer; number of sample points to use

        - ``title`` -- string; title of the graph

        EXAMPLES::

            sage: gnuplot.plot3d_parametric('v^2*sin(u), v*cos(u), v*(1-v)')   # optional - gnuplot, not tested (since something pops up)
        """
        if title is None:
            title = str(f)
        cmd = """
        set key top
        set border 4095
        set samples %s

        set title "%s"
        set pm3d; set palette; set parametric
        splot %s %s %s
        """ % (samples, title, range1, range2, f)
        cmd = cmd.replace('^', '**')
        if interact:
            self.interact(cmd)
        else:
            self(cmd)

    def interact(self, cmd):
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w+t') as f:
            f.write(cmd + '\n pause -1 "Press return to continue (no further rotation possible)"')
            os.system(f'gnuplot -persist {f.name}')

    def console(self):
        gnuplot_console()


# An instance
gnuplot = Gnuplot()


def gnuplot_console():
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%gnuplot magics instead.')
    os.system('gnuplot')
