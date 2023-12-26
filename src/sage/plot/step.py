# sage.doctest: needs sage.symbolic
"""
Step function plots
"""

# *****************************************************************************
#       Copyright (C) 2009 William Stein <wstein@gmail.com>,
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
# *****************************************************************************
from sage.plot.line import line


def plot_step_function(v, vertical_lines=True, **kwds):
    r"""
    Return the line graphics object that gives the plot of the step
    function `f` defined by the list `v` of pairs `(a,b)`.  Here if
    `(a,b)` is in `v`, then `f(a) = b`.  The user does not have to
    worry about sorting the input list `v`.

    INPUT:

    - ``v`` -- list of pairs (a,b)

    - ``vertical_lines`` -- bool (default: ``True``) if ``True``, draw
      vertical risers at each step of this step function.
      Technically these vertical lines are not part of the graph
      of this function, but they look very nice in the plot, so we
      include them by default

    EXAMPLES:

    We plot the prime counting function::

        sage: plot_step_function([(i, prime_pi(i)) for i in range(20)])
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(plot_step_function([(i, prime_pi(i)) for i in range(20)]))

    ::

        sage: plot_step_function([(i, sin(i)) for i in range(5, 20)])
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(plot_step_function([(i, sin(i)) for i in range(5, 20)]))

    We pass in many options and get something that looks like "Space Invaders"::

        sage: v = [(i, sin(i)) for i in range(5, 20)]
        sage: plot_step_function(v, vertical_lines=False, thickness=30,
        ....:                    rgbcolor='purple', axes=False)
        Graphics object consisting of 14 graphics primitives

    .. PLOT::

        v = [(i, sin(i)) for i in range(5, 20)]
        sphinx_plot(plot_step_function(v, vertical_lines=False, thickness=30, rgbcolor='purple', axes=False))
    """
    # make sorted copy of v (don't change in place, since that would be rude).
    v = sorted(v)
    if len(v) <= 1:
        return line([])  # empty line
    if vertical_lines:
        w = []
        for i in range(len(v)):
            w.append(v[i])
            if i+1 < len(v):
                w.append((v[i+1][0], v[i][1]))
        return line(w, **kwds)
    else:
        return sum(line([v[i], (v[i+1][0], v[i][1])], **kwds) for i in range(len(v)-1))
