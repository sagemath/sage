r"""
sine-Gordon Y-system plotter

This class builds the triangulations associated to sine-Gordon and reduced
sine-Gordon Y-systems as constructed in [NS]_.

AUTHORS:

- Salvatore Stella (2014-07-18): initial version

EXAMPLES:

A reduced sine-Gordon example with 3 generations::

    sage: Y = SineGordonYsystem('A',(6,4,3)); Y
    A sine-Gordon Y-system of type A with defining integer tuple (6, 4, 3)
    sage: Y.plot()     #not tested

The same integer tuple but for the non-reduced case::

    sage: Y = SineGordonYsystem('D',(6,4,3)); Y
    A sine-Gordon Y-system of type D with defining integer tuple (6, 4, 3)
    sage: Y.plot()     #not tested

.. TODO::

    The code for plotting is extremely slow.

REFERENCES:

.. [NS] \T. Nakanishi, S. Stella, Wonder of sine-Gordon Y-systems,
   to appear in Trans. Amer. Math. Soc., :arxiv:`1212.6853`
"""
#*****************************************************************************
#       Copyright (C) 2014 Salvatore Stella <sstella@ncsu.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject

from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RR
from sage.rings.all import NN
from sage.functions.trig import cos, sin
from sage.plot.plot import parametric_plot
from sage.plot.graphics import Graphics
from sage.plot.polygon import polygon2d
from sage.plot.circle import circle
from sage.plot.bezier_path import bezier_path
from sage.plot.point import point
from sage.plot.line import line
from sage.symbolic.constants import pi, I
from sage.functions.log import exp
from sage.functions.other import ceil
from sage.misc.flatten import flatten
from sage.symbolic.ring import SR
from sage.functions.other import real_part, imag_part
from sage.misc.cachefunc import cached_method


class SineGordonYsystem(SageObject):
    r"""
    A class to model a (reduced) sine-Gordon Y-system

    Note that the generations, together with all integer tuples, in this
    implementation are numbered from 0 while in [NS]_ they are numbered from 1

    INPUT:

    - ``X`` -- the type of the Y-system to construct (either 'A' or 'D')
    - ``na`` -- the tuple of positive integers defining the Y-system
      with ``na[0] > 2``

    See [NS]_

    EXAMPLES::

        sage: Y = SineGordonYsystem('A',(6,4,3)); Y
        A sine-Gordon Y-system of type A with defining integer tuple (6, 4, 3)
        sage: Y.intervals()
        (((0, 0, 'R'),),
         ((0, 17, 'L'),
          (17, 34, 'L'),
        ...
          (104, 105, 'R'),
          (105, 0, 'R')))
        sage: Y.triangulation()
        ((17, 89),
         (17, 72),
         (34, 72),
        ...
         (102, 105),
         (103, 105))
        sage: Y.plot()     #not tested
    """
    def __init__(self, X, na):
        """
        TESTS::

            sage: Y = SineGordonYsystem('A',(6,4,3)); Y  # indirect doctest
            A sine-Gordon Y-system of type A with defining integer tuple
            (6, 4, 3)

            sage: SineGordonYsystem('E',(6,4,3))
            Traceback (most recent call last):
            ...
            ValueError: the type must be either 'A' of 'D'.
            sage: SineGordonYsystem('A',(2,4,3))
            Traceback (most recent call last):
            ...
            ValueError: the first integer in the defining sequence must be
            greater than 2.
            sage: SineGordonYsystem('A',(6,-4,3))
            Traceback (most recent call last):
            ...
            ValueError: the defining sequence must contain only positive
            integers.
            sage: SineGordonYsystem('A',(3,))
            Traceback (most recent call last):
            ...
            ValueError: the integer sequence (3,) in type 'A' is not allowed
            as input
        """
        if X not in ['A', 'D']:
            raise ValueError("the type must be either 'A' of 'D'.")
        self._type = X
        if na[0] <= 2:
            raise ValueError("the first integer in the defining sequence "
                             "must be greater than 2.")
        if any(x not in NN for x in na):
            raise ValueError("the defining sequence must contain only "
                             "positive integers.")
        self._na = tuple(na)
        if self._na == (3,) and self._type == 'A':
            raise ValueError("the integer sequence (3,) in type 'A'"
                             " is not allowed as input")
        self._F = len(self._na)

    def _repr_(self):
        """
        Return the string representation of ``self``.

        TESTS::

            sage: Y = SineGordonYsystem('A',(6,4,3)); Y  # indirect doctest
            A sine-Gordon Y-system of type A with defining integer tuple
            (6, 4, 3)
            sage: Y = SineGordonYsystem('D',(6,4,3)); Y  # indirect doctest
            A sine-Gordon Y-system of type D with defining integer tuple
            (6, 4, 3)
        """
        msg = "A sine-Gordon Y-system of type {}"
        msg += " with defining integer tuple {}"
        return  msg.format(self._type, self._na)

    def type(self):
        r"""
        Return the type of ``self``.

        EXAMPLES::

            sage: Y = SineGordonYsystem('A',(6,4,3))
            sage: Y.type()
            'A'
        """
        return self._type

    def F(self):
        r"""
        Return the number of generations in ``self``.

        EXAMPLES::

            sage: Y = SineGordonYsystem('A',(6,4,3))
            sage: Y.F()
            3
        """
        return self._F

    def na(self):
        r"""
        Return the sequence of the integers `n_a` defining ``self``.

        EXAMPLES::

            sage: Y = SineGordonYsystem('A',(6,4,3))
            sage: Y.na()
            (6, 4, 3)
        """
        return self._na

    @cached_method
    def rk(self):
        r"""
        Return the sequence of integers ``r^{(k)}``, i.e. the width of
        an interval of type 'L' or 'R' in the ``k``-th generation.

        EXAMPLES::

            sage: Y = SineGordonYsystem('A',(6,4,3))
            sage: Y.rk()
            (106, 17, 4)
        """
        na = self._na
        F = self._F
        rk = [na[F - 1] + 1]
        if F > 1:
            rk.append(na[F - 2] * na[F - 1] + na[F - 2] + 1)
            for k in range(2, F):
                rk.append(na[F - k - 1] * rk[k - 1] + rk[k - 2])
        rk.reverse()
        return tuple(rk)

    @cached_method
    def pa(self):
        r"""
        Return the sequence of integers  ``p_a``, i.e. the total number of
        intervals of types 'NL' and 'NR' in the ``(a+1)``-th generation.

        EXAMPLES::

            sage: Y = SineGordonYsystem('A',(6,4,3))
            sage: Y.pa()
            (1, 6, 25)
        """
        na = self._na
        F = self._F
        pa = [1]
        if F > 1:
            pa.append(na[0])
            for k in range(2, F):
                pa.append(na[k-1] * pa[k-1] + pa[k-2])
        return tuple(pa)

    @cached_method
    def qa(self):
        r"""
        Return the sequence of integers  ``q_a``, i.e. the total number of
        intervals of types 'L' and 'R' in the ``(a+1)``-th generation.

        EXAMPLES::

            sage: Y = SineGordonYsystem('A',(6,4,3))
            sage: Y.qa()
            (6, 25, 81)
        """
        na = self._na
        F = self._F
        qa = [na[0]]
        if F > 1:
            qa.append(na[1] * qa[0] + 1)
            for k in range(2, F):
                qa.append(na[k] * qa[k - 1] + qa[k - 2])
        return tuple(qa)

    @cached_method
    def r(self):
        r"""
        Return the number of vertices in the polygon realizing ``self``.

        EXAMPLES::

            sage: Y = SineGordonYsystem('A',(6,4,3))
            sage: Y.r()
            106
        """
        return self.rk()[0]

    @cached_method
    def vertices(self):
        r"""
        Return the vertices of the polygon realizing ``self`` as the ring of
        integers modulo ``self.r()``.

        EXAMPLES::

            sage: Y = SineGordonYsystem('A',(6,4,3))
            sage: Y.vertices()
            Ring of integers modulo 106
        """
        return ZZ.quotient(self.r())

    @cached_method
    def triangulation(self):
        r"""
        Return the initial triangulation of the polygon realizing
        ``self`` as a tuple of pairs of vertices.

        .. WARNING::

            In type 'D' the returned triangulation does NOT contain the two
            radii.

        ALGORITHM:

        We implement the four cases described by Figure 14 in [NS]_.

        EXAMPLES::

            sage: Y = SineGordonYsystem('A',(6,4,3))
            sage: Y.triangulation()
            ((17, 89),
             (17, 72),
            ...
             (102, 105),
             (103, 105))
        """
        rk = self.rk() + (1, 1)
        na = self.na()
        vert = self.vertices()
        triangulation = []
        intervals = self.intervals()
        for a in range(self.F()):
            for (first, last, typ) in intervals[a]:
                if first - last in [vert(1), vert(-1)]:
                    continue
                if typ == "L":
                    left = True
                    if na[a] % 2 == 0:
                        last_cw = first + vert(na[a] / 2 * rk[a + 1])
                        last_ccw = last - vert(na[a] / 2 * rk[a + 1])
                    else:
                        last_cw = first + vert((na[a] + 1) / 2 * rk[a + 1])
                        last_ccw = last - vert((na[a] - 1) / 2 * rk[a + 1])
                elif typ == "R":
                    left = False
                    if na[a] % 2 == 0:
                        last_cw = first + vert(na[a] / 2 * rk[a + 1])
                        last_ccw = last - vert(na[a] / 2 * rk[a + 1])
                    else:
                        last_cw = first + vert((na[a] - 1) / 2 * rk[a + 1])
                        last_ccw = last - vert((na[a] + 1) / 2 * rk[a + 1])
                else:
                    continue
                if first == last:
                    # this happens only when the interval is the whole disk
                    first = first + vert(rk[a + 1])
                    last = last - vert(rk[a + 1])
                edge = (first, last)
                triangulation.append(edge)
                done = False
                while not done:
                    if left:
                        edge = (edge[0] + vert(rk[a+1]), edge[1])
                    else:
                        edge = (edge[0], edge[1] - vert(rk[a + 1]))
                    left = not left
                    if (edge[1] >= last_ccw and edge[0] < last_cw) or (edge[1] > last_ccw and edge[0] <= last_cw):
                        triangulation.append(edge)
                    else:
                        done = True
        if self.type() == 'D':
            triangulation.append((vert(0), vert(rk[0] - rk[1])))
        return tuple(triangulation)

    @cached_method
    def intervals(self):
        r"""
        Return, divided by generation, the list of intervals used to construct
        the initial triangulation.

        Each such interval is a triple ``(p, q, X)`` where ``p`` and
        ``q`` are the two extremal vertices of the interval and ``X``
        is the type of the interval (one of 'L', 'R', 'NL', 'NR').

        ALGORITHM:

        The algorithm used here is the one described in section 5.1 of [NS]_.
        The only difference is that we get rid of the special case of the first
        generation by treating the whole disk as a type 'R' interval.

        EXAMPLES::

            sage: Y = SineGordonYsystem('A',(6,4,3))
            sage: Y.intervals()
            (((0, 0, 'R'),),
             ((0, 17, 'L'),
              (17, 34, 'L'),
            ...
              (104, 105, 'R'),
              (105, 0, 'R')))
        """
        rk = self.rk() + (1, 1)
        na = self.na()
        vert = self.vertices()
        intervals = [[(vert(0), vert(0), "R")]]
        for a in range(self.F()):
            new_intervals = []
            if na[a] % 2 == 0:
                for (first, last, typ) in intervals[a]:
                    if typ == "NR":
                        new_intervals.append((first, last, "R"))
                    elif typ == "NL":
                        new_intervals.append((first, last, "L"))
                    else:
                        last_cw = first + vert(na[a] / 2 * rk[a + 1])
                        last_ccw = vert(last_cw + rk[a + 2])
                        x = first
                        while x < last_cw:
                            new_intervals.append((vert(x), vert(x + rk[a+1]), "L"))
                            x = vert(x + rk[a + 1])
                        if typ == "L":
                            new_intervals.append((last_cw, last_ccw, "NL"))
                        else:
                            new_intervals.append((last_cw, last_ccw, "NR"))
                        x = last_ccw
                        while x != last:
                            new_intervals.append((vert(x), vert(x+rk[a+1]), "R"))
                            x = vert(x + rk[a + 1])
            else:
                for (first, last, typ) in intervals[a]:
                    if typ == "NR":
                        new_intervals.append((first, last, "R"))
                    elif typ == "NL":
                        new_intervals.append((first, last, "L"))
                    else:
                        if typ == "L":
                            last_cw = first + vert((na[a] + 1) / 2 * rk[a + 1])
                        else:
                            last_cw = first + vert((na[a] - 1) / 2 * rk[a + 1])
                        last_ccw = vert(last_cw + rk[a + 2])
                        x = first
                        while x < last_cw:
                            new_intervals.append((vert(x), vert(x + rk[a+1]), "L"))
                            x = vert(x + rk[a+1])
                        if typ == "L":
                            new_intervals.append((last_cw, last_ccw, "NR"))
                        else:
                            new_intervals.append((last_cw, last_ccw, "NL"))
                        x = last_ccw
                        while x != last:
                            new_intervals.append((vert(x), vert(x + rk[a+1]), "R"))
                            x = vert(x + rk[a + 1])
            intervals.append(new_intervals)
        return tuple(map(tuple, intervals))

    def plot(self, **kwds):
        r"""
        Plot the initial triangulation associated to ``self``.

        INPUT:

        - ``radius`` - the radius of the disk; by default the length of
          the circle is the number of vertices
        - ``points_color`` - the color of the vertices; default 'black'
        - ``points_size`` - the size of the vertices; default 7
        - ``triangulation_color`` - the color of the arcs; default 'black'
        - ``triangulation_thickness`` - the thickness of the arcs; default 0.5
        - ``shading_color`` - the color of the shading used on neuter
          intervals; default 'lightgray'
        - ``reflections_color`` - the color of the reflection axes; default
          'blue'
        - ``reflections_thickness`` - the thickness of the reflection axes;
          default 1

        EXAMPLES::

            sage: Y = SineGordonYsystem('A',(6,4,3))
            sage: Y.plot()  # long time 2s
            Graphics object consisting of 219 graphics primitives
        """
        # Set up plotting options
        if 'radius' in kwds:
            radius = kwds['radius']
        else:
            radius = ceil(self.r() / (2 * pi))
        points_opts = {}
        if 'points_color' in kwds:
            points_opts['color'] = kwds['points_color']
        else:
            points_opts['color'] = 'black'
        if 'points_size' in kwds:
            points_opts['size'] = kwds['points_size']
        else:
            points_opts['size'] = 7
        triangulation_opts = {}
        if 'triangulation_color' in kwds:
            triangulation_opts['color'] = kwds['triangulation_color']
        else:
            triangulation_opts['color'] = 'black'
        if 'triangulation_thickness' in kwds:
            triangulation_opts['thickness'] = kwds['triangulation_thickness']
        else:
            triangulation_opts['thickness'] = 0.5
        shading_opts = {}
        if 'shading_color' in kwds:
            shading_opts['color'] = kwds['shading_color']
        else:
            shading_opts['color'] = 'lightgray'
        reflections_opts = {}
        if 'reflections_color' in kwds:
            reflections_opts['color'] = kwds['reflections_color']
        else:
            reflections_opts['color'] = 'blue'
        if 'reflections_thickness' in kwds:
            reflections_opts['thickness'] = kwds['reflections_thickness']
        else:
            reflections_opts['thickness'] = 1
        # Helper functions

        def triangle(x):
            (a, b) = sorted(x[:2])
            for p in self.vertices():
                if (p, a) in self.triangulation() or (a, p) in self.triangulation():
                    if (p, b) in self.triangulation() or (b, p) in self.triangulation():
                        if p < a or p > b:
                            return sorted((a, b, p))

        def plot_arc(radius, p, q, **opts):
            # TODO: THIS SHOULD USE THE EXISTING PLOT OF ARCS!
            # plot the arc from p to q differently depending on the type of self
            p = ZZ(p)
            q = ZZ(q)
            t = SR.var('t')
            if p - q in [1, -1]:
                def f(t):
                    return (radius * cos(t), radius * sin(t))
                (p, q) = sorted([p, q])
                angle_p = vertex_to_angle(p)
                angle_q = vertex_to_angle(q)
                return parametric_plot(f(t), (t, angle_q, angle_p), **opts)
            if self.type() == 'A':
                angle_p = vertex_to_angle(p)
                angle_q = vertex_to_angle(q)
                if angle_p < angle_q:
                    angle_p += 2 * pi
                internal_angle = angle_p - angle_q
                if internal_angle > pi:
                    (angle_p, angle_q) = (angle_q + 2 * pi, angle_p)
                    internal_angle = angle_p - angle_q
                angle_center = (angle_p+angle_q) / 2
                hypotenuse = radius / cos(internal_angle / 2)
                radius_arc = hypotenuse * sin(internal_angle / 2)
                center = (hypotenuse * cos(angle_center),
                          hypotenuse * sin(angle_center))
                center_angle_p = angle_p + pi / 2
                center_angle_q = angle_q + 3 * pi / 2

                def f(t):
                    return (radius_arc * cos(t) + center[0],
                            radius_arc * sin(t) + center[1])
                return parametric_plot(f(t), (t, center_angle_p,
                                              center_angle_q), **opts)
            elif self.type() == 'D':
                if p >= q:
                    q += self.r()
                px = -2 * pi * p / self.r() + pi / 2
                qx = -2 * pi * q / self.r() + pi / 2
                arc_radius = (px - qx) / 2
                arc_center = qx + arc_radius

                def f(t):
                    return exp(I * ((cos(t) + I * sin(t)) *
                                    arc_radius + arc_center)) * radius
                return parametric_plot((real_part(f(t)), imag_part(f(t))),
                                       (t, 0, pi), **opts)

        def vertex_to_angle(v):
            # v==0 corresponds to pi/2
            return -2 * pi * RR(v) / self.r() + 5 * pi / 2

        # Begin plotting
        P = Graphics()
        # Shade neuter intervals
        neuter_intervals = [x for x in flatten(self.intervals()[:-1],
                                               max_level=1)
                            if x[2] in ["NR", "NL"]]
        shaded_triangles = map(triangle, neuter_intervals)
        for (p, q, r) in shaded_triangles:
            points = list(plot_arc(radius, p, q)[0])
            points += list(plot_arc(radius, q, r)[0])
            points += list(reversed(plot_arc(radius, p, r)[0]))
            P += polygon2d(points, **shading_opts)
        # Disk boundary
        P += circle((0, 0), radius, **triangulation_opts)
        # Triangulation
        for (p, q) in self.triangulation():
            P += plot_arc(radius, p, q, **triangulation_opts)
        if self.type() == 'D':
            s = radius / 50.0
            P += polygon2d([(s, 5 * s), (s, 7 * s),
                            (3 * s, 5 * s), (3 * s, 7 * s)],
                           color=triangulation_opts['color'])
            P += bezier_path([[(0, 0), (2 * s, 1 * s), (2 * s, 6 * s)],
                              [(2 * s, 10 * s), (s, 20 * s)],
                              [(0, 30 * s), (0, radius)]],
                             **triangulation_opts)
            P += bezier_path([[(0, 0), (-2 * s, 1 * s), (-2 * s, 6 * s)],
                              [(-2 * s, 10 * s), (-s, 20 * s)],
                              [(0, 30 * s), (0, radius)]],
                             **triangulation_opts)
            P += point((0, 0), zorder=len(P), **points_opts)
        # Vertices
        v_points = {x: (radius * cos(vertex_to_angle(x)),
                        radius * sin(vertex_to_angle(x)))
                    for x in self.vertices()}
        for v in v_points:
            P += point(v_points[v], zorder=len(P), **points_opts)
        # Reflection axes
        P += line([(0, 1.1 * radius), (0, -1.1 * radius)],
                  zorder=len(P), **reflections_opts)
        axis_angle = vertex_to_angle(-0.5 * (self.rk() + (1, 1))[1])
        (a, b) = (1.1 * radius * cos(axis_angle),
                  1.1 * radius * sin(axis_angle))
        P += line([(a, b), (-a, -b)], zorder=len(P), **reflections_opts)
        # Wrap up
        P.set_aspect_ratio(1)
        P.axes(False)
        return P
