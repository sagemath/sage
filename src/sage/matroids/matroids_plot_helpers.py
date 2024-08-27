# sage.doctest: optional - scipy
r"""
Helper functions for plotting the geometric representation of matroids


AUTHORS:

- Jayant Apte (2014-06-06): initial version

.. NOTE::

    This file provides functions that are called by ``show()`` and ``plot()``
    methods of abstract matroids class. The basic idea is to first decide
    the placement of points in $\mathbb{R}^2$ and then draw lines in geometric
    representation through these points. Point placement procedures such as
    ``addtripts``, ``addnontripts`` together produce ``(x,y)`` tuples
    corresponding to ground set of the matroid in a dictionary.
    These methods provide simple but rigid point placement algorithm.
    Alternatively, one can build the point placement dictionary manually or
    via an optimization that gives aesthetically pleasing point placement (in
    some sense. This is not yet implemented).    One can then use
    ``createline`` function to produce sequence of ``100`` points on a smooth
    curve containing the points in the specified line which in turn uses
    :func:`scipy:scipy.interpolate.splprep` and
    :func:`scipy:scipy.interpolate.splev`.  Then one
    can use sage's graphics primitives ``line``, ``point``, ``text`` and
    ``points`` to produce graphics object containing points (ground set
    elements) and lines (for a rank 3 matroid, these are flats of rank 2 of
    size greater than equal to 3) of the geometric representation of the
    matroid. Loops and parallel elements are added as per conventions in
    [Oxl2011]_ using function ``addlp``. The priority order for point placement
    methods used inside plot() and show() is as follows:

    1. User Specified points dictionary and lineorders
    2. cached point placement dictionary and line orders (a list of ordered
       lists) in M._cached_info (a dictionary)
    3. Internal point placement and orders deciding heuristics
       If a custom point placement and/or line orders is desired, then user
       can simply specify the custom points dictionary as::

            M.cached info = {'plot_positions':<dictionary_of _points>,
                             'plot_lineorders':<list of lists>}


REFERENCES
==========

- [Oxl2011]_ James Oxley, "Matroid Theory, Second Edition". Oxford University
  Press, 2011.

EXAMPLES::

    sage: from sage.matroids import matroids_plot_helpers
    sage: M1 = Matroid(ring=GF(2), matrix=[[1, 0, 0, 0, 1, 1, 1,0,1,0,1],
    ....:                                  [0, 1, 0, 1, 0, 1, 1,0,0,1,0],
    ....:                                  [0, 0, 1, 1, 1, 0, 1,0,0,0,0]])
    sage: pos_dict = {0: (0, 0),  1: (2, 0),  2: (1, 2),  3: (1.5, 1.0),
    ....:             4: (0.5, 1.0),  5: (1.0, 0.0), 6: (1.0, 0.666666666666667),
    ....:             7: (3,3), 8: (4,0), 9: (-1,1), 10: (-2,-2)}
    sage: M1._cached_info = {'plot_positions': pos_dict, 'plot_lineorders': None}
    sage: matroids_plot_helpers.geomrep(M1, sp=True)                                    # needs sage.plot sage.rings.finite_rings
    Graphics object consisting of 22 graphics primitives
"""
# *****************************************************************************
#       Copyright (C) 2013 Jayant Apte <jayant91089@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

import scipy
import scipy.interpolate
import numpy as np
from sage.misc.lazy_import import lazy_import
lazy_import("sage.plot.all", ["Graphics", "line", "text", "polygon2d", "point", "points"])
lazy_import("sage.plot.colors", "Color")
from sage.sets.set import Set
from sage.matroids.advanced import newlabel


def it(M, B1, nB1, lps):
    r"""
    Return points on and off the triangle and lines to be drawn for a rank 3
    matroid.

    INPUT:

    - ``M`` -- matroid
    - ``B1`` -- list of groundset elements of ``M`` that corresponds to a
      basis of matroid ``M``
    - ``nB1`` -- list of elements in the ground set of M that corresponds to
      ``M.simplify.groundset() \ B1``
    - ``lps`` -- list of elements in the ground set of matroid M that are loops

    OUTPUT: a tuple containing 4 elements in this order:

    1. A dictionary containing 2-tuple (x,y) coordinates with
       ``M.simplify.groundset()`` elements that can be placed on the sides of
       the triangle as keys.
    2. A list of 3 lists of elements of ``M.simplify.groundset()`` that can
       be placed on the 3 sides of the triangle.
    3. A list of elements of ``M.simplify.groundset()`` that cane be placed
       inside the triangle in the geometric representation.
    4. A list of lists of elements of ``M.simplify.groundset()`` that
       correspond to lines in the geometric representation other than the sides
       of the triangle.

    EXAMPLES::

        sage: from sage.matroids import matroids_plot_helpers as mph
        sage: M = Matroid(ring=GF(2), matrix=[[1, 0, 0, 0, 1, 1, 1,0],
        ....: [0, 1, 0, 1, 0, 1, 1,0],[0, 0, 1, 1, 1, 0, 1,0]])
        sage: N = M.simplify()
        sage: B1 = list(N.basis())
        sage: nB1 = list(set(M.simplify().groundset())-set(B1))
        sage: pts,trilines,nontripts,curvedlines = mph.it(M,
        ....: B1,nB1,M.loops())
        sage: pts
        {1: (1.0, 0.0), 2: (1.5, 1.0), 3: (0.5, 1.0),
         4: (0, 0), 5: (1, 2), 6: (2, 0)}
        sage: trilines
        [[3, 4, 5], [2, 5, 6], [1, 4, 6]]
        sage: nontripts
        [0]
        sage: curvedlines
        [[0, 1, 5], [0, 2, 4], [0, 3, 6], [1, 2, 3], [1, 4, 6], [2, 5, 6],
         [3, 4, 5]]

    .. NOTE::

        This method does NOT do any checks.
    """
    tripts = [(0, 0), (1, 2), (2, 0)]
    pts = {}
    j = 0
    for i in B1:
        pts[i] = tripts[j]
        j = j + 1
    pairs = [[0, 1], [1, 2], [0, 2]]
    L1 = []
    L2 = []
    L3 = []
    for i in nB1:
        if M.is_dependent([i, B1[pairs[0][0]], B1[pairs[0][1]]]):
            # Add to L1
            L1.append(i)
        elif M.is_dependent([i, B1[pairs[1][0]], B1[pairs[1][1]]]):
            # Add to L2
            L2.append(i)
        elif M.is_dependent([i, B1[pairs[2][0]], B1[pairs[2][1]]]):
            # Add to L3
            L3.append(i)
    L = [L1, L2, L3]  # megalist
    lines = []  # the list of lines
    for i in range(1, len(L)+1):
        lines.append([B1[pairs[i-1][0]]])
        lines[i-1].extend(L[i-1])
        lines[i-1].extend([B1[pairs[i-1][1]]])
    # place triangle and L1,L2,L3
    for i in L:  # loop over megalist
        interval = 1/float(len(i)+1)
        pt1 = list(tripts[pairs[L.index(i)][0]])
        pt2 = list(tripts[pairs[L.index(i)][1]])
        for j in range(1, len(i)+1):
            # loop over L1,L2,L3
            cc = interval*j
            pts[i[j-1]] = (cc*pt1[0]+(1-cc)*pt2[0], cc*pt1[1]+(1-cc)*pt2[1])
    trilines = [set(x) for x in lines if len(x) >= 3]
    set_lps = set(lps)
    curvedlines = [list(sx.difference(set_lps))
                   for x in M.flats(2) if (sx := set(x)) not in trilines
                   and len(list(x)) >= 3]
    nontripts = [i for i in nB1 if i not in pts]
    trilines = [list(s) for s in trilines]
    return pts, trilines, nontripts, curvedlines


def trigrid(tripts):
    """
    Return a grid of 4 points inside given 3 points as a list.

    INPUT:

    - ``tripts`` -- list of 3 lists of the form [x,y] where x and y are the
      Cartesian coordinates of a point

    OUTPUT: list of lists containing 4 points in following order:

    - 1. Barycenter of 3 input points.
    - 2,3,4. Barycenters of 1. with 3 different 2-subsets of input points
      respectively.

    EXAMPLES::

        sage: from sage.matroids import matroids_plot_helpers
        sage: points = matroids_plot_helpers.trigrid([[2,1],[4,5],[5,2]])
        sage: points
        [[3.6666666666666665, 2.6666666666666665],
         [3.222222222222222, 2.888888888888889],
         [4.222222222222222, 3.222222222222222],
         [3.5555555555555554, 1.8888888888888886]]

    .. NOTE::

            This method does NOT do any checks.
    """
    pairs = [[0, 1], [1, 2], [0, 2]]
    cpt = [float(tripts[0][0] + tripts[1][0] + tripts[2][0]) / 3,
           float(tripts[0][1] + tripts[1][1] + tripts[2][1]) / 3]
    grid = [cpt]
    for p, q in pairs:
        pt = [float(tripts[p][0] + tripts[q][0] + cpt[0]) / 3,
              float(tripts[p][1] + tripts[q][1] + cpt[1]) / 3]
        grid.append(pt)
    return grid


def addnontripts(tripts_labels, nontripts_labels, ptsdict):
    """
    Return modified ``ptsdict`` with additional keys and values corresponding
    to ``nontripts``.

    INPUT:

    - ``tripts`` -- list of 3 ground set elements that are to be placed on
      vertices of the triangle
    - ``ptsdict`` -- dictionary (at least) containing ground set elements in
      ``tripts`` as keys and their (x,y) position as values
    - ``nontripts`` -- list of ground set elements whose corresponding points
      are to be placed inside the triangle

    OUTPUT:

    A dictionary containing ground set elements in ``tripts`` as keys and
    their (x,y) position as values along with all keys and respective values
    in ``ptsdict``.

    EXAMPLES::

        sage: from sage.matroids import matroids_plot_helpers
        sage: from sage.matroids.advanced import setprint
        sage: ptsdict={'a':(0,0),'b':(1,2),'c':(2,0)}
        sage: ptsdict_1=matroids_plot_helpers.addnontripts(['a','b','c'],
        ....:         ['d','e','f'],ptsdict)
        sage: setprint(ptsdict_1)
        {'a': [0, 0], 'b': [1, 2], 'c': [0, 2], 'd': [0.6666666666666666, 1.0],
        'e': [0.6666666666666666, 0.8888888888888888],
        'f': [0.8888888888888888, 1.3333333333333333]}
        sage: ptsdict_2=matroids_plot_helpers.addnontripts(['a','b','c'],
        ....:         ['d','e','f','g','h'],ptsdict)
        sage: setprint(ptsdict_2)
        {'a': [0, 0], 'b': [1, 2], 'c': [0, 2], 'd': [0.6666666666666666, 1.0],
        'e': [0.6666666666666666, 0.8888888888888888],
        'f': [0.8888888888888888, 1.3333333333333333],
        'g': [0.2222222222222222, 1.0],
        'h': [0.5185185185185185, 0.5555555555555555]}

    .. NOTE::

            This method does NOT do any checks.
    """
    tripts = [list(ptsdict[p]) for p in tripts_labels]
    pairs = [[0, 1], [1, 2], [0, 2]]
    q = [tripts]
    num = len(nontripts_labels)
    gridpts = [[float((tripts[0][0]+tripts[1][0]+tripts[2][0])/3),
               float(tripts[0][1]+tripts[1][1]+tripts[2][1])/3]]
    n = 0
    while n < num+1:
        g = trigrid(q[0])
        q.extend([[g[0], q[0][pairs[0][0]], q[0][pairs[0][1]]],
                  [g[0], q[0][pairs[1][0]], q[0][pairs[1][1]]],
                  [g[0], q[0][pairs[2][0]], q[0][pairs[2][1]]]])
        q.remove(q[0])
        gridpts.extend(g[1:])
        if n == 0:
            n = n + 4
        else:
            n = n + 3
    for j, p in enumerate(nontripts_labels):
        ptsdict[p] = tuple(gridpts[j])
    return ptsdict


def createline(ptsdict, ll, lineorders2=None):
    """
    Return ordered lists of coordinates of points to be traversed to draw a
    2D line.

    INPUT:

    - ``ptsdict`` -- dictionary containing keys and their (x,y) position as
      values
    - ``ll`` -- list of keys in ``ptsdict`` through which a line is to be
      drawn
    - ``lineorders2`` -- (optional) list of ordered lists of keys in
      ``ptsdict`` such that if ll is setwise same as any of these then points
      corresponding to values of the keys will be traversed in that order thus
      overriding internal order deciding heuristic

    OUTPUT: a tuple containing 4 elements in this order:

    1. Ordered list of x-coordinates of values of keys in ``ll`` specified in
       ptsdict.
    2. Ordered list of y-coordinates of values of keys in ``ll`` specified
       in ptsdict.
    3. Ordered list of interpolated x-coordinates of points through which a
       line can be drawn.
    4. Ordered list of interpolated y-coordinates of points through which a
       line can be drawn.

    EXAMPLES::

        sage: from sage.matroids import matroids_plot_helpers
        sage: ptsdict = {'a':(1,3),'b':(2,1),'c':(4,5),'d':(5,2)}
        sage: x,y,x_i,y_i = matroids_plot_helpers.createline(ptsdict,
        ....: ['a','b','c','d'])
        sage: [len(x), len(y), len(x_i), len(y_i)]
        [4, 4, 100, 100]
        sage: G = line(zip(x_i, y_i), color='black', thickness=3, zorder=1)             # needs sage.plot
        sage: G += points(zip(x, y), color='black', size=300, zorder=2)                 # needs sage.plot
        sage: G.show()                                                                  # needs sage.plot
        sage: x,y,x_i,y_i = matroids_plot_helpers.createline(ptsdict,
        ....: ['a','b','c','d'],lineorders2=[['b','a','c','d'],
        ....: ['p','q','r','s']])
        sage: [len(x), len(y), len(x_i), len(y_i)]
        [4, 4, 100, 100]
        sage: G = line(zip(x_i, y_i), color='black', thickness=3, zorder=1)             # needs sage.plot
        sage: G += points(zip(x, y), color='black', size=300, zorder=2)                 # needs sage.plot
        sage: G.show()                                                                  # needs sage.plot

    .. NOTE::

            This method does NOT do any checks.
    """
    x, lo = line_hasorder(ll, lineorders2)
    flip = False
    if not x:
        # convert dictionary to list of lists
        linepts = [list(ptsdict[i]) for i in ll]
        xpts = [xx[0] for xx in linepts]
        ypts = [yy[1] for yy in linepts]
        xdim = (float(max(xpts))-float(min(xpts)))
        ydim = (float(max(ypts))-float(min(ypts)))
        if xdim > ydim:
            sortedind = sorted(range(len(xpts)), key=lambda k: float(xpts[k]))
        else:
            sortedind = sorted(range(len(ypts)), key=lambda k: float(ypts[k]))
            flip = True
        sortedlinepts = [linepts[i] for i in sortedind]
        sortedx = [k[0] for k in sortedlinepts]
        sortedy = [k[1] for k in sortedlinepts]
    else:
        linepts = [list(ptsdict[i]) for i in lo]
        sortedx = [k[0] for k in linepts]
        sortedy = [k[1] for k in linepts]

    if flip:
        tck, u = scipy.interpolate.splprep([sortedy, sortedx], s=0.0, k=2)
        y_i, x_i = scipy.interpolate.splev(np.linspace(0, 1, 100), tck)
    else:
        tck, u = scipy.interpolate.splprep([sortedx, sortedy], s=0.0, k=2)
        x_i, y_i = scipy.interpolate.splev(np.linspace(0, 1, 100), tck)
    return sortedx, sortedy, x_i, y_i


def slp(M1, pos_dict=None, B=None):
    """
    Return simple matroid, loops and parallel elements of given matroid.

    INPUT:

    - ``M1`` -- matroid
    - ``pos_dict`` -- (optional) dictionary containing non loopy elements of
      ``M`` as keys and their (x,y) positions as keys. While simplifying the
      matroid, all except one element in a parallel class that is also
      specified in ``pos_dict`` will be retained.
    - ``B`` -- (optional) a basis of M1 that has been chosen for placement on
      vertices of triangle

    OUTPUT: a tuple containing 3 elements in this order:

    1. Simple matroid corresponding to ``M1``.
    2. Loops of matroid ``M1``.
    3. Elements that are in `M1.groundset()` but not in ground set of 1 or
       in 2

    EXAMPLES::

        sage: from sage.matroids import matroids_plot_helpers
        sage: from sage.matroids.advanced import setprint
        sage: M1 = Matroid(ring=GF(2), matrix=[[1, 0, 0, 0, 1, 1, 1,0,1,0,1],
        ....:                                  [0, 1, 0, 1, 0, 1, 1,0,0,1,0],
        ....:                                  [0, 0, 1, 1, 1, 0, 1,0,0,0,0]])
        sage: [M,L,P] = matroids_plot_helpers.slp(M1)                                   # needs sage.rings.finite_rings
        sage: M.is_simple()                                                             # needs sage.rings.finite_rings
        True
        sage: setprint([L,P])                                                           # needs sage.rings.finite_rings
        [{10, 8, 9}, {7}]
        sage: M1 = Matroid(ring=GF(2), matrix=[[1, 0, 0, 0, 1, 1, 1,0,1,0,1],
        ....:                                  [0, 1, 0, 1, 0, 1, 1,0,0,1,0],
        ....:                                  [0, 0, 1, 1, 1, 0, 1,0,0,0,0]])
        sage: posdict = {8: (0, 0),  1: (2, 0),  2: (1, 2),  3: (1.5, 1.0),
        ....:            4: (0.5, 1.0),  5: (1.0, 0.0), 6: (1.0, 0.6666666666666666)}
        sage: [M,L,P] = matroids_plot_helpers.slp(M1, pos_dict=posdict)                 # needs sage.rings.finite_rings
        sage: M.is_simple()                                                             # needs sage.rings.finite_rings
        True
        sage: setprint([L,P])                                                           # needs sage.rings.finite_rings
        [{0, 10, 9}, {7}]

    .. NOTE::

            This method does NOT do any checks.
    """
    L = set(M1.loops())
    nP = L | set(M1.simplify().groundset())
    P = set(M1.groundset())-nP
    if P:
        if pos_dict is not None:
            pcls = list(set([frozenset(set(M1.closure([p])) - L)
                             for p in list(P)]))
            newP = []
            for pcl in pcls:
                pcl_in_dict = [p for p in list(pcl) if p in pos_dict.keys()]
                newP.extend(list(pcl-set([pcl_in_dict[0]])))
            return [M1.delete(L | set(newP)), L, set(newP)]
        elif B is not None:
            pcls = list(set([frozenset(set(M1.closure([p])) - L)
                             for p in list(P)]))
            newP = []
            for pcl in pcls:
                pcl_list = list(pcl)
                pcl_in_basis = [p for p in pcl_list if p in B]
                if pcl_in_basis:
                    newP.extend(list(pcl - set([pcl_in_basis[0]])))
                else:
                    newP.extend(list(pcl - set([pcl_list[0]])))
            return [M1.delete(L | set(newP)), L, set(newP)]
        else:
            return [M1.delete(L | P), L, P]
    else:
        return [M1.delete(L | P), L, P]


def addlp(M, M1, L, P, ptsdict, G=None, limits=None):
    """
    Return a graphics object containing loops (in inset) and parallel elements
    of matroid.

    INPUT:

    - ``M`` -- matroid
    - ``M1`` -- a simple matroid corresponding to ``M``
    - ``L`` -- list of elements in ``M.groundset()`` that are loops of matroid
      ``M``
    - ``P`` -- list of elements in ``M.groundset()`` not in
      ``M.simplify.groundset()`` or ``L``
    - ``ptsdict`` -- dictionary containing elements in ``M.groundset()`` not
      necessarily containing elements of ``L``
    - ``G`` -- (optional) a sage graphics object to which loops and parallel
      elements of matroid `M` added
    - ``limits`` -- (optional) current axes limits [xmin,xmax,ymin,ymax]

    OUTPUT: a 2-tuple containing:

    1. A sage graphics object containing loops and parallel elements of
       matroid ``M``
    2. axes limits array

    EXAMPLES::

        sage: from sage.matroids import matroids_plot_helpers
        sage: M = Matroid(ring=GF(2), matrix=[[1, 0, 0, 0, 1, 1, 1,0,1],
        ....:                                 [0, 1, 0, 1, 0, 1, 1,0,0],
        ....:                                 [0, 0, 1, 1, 1, 0, 1,0,0]])
        sage: [M1,L,P] = matroids_plot_helpers.slp(M)                                   # needs sage.rings.finite_rings
        sage: G, lims = matroids_plot_helpers.addlp(M,M1,L,P,{0:(0,0)})                 # needs sage.plot sage.rings.finite_rings
        sage: G.show(axes=False)                                                        # needs sage.plot sage.rings.finite_rings

    .. NOTE::

            This method does NOT do any checks.
    """
    if G is None:
        G = Graphics()
    # deal with loops
    if L:
        loops = L
        looptext = ", ".join(str(l) for l in loops)
        if limits is None:
            rectx = -1
            recty = -1
        else:
            rectx = limits[0]
            recty = limits[2]-1
        rectw = 0.5 + 0.4*len(loops) + 0.5  # controlled based on len(loops)
        recth = 0.6
        G += polygon2d([[rectx, recty], [rectx, recty+recth],
                        [rectx+rectw, recty+recth], [rectx+rectw, recty]],
                       color='black', fill=False, thickness=4)
        G += text(looptext, (rectx+0.5, recty+0.3), color='black',
                  fontsize=13)
        G += point((rectx+0.2, recty+0.3), color=Color('#BDBDBD'), size=300,
                   zorder=2)
        G += text('Loop(s)', (rectx+0.5+0.4*len(loops)+0.1, recty+0.3),
                  fontsize=13, color='black')
        limits = tracklims(limits, [rectx, rectx+rectw], [recty, recty+recth])
    # deal with parallel elements
    if P:
        # create list of lists where inner lists are parallel classes
        pcls = []
        gnd = sorted(M1.groundset())
        for g in gnd:
            pcl = [g]
            for p in P:
                if M.rank([g, p]) == 1:
                    pcl.extend([p])
            pcls.append(pcl)
        ext_gnd = list(M.groundset())
        for pcl in pcls:
            if len(pcl) > 1:
                basept = list(ptsdict[pcl[0]])
                if len(pcl) <= 2:
                    # add side by side
                    ptsdict[pcl[1]] = (basept[0], basept[1]-0.13)
                    G += points(zip([basept[0]], [basept[1]-0.13]),
                                color=Color('#BDBDBD'), size=300, zorder=2)
                    G += text(pcl[0], (float(basept[0]),
                              float(basept[1])), color='black',
                              fontsize=13)
                    G += text(pcl[1], (float(basept[0]),
                              float(basept[1])-0.13), color='black',
                              fontsize=13)
                    limits = tracklims(limits, [basept[0]], [basept[1]-0.13])
                else:
                    # add in a bracket
                    pce = sorted([str(kk) for kk in pcl])
                    l = newlabel(set(ext_gnd))
                    ext_gnd.append(l)
                    G += text(l+'={ '+", ".join(pce)+' }', (float(basept[0]),
                              float(basept[1]-0.2)-0.034), color='black',
                              fontsize=13)
                    G += text(l, (float(basept[0]),
                              float(basept[1])), color='black',
                              fontsize=13)
                    limits = tracklims(limits, [basept[0]],
                                       [(basept[1]-0.2)-0.034])
    return G, limits


def line_hasorder(l, lodrs=None):
    """
    Determine if an order is specified for a line.

    INPUT:

    - ``l`` -- a line specified as a list of ground set elements
    - ``lordrs`` -- (optional) list of lists each specifying an order on
      a subset of ground set elements that may or may not correspond to a
      line in geometric representation

    OUTPUT: a tuple containing 2 elements in this order:

    1. A boolean indicating whether there is any list in ``lordrs`` that is
       setwise equal to ``l``.
    2. A list specifying an order on ``set(l)`` if 1. is True, otherwise
       an empty list.

    EXAMPLES::

        sage: from sage.matroids import matroids_plot_helpers
        sage: matroids_plot_helpers.line_hasorder(['a','b','c','d'],
        ....: [['a','c','d','b'],['p','q','r']])
        (True, ['a', 'c', 'd', 'b'])
        sage: matroids_plot_helpers.line_hasorder(['a','b','c','d'],
        ....: [['p','q','r'],['l','m','n','o']])
        (False, [])

    .. NOTE::

            This method does NOT do any checks.
    """
    if lodrs is not None:
        set_l = Set(l)
        if lodrs:
            for i in lodrs:
                if Set(i) == set_l:
                    return True, i
    return False, []


def lineorders_union(lineorders1, lineorders2):
    """
    Return a list of ordered lists of ground set elements that corresponds to
    union of two sets of ordered lists of ground set elements in a sense.

    INPUT:

    - ``lineorders1`` -- list of ordered lists specifying orders on subsets
      of ground set
    - ``lineorders2`` -- list of ordered lists specifying orders subsets of
      ground set

    OUTPUT:

    A list of ordered lists of ground set elements that are (setwise) in only
    one of ``lineorders1`` or ``lineorders2`` along with the ones in
    lineorder2 that are setwise equal to any list in lineorders1.

    EXAMPLES::

        sage: from sage.matroids import matroids_plot_helpers
        sage: matroids_plot_helpers.lineorders_union([['a','b','c'],
        ....: ['p','q','r'],['i','j','k','l']],[['r','p','q']])
        [['a', 'b', 'c'], ['p', 'q', 'r'], ['i', 'j', 'k', 'l']]
    """
    if lineorders1 is not None and lineorders2 is not None:
        lineorders = lineorders1
        for order in lineorders2:
            x, lo = line_hasorder(order, lineorders1)
            if not x:
                lineorders.append(order)
                lineorders.remove(lo)
        return lineorders
    elif lineorders1 is None and lineorders2 is not None:
        return lineorders2
    elif lineorders1 is not None:
        return lineorders1
    else:
        return None


def posdict_is_sane(M1, pos_dict):
    """
    Return a boolean establishing sanity of ``posdict`` wrt matroid ``M``.

    INPUT:

    - ``M1`` -- matroid
    - ``posdict`` -- dictionary mapping ground set elements to (x,y)
      positions

    OUTPUT:

    A boolean that is ``True`` if posdict indeed has all the required elements
    to plot the geometric elements, otherwise ``False``.

    EXAMPLES::

        sage: from sage.matroids import matroids_plot_helpers
        sage: M1 = Matroid(ring=GF(2), matrix=[[1, 0, 0, 0, 1, 1, 1,0,1,0,1],
        ....:                                  [0, 1, 0, 1, 0, 1, 1,0,0,1,0],
        ....:                                  [0, 0, 1, 1, 1, 0, 1,0,0,0,0]])
        sage: pos_dict = {0: (0, 0),  1: (2, 0),  2: (1, 2),  3: (1.5, 1.0),
        ....: 4: (0.5, 1.0),  5: (1.0, 0.0), 6: (1.0, 0.6666666666666666)}
        sage: matroids_plot_helpers.posdict_is_sane(M1,pos_dict)                        # needs sage.rings.finite_rings
        True
        sage: pos_dict = {1: (2, 0),  2: (1, 2),  3: (1.5, 1.0),
        ....:             4: (0.5, 1.0), 5: (1.0, 0.0), 6: (1.0, 0.6666666666666666)}
        sage: matroids_plot_helpers.posdict_is_sane(M1,pos_dict)                        # needs sage.rings.finite_rings
        False

    .. NOTE::

            This method does NOT do any checks. ``M1`` is assumed to be a
            matroid and ``posdict`` is assumed to be a dictionary.
    """
    L = set(M1.loops())
    nP = L | set(M1.simplify().groundset())
    P = set(M1.groundset())-nP
    pcls = list(set([frozenset(set(M1.closure([p])) - L) for p in list(P)]))
    for pcl in pcls:
        pcl_list = list(pcl)
        if not any(x in pos_dict for x in pcl_list):
            return False
    allP = []
    for pcl in pcls:
        allP.extend(list(pcl))
    return all(x in pos_dict
               for x in list(set(M1.groundset()) - (L | set(allP))))


def tracklims(lims, x_i=[], y_i=[]):
    """
    Return modified limits list.

    INPUT:

    - ``lims`` -- list with 4 elements ``[xmin,xmax,ymin,ymax]``
    - ``x_i`` -- new x values to track
    - ``y_i`` -- new y values to track

    OUTPUT: list with 4 elements ``[xmin,xmax,ymin,ymax]``

    EXAMPLES::

        sage: from sage.matroids import matroids_plot_helpers
        sage: matroids_plot_helpers.tracklims([0,5,-1,7],[1,2,3,6,-1],
        ....: [-1,2,3,6])
        [-1, 6, -1, 7]

    .. NOTE::

            This method does NOT do any checks.
    """
    if lims is not None and lims[0] is not None and lims[1] is not None and \
       lims[2] is not None and lims[3] is not None:
        lims = [min(min(x_i), lims[0]), max(max(x_i), lims[1]),
                min(min(y_i), lims[2]), max(max(y_i), lims[3])]
    else:
        lims = [min(x_i), max(x_i), min(y_i), max(y_i)]
    return lims


def geomrep(M1, B1=None, lineorders1=None, pd=None, sp=False):
    """
    Return a sage graphics object containing geometric representation of
    matroid M1.

    INPUT:

    - ``M1`` -- matroid
    - ``B1`` -- (optional) list of elements in ``M1.groundset()`` that
      correspond to a basis of ``M1`` and will be placed as vertices of the
      triangle in the geometric representation of ``M1``
    - ``lineorders1`` -- (optional) list of ordered lists of elements of
      ``M1.groundset()`` such that if a line in geometric representation is
      setwise same as any of these then points contained will be traversed in
      that order thus overriding internal order deciding heuristic
    - ``pd`` -- (optional) dictionary mapping ground set elements to their
      (x,y) positions
    - ``sp`` -- (optional) if ``True``, a positioning dictionary and line orders
      will be placed in ``M._cached_info``

    OUTPUT:

    A sage graphics object of type <class 'sage.plot.graphics.Graphics'> that
    corresponds to the geometric representation of the matroid.

    EXAMPLES::

        sage: from sage.matroids import matroids_plot_helpers
        sage: M = matroids.catalog.P7()
        sage: G = matroids_plot_helpers.geomrep(M)                                      # needs sage.plot
        sage: G.show(xmin=-2, xmax=3, ymin=-2, ymax=3)                                  # needs sage.plot
        sage: M = matroids.catalog.P7()
        sage: G = matroids_plot_helpers.geomrep(M, lineorders1=[['f','e','d']])         # needs sage.plot
        sage: G.show(xmin=-2, xmax=3, ymin=-2, ymax=3)                                  # needs sage.plot

    .. NOTE::

            This method does NOT do any checks.
    """
    G = Graphics()
    # create lists of loops and parallel elements and simplify given matroid
    [M, L, P] = slp(M1, pos_dict=pd, B=B1)
    if B1 is None:
        B1 = list(M.basis())
    M._cached_info = M1._cached_info

    if M.rank() == 0:
        limits = None
        loops = L
        looptext = ", ".join(str(l) for l in loops)
        rectx = -1
        recty = -1
        rectw = 0.5 + 0.4 * len(loops) + 0.5  # controlled based on len(loops)
        recth = 0.6
        G += polygon2d([[rectx, recty], [rectx, recty+recth],
                        [rectx+rectw, recty+recth], [rectx+rectw, recty]],
                       color='black', fill=False, thickness=4)
        G += text(looptext, (rectx+0.5, recty+0.3), color='black',
                  fontsize=13)
        G += point((rectx+0.2, recty+0.3), color=Color('#BDBDBD'), size=300,
                   zorder=2)
        G += text('Loop(s)', (rectx+0.5+0.4*len(loops)+0.1, recty+0.3),
                  fontsize=13, color='black')
        limits = tracklims(limits, [rectx, rectx+rectw], [recty, recty+recth])
        G.axes(False)
        G.axes_range(xmin=limits[0]-0.5, xmax=limits[1]+0.5,
                     ymin=limits[2]-0.5, ymax=limits[3]+0.5)
        return G
    elif M.rank() == 1:
        if M._cached_info is not None and \
           'plot_positions' in M._cached_info.keys() and \
           M._cached_info['plot_positions'] is not None:
            pts = M._cached_info['plot_positions']
        else:
            pts = {}
            gnd = sorted(M.groundset())
        pts[gnd[0]] = (1, float(2)/3)
        G += point((1, float(2)/3), size=300, color=Color('#BDBDBD'), zorder=2)
        pt = [1, float(2)/3]
        if not P:
            G += text(gnd[0], (float(pt[0]), float(pt[1])), color='black',
                      fontsize=13)
        pts2 = pts
        # track limits [xmin,xmax,ymin,ymax]
        pl = [list(x) for x in pts2.values()]
        lims = tracklims([None, None, None, None], [pnt[0] for pnt in pl],
                         [pnt[1] for pnt in pl])
    elif M.rank() == 2:
        nB1 = set(M.groundset()) - set(B1)
        bline = [j for j in nB1 if M.is_dependent([j, B1[0], B1[1]])]
        interval = len(bline)+1
        if M._cached_info is not None and \
           'plot_positions' in M._cached_info.keys() and \
           M._cached_info['plot_positions'] is not None:
            pts2 = M._cached_info['plot_positions']
        else:
            pts2 = {}
            pts2[B1[0]] = (0, 0)
            pts2[B1[1]] = (2, 0)
            lpt = list(pts2[B1[0]])
            rpt = list(pts2[B1[1]])
            for k in range(len(bline)):
                cc = (float(1)/interval)*(k+1)
                pts2[bline[k]] = (cc*lpt[0]+(1-cc)*rpt[0],
                                  cc*lpt[1]+(1-cc)*rpt[1])
            if sp:
                M._cached_info['plot_positions'] = pts2
        # track limits [xmin,xmax,ymin,ymax]
        pl = [list(x) for x in pts2.values()]
        lims = tracklims([None, None, None, None], [pt[0] for pt in pl],
                         [pt[1] for pt in pl])
        bline.extend(B1)
        ptsx, ptsy, x_i, y_i = createline(pts2, bline, lineorders1)
        lims = tracklims(lims, x_i, y_i)
        G += line(zip(x_i, y_i), color='black', thickness=3, zorder=1)
        pels = [p for p in pts2
                if any(M1.rank([p, q]) == 1 for q in P)]
        allpts = [list(pts2[i]) for i in M.groundset()]
        xpts = [float(k[0]) for k in allpts]
        ypts = [float(k[1]) for k in allpts]
        G += points(zip(xpts, ypts), color=Color('#BDBDBD'), size=300,
                    zorder=2)
        for i in pts2:
            if i not in pels:
                pt = list(pts2[i])
                G += text(i, (float(pt[0]), float(pt[1])), color='black',
                          fontsize=13)
    else:
        if M._cached_info is None or \
           'plot_positions' not in M._cached_info.keys() or \
           M._cached_info['plot_positions'] is None:
            (pts, trilines,
             nontripts, curvedlines) = it(M1, B1,
                                          list(set(M.groundset())-set(B1)),
                                          list(set(L) | set(P)))
            pts2 = addnontripts([B1[0], B1[1], B1[2]], nontripts, pts)
            trilines.extend(curvedlines)
        else:
            pts2 = M._cached_info['plot_positions']
            trilines = [list(set(x).difference(L | P))
                        for x in M1.flats(2) if len(list(x)) >= 3]
        pl = [list(x) for x in pts2.values()]
        lims = tracklims([None, None, None, None], [pt[0] for pt in pl],
                         [pt[1] for pt in pl])
        j = 0
        for ll in trilines:
            if len(ll) >= 3:
                ptsx, ptsy, x_i, y_i = createline(pts2, ll, lineorders1)
                lims = tracklims(lims, x_i, y_i)
                G += line(zip(x_i, y_i), color='black', thickness=3, zorder=1)
        pels = [p for p in pts2
                if any(M1.rank([p, q]) == 1 for q in P)]
        allpts = [list(pts2[i]) for i in M.groundset()]
        xpts = [float(k[0]) for k in allpts]
        ypts = [float(k[1]) for k in allpts]
        G += points(zip(xpts, ypts), color=Color('#BDBDBD'), size=300,
                    zorder=2)
        for i in pts2:
            if i not in pels:
                pt = list(pts2[i])
                G += text(i, (float(pt[0]), float(pt[1])), color='black',
                          fontsize=13)
        if sp:
            M1._cached_info['plot_positions'] = pts2
            M1._cached_info['plot_lineorders'] = lineorders1
    # deal with loops and parallel elements
    G, lims = addlp(M1, M, L, P, pts2, G, lims)
    G.axes(False)
    G.axes_range(xmin=lims[0]-0.5, xmax=lims[1]+0.5, ymin=lims[2]-0.5,
                 ymax=lims[3]+0.5)
    return G
