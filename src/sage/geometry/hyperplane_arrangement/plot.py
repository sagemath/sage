# sage_setup: distribution = sagemath-polyhedra
"""
Plotting of Hyperplane Arrangements

PLOT OPTIONS:

Beside the usual plot options (enter ``plot?``), the plot command for
hyperplane arrangements includes the following:

- ``hyperplane_colors`` -- color or list of colors, one for each
  hyperplane (default: equally spread range of hues)

- ``hyperplane_labels`` -- boolean, ``'short'``, ``'long'`` (default:
  ``False``).  If ``False``, no labels are shown; if 'short' or 'long',
  the hyperplanes are given short or long labels, respectively.  If
  ``True``, the hyperplanes are given long labels.

- ``label_colors`` -- color or list of colors, one for each hyperplane
  (default: black)

- ``label_fontsize`` -- size for hyperplane_label font (default:
  ``14``); this does not work for 3d plots

- ``label_offsets`` -- amount be which labels are offset from
  h.point() for each hyperplane h.  The format is different for each
  dimension: if the hyperplanes have dimension 0, the offset can be a
  single number or a list of numbers, one for each hyperplane; if the
  hyperplanes have dimension 1, the offset can be a single 2-tuple, or
  a list of 2-tuples, one for each hyperplane; if the hyperplanes have
  dimension 2, the offset can be a single 3-tuple or a list of
  3-tuples, one for each hyperplane.  (Defaults: 0-dim: ``0.1``,
  1-dim: ``(0,1)``, 2-dim: ``(0,0,0.2)``).

- ``hyperplane_legend`` -- boolean, ``'short'``, ``'long'`` (default:
  ``'long'``; in 3-d: ``False``).  If ``False``, no legend is shown;
  if ``True``, ``'short'``, or ``'long'``, the legend is shown with
  the default, long, or short labeling, respectively. (For
  arrangements of lines or planes, only.)

- ``hyperplane_opacities`` -- a number or list of numbers, one for each
  hyperplane, between 0 and 1;  only applies to 3d plots

- ``point_sizes`` -- number or list of numbers, one for each hyperplane
  giving the sizes of points in a zero-dimensional arrangement
  (default: ``50``)

- ``ranges`` -- range for the parameters or a list of ranges of
  parameters, one for each hyperplane, for the parametric plots of the
  hyperplanes.  If a single positive number `r` is given for
  ``ranges``, then all parameters run from -r to r.  Otherwise, for a
  line in the plane, the range has the form ``[a,b]`` (default:
  ``[-3,3]``), and for a plane in 3-space, the range has the form
  ``[[a,b],[c,d]]`` (default: ``[[-3,3],[-3,3]]``). The ranges are
  centered around ``hyperplane_arrangement.point()``.

EXAMPLES::

    sage: H3.<x,y,z> = HyperplaneArrangements(QQ)
    sage: A = H3([(1,0,0), 0], [(0,0,1), 5])
    sage: A.plot(hyperplane_opacities=0.5, hyperplane_labels=True,                      # needs sage.plot
    ....:        hyperplane_legend=False)
    Graphics3d Object

    sage: c = H3([(1,0,0),0], [(0,0,1),5])
    sage: c.plot(ranges=10)                                                             # needs sage.plot
    Graphics3d Object
    sage: c.plot(ranges=[[9.5,10], [-3,3]])                                             # needs sage.plot
    Graphics3d Object
    sage: c.plot(ranges=[[[9.5,10], [-3,3]], [[-6,6], [-5,5]]])                         # needs sage.plot
    Graphics3d Object


    sage: H2.<s,t> = HyperplaneArrangements(QQ)
    sage: h = H2([(1,1),0], [(1,-1),0], [(0,1),2])
    sage: h.plot(ranges=20)                                                             # needs sage.plot
    Graphics object consisting of 3 graphics primitives
    sage: h.plot(ranges=[-1, 10])                                                       # needs sage.plot
    Graphics object consisting of 3 graphics primitives
    sage: h.plot(ranges=[[-1, 1], [-5, 5], [-1, 10]])                                   # needs sage.plot
    Graphics object consisting of 3 graphics primitives

    sage: a = hyperplane_arrangements.coordinate(3)
    sage: opts = {'hyperplane_colors':['yellow', 'green', 'blue']}
    sage: opts['hyperplane_labels'] = True
    sage: opts['label_offsets'] = [(0,2,2), (2,0,2), (2,2,0)]
    sage: opts['hyperplane_legend'] = False
    sage: opts['hyperplane_opacities'] = 0.7
    sage: a.plot(**opts)                                                                # needs sage.plot
    Graphics3d Object
    sage: opts['hyperplane_labels'] = 'short'
    sage: a.plot(**opts)                                                                # needs sage.plot
    Graphics3d Object

    sage: H.<u> = HyperplaneArrangements(QQ)
    sage: pts = H(3*u+4, 2*u+5, 7*u+1)
    sage: pts.plot(hyperplane_colors=['yellow','black','blue'])                         # needs sage.plot
    Graphics object consisting of 3 graphics primitives
    sage: pts.plot(point_sizes=[50,100,200], hyperplane_colors='blue')                  # needs sage.plot
    Graphics object consisting of 3 graphics primitives

    sage: H.<x,y,z> = HyperplaneArrangements(QQ)
    sage: a = H(x, y+1, y+2)
    sage: a.plot(hyperplane_labels=True, label_colors='blue', label_fontsize=18)        # needs sage.plot
    Graphics3d Object
    sage: a.plot(hyperplane_labels=True, label_colors=['red','green','black'])          # needs sage.plot
    Graphics3d Object
"""
from copy import copy
from colorsys import hsv_to_rgb

from sage.misc.lazy_import import lazy_import
lazy_import("sage.plot.plot3d.parametric_plot3d", "parametric_plot3d")
lazy_import("sage.plot.plot3d.shapes2", "text3d")
lazy_import("sage.plot.graphics", "Graphics")
lazy_import("sage.plot.line", "line")
lazy_import("sage.plot.text", "text")
lazy_import("sage.plot.point", "point")
lazy_import("sage.plot.plot", "parametric_plot")


def plot(hyperplane_arrangement, **kwds):
    r"""
    Return a plot of the hyperplane arrangement.

    If the arrangement is in 4 dimensions but inessential, a plot of
    the essentialization is returned.

    .. NOTE::

        This function is available as the
        :meth:`~sage.geometry.hyperplane_arrangement.arrangement.HyperplaneArrangementElement.plot`
        method of hyperplane arrangements. You should not call this
        function directly, only through the method.

    INPUT:

    - ``hyperplane_arrangement`` -- the hyperplane arrangement to plot

    - ``**kwds`` -- plot options: see
      :mod:`sage.geometry.hyperplane_arrangement.plot`

    OUTPUT: a graphics object of the plot

    EXAMPLES::

        sage: B = hyperplane_arrangements.semiorder(4)
        sage: B.plot()                                                                  # needs sage.combinat sage.plot
        Displaying the essentialization.
        Graphics3d Object
    """
    N = len(hyperplane_arrangement)
    dim = hyperplane_arrangement.dimension()
    if hyperplane_arrangement.base_ring().characteristic() != 0:
        raise NotImplementedError('must be a field of characteristic 0')
    elif dim == 4:
        if not hyperplane_arrangement.is_essential():
            print('Displaying the essentialization.')
            hyperplane_arrangement = hyperplane_arrangement.essentialization()
    elif dim not in [1,2,3]: # revise to handle 4d
        return # silently
    # handle extra keywords
    if 'hyperplane_colors' in kwds:
        hyp_colors = kwds.pop('hyperplane_colors')
        if not isinstance(hyp_colors, list): # we assume its a single color then
            hyp_colors = [hyp_colors] * N
    else:
        HSV_tuples = [(i*1.0/N, 0.8, 0.9) for i in range(N)]
        hyp_colors = [hsv_to_rgb(*x) for x in HSV_tuples]
    if 'hyperplane_labels' in kwds:
        hyp_labels = kwds.pop('hyperplane_labels')
        has_hyp_label = True
        if not isinstance(hyp_labels, list): # we assume its a boolean then
            hyp_labels = [hyp_labels] * N
        relabeled = []
        for i in range(N):
            if hyp_labels[i] in [True,'long']:
                relabeled.append(True)
            else:
                relabeled.append(str(i))
        hyp_labels = relabeled
    else:
        has_hyp_label = False
    if 'label_colors' in kwds:
        label_colors = kwds.pop('label_colors')
        has_label_color = True
        if not isinstance(label_colors, list): # we assume its a single color then
            label_colors = [label_colors] * N
    else:
        has_label_color = False
    if 'label_fontsize' in kwds:
        label_fontsize = kwds.pop('label_fontsize')
        has_label_fontsize = True
        if not isinstance(label_fontsize, list): # we assume its a single size then
            label_fontsize = [label_fontsize] * N
    else:
        has_label_fontsize = False
    if 'label_offsets' in kwds:
        has_offsets = True
        offsets = kwds.pop('label_offsets')
    else:
        has_offsets = False # give default values below
    hyperplane_legend = kwds.pop('hyperplane_legend', 'long' if dim < 3 else False)
    if 'hyperplane_opacities' in kwds:
        hyperplane_opacities = kwds.pop('hyperplane_opacities')
        has_opacity = True
        if not isinstance(hyperplane_opacities, list): # we assume a single number then
            hyperplane_opacities = [hyperplane_opacities] * N
    else:
        has_opacity = False
    point_sizes = kwds.pop('point_sizes', 50)
    if not isinstance(point_sizes, list):
        point_sizes = [point_sizes] * N
    if 'ranges' in kwds:
        ranges_set = True
        ranges = kwds.pop('ranges')
        if type(ranges) not in [list,tuple]: # ranges is a single number
            ranges = [ranges] * N
        # So ranges is some type of list.
        elif dim == 2: # arrangement of lines in the plane
            if type(ranges[0]) not in [list,tuple]: # a single interval
                ranges = [ranges] * N
        elif dim == 3: # arrangement of planes in 3-space
            if type(ranges[0][0]) not in [list,tuple]:
                ranges = [ranges] * N
        elif dim not in [2,3]: # ranges is not an option unless dim is 2 or 3
            ranges_set = False
        else: # a list of intervals, one for each hyperplane is given
            pass # ranges does not need to be modified
    else:
        ranges_set = False # give default values below
    # the extra keywords have now been handled
    # now handle the legend
    if dim in [1,2]: # points on a line or lines in the plane
        if hyperplane_legend in [True,'long']:
            hyps = hyperplane_arrangement.hyperplanes()
            legend_labels = [hyps[i]._latex_() for i in range(N)]
        elif hyperplane_legend == 'short' :
            legend_labels = [str(i) for i in range(N)]
    else: # dim==3,  arrangement of planes in 3-space
        if hyperplane_legend in [True, 'long']:
            legend3d = legend_3d(hyperplane_arrangement, hyp_colors, 'long')
        elif hyperplane_legend == 'short':
            legend3d = legend_3d(hyperplane_arrangement, hyp_colors, 'short')
    ## done handling the legend
    ## now create the plot
    p = Graphics()
    for i in range(N):
        newk = copy(kwds)
        if has_hyp_label:
            newk['hyperplane_label'] = hyp_labels[i]
            if has_offsets:
                if not isinstance(offsets, list):
                    newk['label_offset'] = offsets
                else:
                    newk['label_offset'] = offsets[i]
        else:
            newk['hyperplane_label'] = False
        if has_label_color:
            newk['label_color'] = label_colors[i]
        if has_label_fontsize:
            newk['label_fontsize'] = label_fontsize[i]
        if has_opacity:
            newk['opacity'] = hyperplane_opacities[i]
        if dim == 1:
            newk['point_size'] = point_sizes[i]
        if dim in [1,2] and hyperplane_legend: # more options than T/F
            newk['legend_label'] = legend_labels[i]
        if ranges_set:
            newk['ranges'] = ranges[i]
        p += plot_hyperplane(hyperplane_arrangement[i], rgbcolor=hyp_colors[i], **newk)
    if dim == 1:
        if hyperplane_legend: # there are more options than T/F
            p.legend(True)
        return p
    elif dim == 2:
        if hyperplane_legend: # there are more options than T/F
            p.legend(True)
        return p
    else: # dim==3
        if hyperplane_legend: # there are more options than T/F
            return p, legend3d
        else:
            return p


def plot_hyperplane(hyperplane, **kwds):
    r"""
    Return the plot of a single hyperplane.

    INPUT:

    - ``**kwds`` -- plot options: see below

    OUTPUT: a graphics object of the plot

    .. RUBRIC:: Plot Options

    Beside the usual plot options (enter ``plot?``), the plot command for
    hyperplanes includes the following:

    - ``hyperplane_label`` -- boolean value or string (default: ``True``);
      if ``True``, the hyperplane is labeled with its equation, if a
      string, it is labeled by that string, otherwise it is not
      labeled

    - ``label_color`` -- (default: ``'black'``) color for hyperplane_label

    - ``label_fontsize`` -- size for ``hyperplane_label`` font (default: 14)
      (does not work in 3d, yet)

    - ``label_offset`` -- (default: 0-dim: 0.1, 1-dim: (0,1),
      2-dim: (0,0,0.2)); amount by which label is offset from
      ``hyperplane.point()``

    - ``point_size`` -- (default: 50) size of points in a zero-dimensional
      arrangement or of an arrangement over a finite field

    - ``ranges`` -- range for the parameters for the parametric plot of the
      hyperplane. If a single positive number ``r`` is given for the
      value of ``ranges``, then the ranges for all parameters are set to
      `[-r, r]`.  Otherwise, for a line in the plane, ``ranges`` has the
      form ``[a, b]`` (default: [-3,3]), and for a plane in 3-space, the
      ``ranges`` has the form ``[[a, b], [c, d]]`` (default: [[-3,3],[-3,3]]).
      (The ranges are centered around ``hyperplane.point()``.)

    EXAMPLES::

        sage: H1.<x> = HyperplaneArrangements(QQ)
        sage: a = 3*x + 4
        sage: a.plot()    # indirect doctest                                            # needs sage.plot
        Graphics object consisting of 3 graphics primitives
        sage: a.plot(point_size=100, hyperplane_label='hello')                          # needs sage.plot
        Graphics object consisting of 3 graphics primitives

        sage: H2.<x,y> = HyperplaneArrangements(QQ)
        sage: b = 3*x + 4*y + 5
        sage: b.plot()                                                                  # needs sage.plot
        Graphics object consisting of 2 graphics primitives
        sage: b.plot(ranges=(1,5), label_offset=(2,-1))                                 # needs sage.plot
        Graphics object consisting of 2 graphics primitives
        sage: opts = {'hyperplane_label': True, 'label_color': 'green',
        ....:         'label_fontsize': 24, 'label_offset': (0,1.5)}
        sage: b.plot(**opts)                                                            # needs sage.plot
        Graphics object consisting of 2 graphics primitives

        sage: # needs sage.plot
        sage: H3.<x,y,z> = HyperplaneArrangements(QQ)
        sage: c = 2*x + 3*y + 4*z + 5
        sage: c.plot()
        Graphics3d Object
        sage: c.plot(label_offset=(1,0,1), color='green', label_color='red',
        ....:        frame=False)
        Graphics3d Object
        sage: d = -3*x + 2*y + 2*z + 3
        sage: d.plot(opacity=0.8)
        Graphics3d Object
        sage: e = 4*x + 2*z + 3
        sage: e.plot(ranges=[[-1,1],[0,8]], label_offset=(2,2,1), aspect_ratio=1)
        Graphics3d Object
    """
    if hyperplane.base_ring().characteristic():
        raise NotImplementedError('base field must have characteristic zero')
    elif hyperplane.dimension() not in [0, 1, 2]: # dimension of hyperplane, not ambient space
        raise ValueError('can only plot hyperplanes in dimensions 1, 2, 3')
    # handle extra keywords
    if 'hyperplane_label' in kwds:
        hyp_label = kwds.pop('hyperplane_label')
        if not hyp_label:
            has_hyp_label = False
        else:
            has_hyp_label = True
    else: # default
        hyp_label = True
        has_hyp_label = True
    if has_hyp_label:
        if hyp_label: # then label hyperplane with its equation
            if hyperplane.dimension() == 2: # jmol does not like latex
                label = hyperplane._repr_linear(include_zero=False)
            else:
                label = hyperplane._latex_()
        else:
            label = hyp_label # a string
    if 'label_color' in kwds:
        label_color = kwds.pop('label_color')
    else:
        label_color = 'black'
    if 'label_fontsize' in kwds:
        label_fontsize = kwds.pop('label_fontsize')
    else:
        label_fontsize = 14
    if 'label_offset' in kwds:
        has_offset = True
        label_offset = kwds.pop('label_offset')
    else:
        has_offset = False # give default values below
    if 'point_size' in kwds:
        pt_size = kwds.pop('point_size')
    else:
        pt_size = 50
    if 'ranges' in kwds:
        ranges_set = True
        ranges = kwds.pop('ranges')
    else:
        ranges_set = False  # give default values below
    # the extra keywords have now been handled
    # now create the plot
    if hyperplane.dimension() == 0:  # a point on a line
        x, = hyperplane.A()
        d = hyperplane.b()
        p = point((d/x,0), size=pt_size, **kwds)
        if has_hyp_label:
            if not has_offset:
                label_offset = 0.1
            p += text(label, (d/x,label_offset),
                    color=label_color,fontsize=label_fontsize)
            p += text('',(d/x,label_offset+0.4)) # add space at top
        if 'ymax' not in kwds:
            kwds['ymax'] = 0.5
    elif hyperplane.dimension() == 1: # a line in the plane
        pnt = hyperplane.point()
        w = hyperplane.linear_part().matrix()
        from sage.symbolic.ring import SR
        t = SR.var('t')
        if ranges_set:
            if isinstance(ranges, (list, tuple)):
                t0, t1 = ranges
            else:  # ranges should be a single positive number
                t0, t1 = -ranges, ranges
        else: # default
            t0, t1 = -3, 3
        p = parametric_plot(pnt + t * w[0], (t, t0, t1), **kwds)
        if has_hyp_label:
            if has_offset:
                b0, b1 = label_offset
            else:
                b0, b1 = 0, 0.2
            label = text(label,(pnt[0] + b0, pnt[1] + b1),
                    color=label_color,fontsize=label_fontsize)
            p += label
    elif hyperplane.dimension() == 2: # a plane in 3-space
        pnt = hyperplane.point()
        w = hyperplane.linear_part().matrix()
        from sage.symbolic.ring import SR
        s, t = SR.var('s t')
        if ranges_set:
            if isinstance(ranges, (list, tuple)):
                s0, s1 = ranges[0]
                t0, t1 = ranges[1]
            else:  # ranges should be a single positive integers
                s0, s1 = -ranges, ranges
                t0, t1 = -ranges, ranges
        else:  # default
            s0, s1 = -3, 3
            t0, t1 = -3, 3
        p = parametric_plot3d(pnt+s*w[0]+t*w[1], (s,s0,s1), (t,t0,t1), **kwds)
        if has_hyp_label:
            if has_offset:
                b0, b1, b2 = label_offset
            else:
                b0, b1, b2 = 0, 0, 0
            label = text3d(label,(pnt[0]+b0, pnt[1]+b1, pnt[2]+b2),
                    color=label_color, fontsize=label_fontsize)
            p += label
    return p


def legend_3d(hyperplane_arrangement, hyperplane_colors, length):
    r"""
    Create plot of a 3d legend for an arrangement of planes in 3-space.

    The ``length`` parameter determines whether short or long labels
    are used in the legend.

    INPUT:

    - ``hyperplane_arrangement`` -- a hyperplane arrangement

    - ``hyperplane_colors`` -- list of colors

    - ``length`` -- either ``'short'`` or ``'long'``

    OUTPUT: a graphics object

    EXAMPLES::

        sage: a = hyperplane_arrangements.semiorder(3)
        sage: from sage.geometry.hyperplane_arrangement.plot import legend_3d
        sage: legend_3d(a, list(colors.values())[:6], length='long')                    # needs sage.combinat sage.plot
        Graphics object consisting of 6 graphics primitives

        sage: b = hyperplane_arrangements.semiorder(4)
        sage: c = b.essentialization()
        sage: legend_3d(c, list(colors.values())[:12], length='long')                   # needs sage.combinat sage.plot
        Graphics object consisting of 12 graphics primitives

        sage: legend_3d(c, list(colors.values())[:12], length='short')                  # needs sage.combinat sage.plot
        Graphics object consisting of 12 graphics primitives

        sage: p = legend_3d(c, list(colors.values())[:12], length='short')              # needs sage.combinat sage.plot
        sage: p.set_legend_options(ncol=4)                                              # needs sage.combinat sage.plot
        sage: type(p)                                                                   # needs sage.combinat sage.plot
        <class 'sage.plot.graphics.Graphics'>
    """
    if hyperplane_arrangement.dimension() != 3:
        raise ValueError('arrangements must be in 3-space')
    hyps = hyperplane_arrangement.hyperplanes()
    N = len(hyperplane_arrangement)
    if length == 'short':
        labels = ['  ' + str(i) for i in range(N)]
    else:
        labels = ['  ' + hyps[i]._repr_linear(include_zero=False) for i in
                  range(N)]
    p = Graphics()
    for i in range(N):
        p += line([(0,0),(0,0)], color=hyperplane_colors[i], thickness=8,
                legend_label=labels[i], axes=False)
    p.set_legend_options(title='Hyperplanes', loc='center', labelspacing=0.4,
            fancybox=True, font_size='x-large', ncol=2)
    p.legend(True)
    return p
