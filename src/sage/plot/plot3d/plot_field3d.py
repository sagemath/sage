"""
Plotting 3D Fields
"""
#*****************************************************************************
#       Copyright (C) 2009 Jason Grout <jason-sage@creativetrax.com>
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

from sage.arith.srange import srange
from sage.plot.misc import setup_for_eval_on_grid
from sage.modules.free_module_element import vector
from sage.plot.plot import plot

def plot_vector_field3d(functions, xrange, yrange, zrange,
                        plot_points=5, colors='jet', center_arrows=False, **kwds):
    r"""
    Plot a 3d vector field

    INPUT:

    - ``functions`` - a list of three functions, representing the x-,
      y-, and z-coordinates of a vector

    - ``xrange``, ``yrange``, and ``zrange`` - three tuples of the
      form (var, start, stop), giving the variables and ranges for each axis

    - ``plot_points`` (default 5) - either a number or list of three
      numbers, specifying how many points to plot for each axis

    - ``colors`` (default 'jet') - a color, list of colors (which are
      interpolated between), or matplotlib colormap name, giving the coloring
      of the arrows.  If a list of colors or a colormap is given,
      coloring is done as a function of length of the vector

    - ``center_arrows`` (default False) - If True, draw the arrows
      centered on the points; otherwise, draw the arrows with the tail
      at the point

    - any other keywords are passed on to the plot command for each arrow

    EXAMPLES::

        sage: x,y,z=var('x y z')
        sage: plot_vector_field3d((x*cos(z),-y*cos(z),sin(z)), (x,0,pi), (y,0,pi), (z,0,pi))
        Graphics3d Object
        sage: plot_vector_field3d((x*cos(z),-y*cos(z),sin(z)), (x,0,pi), (y,0,pi), (z,0,pi),colors=['red','green','blue'])
        Graphics3d Object
        sage: plot_vector_field3d((x*cos(z),-y*cos(z),sin(z)), (x,0,pi), (y,0,pi), (z,0,pi),colors='red')
        Graphics3d Object
        sage: plot_vector_field3d((x*cos(z),-y*cos(z),sin(z)), (x,0,pi), (y,0,pi), (z,0,pi),plot_points=4)
        Graphics3d Object
        sage: plot_vector_field3d((x*cos(z),-y*cos(z),sin(z)), (x,0,pi), (y,0,pi), (z,0,pi),plot_points=[3,5,7])
        Graphics3d Object
        sage: plot_vector_field3d((x*cos(z),-y*cos(z),sin(z)), (x,0,pi), (y,0,pi), (z,0,pi),center_arrows=True)
        Graphics3d Object

    TESTS:

    This tests that :trac:`2100` is fixed in a way compatible with this command::

        sage: plot_vector_field3d((x*cos(z),-y*cos(z),sin(z)), (x,0,pi), (y,0,pi), (z,0,pi),center_arrows=True,aspect_ratio=(1,2,1))
        Graphics3d Object
    """
    (ff,gg,hh), ranges = setup_for_eval_on_grid(functions, [xrange, yrange, zrange], plot_points)
    xpoints, ypoints, zpoints = [srange(*r, include_endpoint=True) for r in ranges]
    points = [vector((i,j,k)) for i in xpoints for j in ypoints for k in zpoints]
    vectors = [vector((ff(*point), gg(*point), hh(*point))) for point in points]

    try:
        from matplotlib.cm import get_cmap
        cm = get_cmap(colors)
    except (TypeError, ValueError):
        cm = None
    if cm is None:
        if isinstance(colors, (list, tuple)):
            from matplotlib.colors import LinearSegmentedColormap
            cm = LinearSegmentedColormap.from_list('mymap',colors)
        else:
            cm = lambda x: colors

    max_len = max(v.norm() for v in vectors)
    scaled_vectors = [v/max_len for v in vectors]

    if center_arrows:
        G = sum([plot(v,color=cm(v.norm()),**kwds).translate(p-v/2) for v,p in zip(scaled_vectors, points)])
        G._set_extra_kwds(kwds)
        return G
    else:
        G = sum([plot(v,color=cm(v.norm()),**kwds).translate(p) for v,p in zip(scaled_vectors, points)])
        G._set_extra_kwds(kwds)
        return G
