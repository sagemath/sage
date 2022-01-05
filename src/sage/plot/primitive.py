"""
Plotting primitives
"""
#*****************************************************************************
#       Copyright (C) 2006 Alex Clemesha <clemesha@gmail.com>,
#                          William Stein <wstein@gmail.com>,
#                     2008 Mike Hansen <mhansen@gmail.com>,
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
from sage.misc.fast_methods import WithEqualityById
from sage.structure.sage_object import SageObject
from sage.misc.verbose import verbose

class GraphicPrimitive(WithEqualityById, SageObject):
    """
    Base class for graphics primitives, e.g., things that knows how to draw
    themselves in 2D.

    EXAMPLES:

    We create an object that derives from GraphicPrimitive::

        sage: P = line([(-1,-2), (3,5)])
        sage: P[0]
        Line defined by 2 points
        sage: type(P[0])
        <class 'sage.plot.line.Line'>

    TESTS::

        sage: hash(circle((0,0),1))  # random
        42
    """
    def __init__(self, options):
        """
        Create a base class GraphicsPrimitive.  All this does is
        set the options.

        EXAMPLES:

        We indirectly test this function::

            sage: from sage.plot.primitive import GraphicPrimitive
            sage: GraphicPrimitive({})
            Graphics primitive
        """
        self._options = options


    def _allowed_options(self):
        """
        Return the allowed options for a graphics primitive.

        OUTPUT:

            - a reference to a dictionary.

        EXAMPLES::

            sage: from sage.plot.primitive import GraphicPrimitive
            sage: GraphicPrimitive({})._allowed_options()
            {}
        """
        return {}

    def plot3d(self, **kwds):
        """
        Plots 3D version of 2D graphics object.  Not implemented
        for base class.

        EXAMPLES::

            sage: from sage.plot.primitive import GraphicPrimitive
            sage: G=GraphicPrimitive({})
            sage: G.plot3d()
            Traceback (most recent call last):
            ...
            NotImplementedError: 3D plotting not implemented for Graphics primitive
        """
        raise NotImplementedError("3D plotting not implemented for %s" % self._repr_())

    def _plot3d_options(self, options=None):
        """
        Translate 2D plot options into 3D plot options.

        EXAMPLES::

            sage: P = line([(-1,-2), (3,5)], alpha=.5, thickness=4)
            sage: p = P[0]; p
            Line defined by 2 points
            sage: q=p.plot3d()
            sage: q.thickness
            4
            sage: q.texture.opacity
            0.5
        """
        if options is None:
            options = self.options()
        options_3d = {}
        if 'rgbcolor' in options:
            options_3d['rgbcolor'] = options['rgbcolor']
            del options['rgbcolor']
        if 'alpha' in options:
            options_3d['opacity'] = options['alpha']
            del options['alpha']

        for o in ('legend_color', 'legend_label', 'zorder'):
            if o in options:
                del options[o]

        if len(options) != 0:
            raise NotImplementedError("Unknown plot3d equivalent for {}".format(
                                      ", ".join(options.keys())))
        return options_3d

    def set_zorder(self, zorder):
        """
        Set the layer in which to draw the object.

        EXAMPLES::

            sage: P = line([(-2,-3), (3,4)], thickness=4)
            sage: p=P[0]
            sage: p.set_zorder(2)
            sage: p.options()['zorder']
            2
            sage: Q = line([(-2,-4), (3,5)], thickness=4,zorder=1,hue=.5)
            sage: P+Q # blue line on top
            Graphics object consisting of 2 graphics primitives
            sage: q=Q[0]
            sage: q.set_zorder(3)
            sage: P+Q # teal line on top
            Graphics object consisting of 2 graphics primitives
            sage: q.options()['zorder']
            3
        """
        self._options['zorder'] = zorder

    def set_options(self, new_options):
        """
        Change the options to `new_options`.

        EXAMPLES::

            sage: from sage.plot.circle import Circle
            sage: c = Circle(0,0,1,{})
            sage: c.set_options({'thickness': 0.6})
            sage: c.options()
            {'thickness': 0.6...}
        """
        if new_options is not None:
            self._options = new_options

    def options(self):
        """
        Return the dictionary of options for this graphics primitive.

        By default this function verifies that the options are all
        valid; if any aren't, then a verbose message is printed with level 0.

        EXAMPLES::

            sage: from sage.plot.primitive import GraphicPrimitive
            sage: GraphicPrimitive({}).options()
            {}
        """
        from sage.plot.graphics import do_verify
        from sage.plot.colors import hue
        O = dict(self._options)
        if do_verify:
            A = self._allowed_options()
            t = False
            K = list(A) + ['xmin', 'xmax', 'ymin', 'ymax', 'axes']
            for k in O.keys():
                if not k in K:
                    do_verify = False
                    verbose("WARNING: Ignoring option '%s'=%s"%(k,O[k]), level=0)
                    t = True
            if t:
                s = "\nThe allowed options for %s are:\n"%self
                K.sort()
                for k in K:
                    if k in A:
                        s += "    %-15s%-60s\n"%(k,A[k])
                verbose(s, level=0)


        if 'hue' in O:
            t = O['hue']
            if not isinstance(t, (tuple,list)):
                t = [t,1,1]
            O['rgbcolor'] = hue(*t)
            del O['hue']
        return O

    def _repr_(self):
        """
        String representation of this graphics primitive.

        EXAMPLES::

            sage: from sage.plot.primitive import GraphicPrimitive
            sage: GraphicPrimitive({})._repr_()
            'Graphics primitive'
        """
        return "Graphics primitive"



class GraphicPrimitive_xydata(GraphicPrimitive):
    def get_minmax_data(self):
        """
        Returns a dictionary with the bounding box data.

        EXAMPLES::

            sage: d = polygon([[1,2], [5,6], [5,0]], rgbcolor=(1,0,1))[0].get_minmax_data()
            sage: d['ymin']
            0.0
            sage: d['xmin']
            1.0

        ::

            sage: d = point((3, 3), rgbcolor=hue(0.75))[0].get_minmax_data()
            sage: d['xmin']
            3.0
            sage: d['ymin']
            3.0

        ::

            sage: l = line([(100, 100), (120, 120)])[0]
            sage: d = l.get_minmax_data()
            sage: d['xmin']
            100.0
            sage: d['xmax']
            120.0

        """
        from sage.plot.plot import minmax_data
        return minmax_data(self.xdata, self.ydata, dict=True)

