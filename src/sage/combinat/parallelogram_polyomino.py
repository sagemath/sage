# sage.doctest: needs sage.combinat sage.modules
r"""
Parallelogram Polyominoes
=========================

The goal of this module is to give some tools to manipulate the
parallelogram polyominoes.
"""
# *****************************************************************************
#  Copyright (C) 2014,2015 Adrien Boussicault (boussica@labri.fr),
#  Copyright (C) 2016 Patxi Laborde-Zubieta (plaborde@labri.fr),
#  Copyright (C) 2019 Henri Derycke (henri.derycke@ens-lyon.org),
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************
from __future__ import annotations

from sage.structure.list_clone import ClonableList
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.set_factories import (SetFactory, ParentWithSetFactory,
                                          TopMostParentPolicy)
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.sets.set import Set
from sage.misc.lazy_attribute import lazy_class_attribute
from sage.misc.lazy_attribute import lazy_attribute
from sage.sets.disjoint_union_enumerated_sets import (
    DisjointUnionEnumeratedSets
)
from sage.rings.integer import Integer
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.sets.family import Family
from sage.sets.positive_integers import PositiveIntegers
from sage.misc.cachefunc import cached_method
from sage.misc.latex import latex
from copy import deepcopy
from sage.matrix.constructor import matrix
from sage.combinat.combinat import catalan_number
from sage.combinat.combinatorial_map import combinatorial_map
from sage.functions.trig import cos, sin
from sage.misc.functional import sqrt

from sage.misc.lazy_import import lazy_import
lazy_import("sage.plot.graphics", "Graphics")
lazy_import("sage.plot.line", "line")
lazy_import("sage.plot.text", "text")
lazy_import("sage.plot.point", "point")

import pprint


class LocalOptions:
    r"""
    This class allow to add local options to an object.
    LocalOptions is like a dictionary, it has keys and values that represent
    options and the values associated to the option. This is useful to
    decorate an object with some optional informations.

    :class:`LocalOptions` should be used as follow.

    INPUT:

    - ``name`` -- the name of the LocalOptions

    - ``<options>=dict(...)`` -- dictionary specifying an option

    The options are specified by keyword arguments with their values
    being a dictionary which describes the option. The
    allowed/expected keys in the dictionary are:

    - ``checker`` -- a function for checking whether a particular value for
      the option is valid
    - ``default`` -- the default value of the option
    - ``values`` -- dictionary of the legal values for this option (this
      automatically defines the corresponding ``checker``); this dictionary
      gives the possible options, as keys, together with a brief description
      of them

    ::

        sage: from sage.combinat.parallelogram_polyomino import LocalOptions
        sage: o = LocalOptions(
        ....:     'Name Example',
        ....:     delim=dict(
        ....:         default='b',
        ....:         values={'b':'the option b', 'p':'the option p'}
        ....:     )
        ....: )
        sage: class Ex:
        ....:     options=o
        ....:     def _repr_b(self): return "b"
        ....:     def _repr_p(self): return "p"
        ....:     def __repr__(self): return self.options._dispatch(
        ....:         self, '_repr_','delim'
        ....:     )
        sage: e = Ex(); e
        b
        sage: e.options(delim='p'); e
        p

    This class is temporary, in the future, this class should be integrated in
    sage.structure.global_options.py. We should split global_option in two
    classes LocalOptions and GlobalOptions.
    """

    def __init__(self, name='', **options):
        r"""
        Construct a new LocalOptions.

        INPUT:

        - ``name`` -- the name of the LocalOptions

        - ``<options>=dict(...)`` -- dictionary specifying an option

        The options are specified by keyword arguments with their values
        being a dictionary which describes the option. The
        allowed/expected keys in the dictionary are:

        - ``checker`` -- a function for checking whether a particular value for
          the option is valid
        - ``default`` -- the default value of the option
        - ``values`` -- dictionary of the legal values for this option (this
          automatically defines the corresponding ``checker``); this dictionary
          gives the possible options, as keys, together with a brief
          description of them.

        EXAMPLES::

            sage: from sage.combinat.parallelogram_polyomino import (
            ....:     LocalOptions
            ....: )
            sage: o = LocalOptions(
            ....:     "Name Example",
            ....:     tikz_options=dict(
            ....:         default='toto',
            ....:         values=dict(
            ....:             toto='name',
            ....:             x="3"
            ....:         )
            ....:     ),
            ....:     display=dict(
            ....:         default='list',
            ....:         values=dict(
            ....:             list="list representation",
            ....:             diagram="diagram representation"
            ....:         )
            ....:     )
            ....: )
        """
        self._name = name
        self._available_options = options
        self._options = {}
        for key in self._available_options:
            self._options[key] = self._available_options[key]["default"]

    def __repr__(self) -> str:
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.parallelogram_polyomino import (
            ....:     LocalOptions
            ....: )
            sage: o = LocalOptions(
            ....:     "Name Example",
            ....:     tikz_options=dict(
            ....:         default='toto',
            ....:         values=dict(
            ....:             toto='name',
            ....:             x="3"
            ....:         )
            ....:     ),
            ....:     display=dict(
            ....:         default='list',
            ....:         values=dict(
            ....:             list="list representation",
            ....:             diagram="diagram representation"
            ....:         )
            ....:     )
            ....: )
            sage: o
            Current options for Name Example
              - display:      'list'
              - tikz_options: 'toto'
        """
        options = list(self._options)
        if options == []:
            return 'Current options for {}'.format(self._name)

        options.sort()
        width = 1 + max(len(key) for key in options)
        txt = '\n'.join('  - {:{}} {}'.format(key + ':', width, pprint.pformat(self[key]))
                        for key in options)
        return 'Current options for {}\n{}'.format(self._name, txt)

    def __setitem__(self, key, value):
        r"""
        The ``__setitem__`` method is used to change the current values of the
        options. It also checks that the supplied options are valid and changes
        any alias to its generic value.

        INPUT:

        - ``key`` -- an option

        - ``value`` -- the value

        EXAMPLES::

            sage: from sage.combinat.parallelogram_polyomino import (
            ....:     LocalOptions
            ....: )
            sage: o = LocalOptions(
            ....:     "Name Example",
            ....:     tikz_options=dict(
            ....:         default='toto',
            ....:         values=dict(
            ....:             toto='name',
            ....:             x="3"
            ....:         )
            ....:     ),
            ....:     display=dict(
            ....:         default='list',
            ....:         values=dict(
            ....:             list="list representation",
            ....:             diagram="diagram representation"
            ....:         )
            ....:     ),
            ....:     size=dict(
            ....:         default=1,
            ....:         checker=lambda x : x in NN
            ....:     )
            ....: )
            sage: o("display")
            'list'
            sage: o["display"]="diagram"
            sage: o("display")
            'diagram'
            sage: o["display"]="?"
            Current value : diagram
            {'default': 'list', 'values':
            {'diagram': 'diagram representation',
            'list': 'list representation'}}
            sage: o("size")
            1
            sage: o["size"]=3
            sage: o("size")
            3
            sage: o["size"]=-6
        """
        assert (key in self._available_options)
        if value == "?":
            res = "Current value : " + str(self._options[key])
            option_key = self._available_options[key]
            if "values" in option_key:
                res += "\n" + pprint.pformat(self._available_options[key])
            print(res)
        else:
            available_options = self._available_options
            if "values" in available_options:
                assert (value in self._available_options[key]["values"])
            if "checker" in available_options:
                assert (available_options["checker"](value))
            self._options[key] = value

    def __call__(self, *get_values, **options):
        r"""
        Get or set value of the option ``option``.

        INPUT:

        - ``get_values`` -- the options to be printed

        - ``<options>=dict(...)`` -- dictionary specifying an option see
          :class:`LocalOptions` for more details

        EXAMPLES::

            sage: from sage.combinat.parallelogram_polyomino import (
            ....:     LocalOptions
            ....: )
            sage: o = LocalOptions(
            ....:     "Name Example",
            ....:     tikz_options=dict(
            ....:         default='toto',
            ....:         values=dict(
            ....:             toto='name',
            ....:             x="3"
            ....:         )
            ....:     ),
            ....:     display=dict(
            ....:         default='list',
            ....:         values=dict(
            ....:             list="list representation",
            ....:             diagram="diagram representation"
            ....:         )
            ....:     )
            ....: )
            sage: o("display")
            'list'
            sage: o(display='diagram')
            sage: o("display")
            'diagram'
            sage: o(display='?')
            Current value : diagram
            {'default': 'list', 'values':
            {'diagram': 'diagram representation',
            'list': 'list representation'}}
        """
        for key in options:
            value = options[key]
            self.__setitem__(key, value)
        for key in get_values:
            return self.__getitem__(key)

    def __getitem__(self, key):
        r"""
        Return the current value of the option ``key``.

        INPUT:

        - ``key`` -- an option

        EXAMPLES::

            sage: from sage.combinat.parallelogram_polyomino import (
            ....:     LocalOptions
            ....: )
            sage: o = LocalOptions(
            ....:     "Name Example",
            ....:     tikz_options=dict(
            ....:         default='toto',
            ....:         values=dict(
            ....:             toto='name',
            ....:             x="3"
            ....:         )
            ....:     ),
            ....:     display=dict(
            ....:         default='list',
            ....:         values=dict(
            ....:             list="list representation",
            ....:             diagram="diagram representation"
            ....:         )
            ....:     )
            ....: )
            sage: o["display"]
            'list'
        """
        return self._options[key]

    def __iter__(self):
        r"""
        A generator for the options in ``self``.

        EXAMPLES::

            sage: from sage.combinat.parallelogram_polyomino import (
            ....:     LocalOptions
            ....: )
            sage: o = LocalOptions(
            ....:     "Name Example",
            ....:     tikz_options=dict(
            ....:         default='toto',
            ....:         values=dict(
            ....:             toto='name',
            ....:             x="3"
            ....:         )
            ....:     ),
            ....:     display=dict(
            ....:         default='list',
            ....:         values=dict(
            ....:             list="list representation",
            ....:             diagram="diagram representation"
            ....:         )
            ....:     )
            ....: )
            sage: all(key in ['tikz_options','display'] for key in o)
            True
        """
        return self._available_options.__iter__()

    def keys(self) -> list[str]:
        r"""
        Return the list of the options in ``self``.

        EXAMPLES::

            sage: from sage.combinat.parallelogram_polyomino import (
            ....:     LocalOptions
            ....: )
            sage: o = LocalOptions(
            ....:     "Name Example",
            ....:     tikz_options=dict(
            ....:         default='toto',
            ....:         values=dict(
            ....:             toto='name',
            ....:             x="3"
            ....:         )
            ....:     ),
            ....:     display=dict(
            ....:         default='list',
            ....:         values=dict(
            ....:             list="list representation",
            ....:             diagram="diagram representation"
            ....:         )
            ....:     )
            ....: )
            sage: keys=o.keys()
            sage: keys.sort()
            sage: keys
            ['display', 'tikz_options']
        """
        return list(self)

    def _dispatch(self, obj, dispatch_to, option, *get_values, **set_values):
        r"""
        The *dispatchable* options are options which dispatch related methods
        of the corresponding class. The format for specifying a dispatchable
        option is to include ``dispatch_to = <option name>`` in the
        specifications for the options and then to add the options to the
        class.

        The ``_dispatch`` method will then call the method named
        ``<option name> + '_' + <current value of option>`` of ``obj``
        with arguments ``*get_values, **set_values``.

        Note that the argument ``self`` is necessary here because the
        dispatcher is a method of the options class and not of ``self``.

        EXAMPLES::

            sage: from sage.combinat.parallelogram_polyomino import LocalOptions
            sage: delim = {'default': 'b',
            ....:          'values': {'b': 'option b', 'p': 'option p'}}
            sage: o = LocalOptions('Name example', delim=delim)
            sage: class Ex:
            ....:     options=o
            ....:     def _repr_b(self): return "b"
            ....:     def _repr_p(self): return "p"
            ....:     def __repr__(self): return self.options._dispatch(
            ....:         self, '_repr_','delim')
            sage: e = Ex(); e
            b
            sage: e.options(delim='p'); e
            p
        """
        assert (option in self._available_options)
        if dispatch_to[-1] == "_":
            dispatch_to = dispatch_to[:-1]
        f = getattr(obj, dispatch_to + "_" + str(self._options[option]))
        return f(*get_values, **set_values)


default_tikz_options = dict(
    scale=1, line_size=1, point_size=3.5, color_line='black',
    color_point='black', color_bounce_0='red', color_bounce_1='blue',
    translation=[0, 0], rotation=0, mirror=None
)
r"""
This is the default TIKZ options.

This option is used to configurate element of a drawing to allow
TIKZ code generation.
"""

ParallelogramPolyominoesOptions = LocalOptions(
    'ParallelogramPolyominoes_size',
    #    module='sage.combinat.parallelogram_polyomino',
    #    doc=r"""
    #    """,
    #    end_doc=r"""
    #    """,
    tikz_options=dict(
        default=default_tikz_options,
        description='the tikz options',
        checker=lambda x: Set(x.keys()).issubset(
            Set(
                [
                    'scale', 'line_size', 'point_size',
                    'color_line', 'color_point', 'translation', 'mirror',
                    'rotation', 'color_bounce_0', 'color_bounce_1',
                ]
            )
        )
    ),
    drawing_components=dict(
        default=dict(diagram=True, tree=False, bounce_0=False, bounce_1=False, bounce_values=False),
        description='Different tree-like tableaux components to draw',
        checker=lambda x: Set(x.keys()).issubset(
            Set(['diagram', 'tree', 'bounce_0', 'bounce_1', 'bounce_values', ])
        )
    ),
    display=dict(
        default='list',
        values=dict(
            list='displayed as list',
            drawing='as a drawing',
        )
    ),
    latex=dict(
        default='drawing',
        values=dict(
            list='displayed as list',
            drawing='as a drawing',
        )
    )
)
r"""
This global option contains all the data needed by the Parallelogram classes
to draw, display in ASCII, compile in latex a parallelogram polyomino.

The available options are:

- tikz_options : this option configurate all the information useful to
  generate TIKZ code. For example, color, line size, etc ...

- drawing_components : this option is used to explain to the system
  which component of the drawing you want to draw. For example,
  you can ask to draw some elements of the following list:
  - the diagram,
  - the tree inside the parallelogram polyomino,
  - the bounce paths inside the parallelogram polyomino,
  - the value of the bounce on each square of a bounce path.

- display : this option is used to configurate the ASCII display.
  The available options are:
  - list : (this is the default value) is used to represent PP as a list
  containing the upper and lower path.
  - drawing : this value is used to explain we want to display an array with
  the PP drawn inside (with connected 1).

- latex : Same as display. The default is "drawing".

See :meth:`ParallelogramPolyomino.get_options` for more details and for an
user use of options.

EXAMPLES::

    sage: from sage.combinat.parallelogram_polyomino import (
    ....:     ParallelogramPolyominoesOptions
    ....: )
    sage: opt = ParallelogramPolyominoesOptions['tikz_options']
    sage: opt
    {'color_bounce_0': 'red',
     'color_bounce_1': 'blue',
     'color_line': 'black',
     'color_point': 'black',
     'line_size': 1,
     'mirror': None,
     'point_size': 3.5,
     'rotation': 0,
     'scale': 1,
     'translation': [0, 0]}
"""


class _drawing_tool:
    r"""
    Technical class to produce TIKZ drawing.

    This class contains some 2D geometric tools to produce some TIKZ
    drawings.

    With that classes you can use options to set up drawing informations.
    Then the class will produce a drawing by using those informations.

    EXAMPLES::

        sage: from sage.combinat.parallelogram_polyomino import (
        ....:     _drawing_tool, default_tikz_options,
        ....:     ParallelogramPolyominoesOptions
        ....: )
        sage: opt = ParallelogramPolyominoesOptions['tikz_options']
        sage: dt = _drawing_tool(opt)
        sage: dt.draw_line([1, 1], [-1, -1])
        '\n  \\draw[color=black, line width=1] (1.000000, 1.000000) --
        (-1.000000, -1.000000);'

        sage: fct = lambda vec: [2*vec[0], vec[1]]
        sage: dt = _drawing_tool(opt, fct)
        sage: dt.draw_line([1, 1], [-1, -1])
        '\n  \\draw[color=black, line width=1] (2.000000, 1.000000) --
        (-2.000000, -1.000000);'

        sage: import copy
        sage: opt = copy.deepcopy(opt)
        sage: opt['mirror'] = [0,1]
        sage: dt = _drawing_tool(opt)
        sage: dt.draw_line([1, 1], [-1, -1])
        '\n  \\draw[color=black, line width=1] (-1.000000, 1.000000) --
        (1.000000, -1.000000);'
    """

    def __init__(self, options, XY=lambda v: v):
        r"""
        Construct a drawing tools to produce some TIKZ drawing.

        INPUT:

        - ``options`` -- drawing options

        - ``XY`` -- a user function to convert vector in other vector
          (default: identity function)

        EXAMPLES::

            sage: from sage.combinat.parallelogram_polyomino import (
            ....:     _drawing_tool, default_tikz_options,
            ....:     ParallelogramPolyominoesOptions
            ....: )
            sage: opt = ParallelogramPolyominoesOptions['tikz_options']
            sage: dt = _drawing_tool(opt)
            sage: dt.draw_line([1, 1], [-1, -1])
            '\n  \\draw[color=black, line width=1] (1.000000, 1.000000) --
            (-1.000000, -1.000000);'
        """
        self._XY = lambda v: XY([float(v[0]), float(v[1])])
        self._translation = options['translation']
        self._mirror = options['mirror']
        self._rotation = options['rotation']
        self._color_line = options['color_line']
        self._line_size = options['line_size']
        self._point_size = options['point_size']
        self._color_point = options['color_point']

    def XY(self, v):
        r"""
        This function give the image of v by some transformation given by the
        drawing option of ``_drawing_tool``.

        The transformation is the composition of rotation, mirror, translation
        and XY user function.

        First we apply XY function, then the translation, then the mirror and
        finally the rotation.

        INPUT:

        - ``v`` -- the vector to transform

        OUTPUT: list of 2 floats encoding a vector

        EXAMPLES::

            sage: from sage.combinat.parallelogram_polyomino import (
            ....:     _drawing_tool, ParallelogramPolyominoesOptions
            ....: )
            sage: opt = ParallelogramPolyominoesOptions['tikz_options']
            sage: dt = _drawing_tool(opt)
            sage: dt.XY([1, 1])
            [1.0, 1.0]

            sage: fct = lambda vec: [2*vec[0], vec[1]]
            sage: dt = _drawing_tool(opt, fct)
            sage: dt.XY([1, 1])
            [2.0, 1.0]

            sage: import copy
            sage: opt = copy.deepcopy(opt)
            sage: opt['mirror'] = [0, 1]
            sage: dt = _drawing_tool(opt)
            sage: dt.XY([1, 1])
            [-1.0, 1.0]
        """
        def translate(pos, v):
            r"""
            Translate a position with a vector.

            INPUT:

            - ``pos`` -- the position to translate

            - ``v`` -- the translation vector

            OUTPUT:

            The translated position.
            """
            return [pos[0] + v[0], pos[1] + v[1]]

        def rotate(pos, angle):
            r"""
            Rotate by `angle` a position around the origin.

            INPUT:

            - ``pos`` -- the position to rotate

            - ``angle`` -- the angle of rotation

            OUTPUT:

            The rotated position.
            """
            [x, y] = pos
            return [x*cos(angle) - y*sin(angle), x*sin(angle) + y*cos(angle)]

        def mirror(pos, axe):
            r"""
            Return the mirror of a position according to a given axe.

            INPUT:

            - ``pos`` -- the position to mirror

            - ``axe`` -- the axe vector

            OUTPUT:

            The mirrored position.
            """
            if axe is None:
                return pos
            if not isinstance(axe, (list, tuple)):
                raise ValueError(
                    "mirror option should be None or a list of two real" +
                    " encoding a 2D vector."
                )
            n = float(sqrt(axe[0]**2 + axe[1]**2))
            axe[0] = float(axe[0]/n)
            axe[1] = float(axe[1]/n)
            sp = (pos[0]*axe[0] + pos[1]*axe[1])
            sn = (- pos[0]*axe[1] + pos[1]*axe[0])
            return [
                sp*axe[0] + sn*axe[1],
                sp*axe[1] - sn*axe[0]
            ]
        return rotate(
            mirror(
                translate(self._XY(v), self._translation),
                self._mirror
            ), self._rotation
        )

    def draw_line(self, v1, v2, color=None, size=None):
        r"""
        Return the TIKZ code for a line.

        INPUT:

        - ``v1`` -- point, The first point of the line

        - ``v2`` -- point, The second point of the line

        - ``color`` -- string (default: ``None``); the color of the line.
          If set to ``None``, the color is chosen according the
          drawing option given by ``_drawing_tool``.

        - ``size`` -- integer (default: ``None``); the size of the line.
          If set to ``None``, the size is chosen according the
          drawing option given by ``_drawing_tool``.

        OUTPUT: the code of a line in TIKZ

        EXAMPLES::

            sage: from sage.combinat.parallelogram_polyomino import (
            ....:     _drawing_tool, ParallelogramPolyominoesOptions
            ....: )
            sage: opt = ParallelogramPolyominoesOptions['tikz_options']
            sage: dt = _drawing_tool(opt)
            sage: dt.draw_line([1, 1], [-1, -1])
            '\n  \\draw[color=black, line width=1] (1.000000, 1.000000) --
            (-1.000000, -1.000000);'
        """
        if color is None:
            color = self._color_line
        if size is None:
            size = self._line_size
        [x1, y1] = self.XY(v1)
        [x2, y2] = self.XY(v2)
        return "\n  \\draw[color=%s, line width=%s] (%f, %f) -- (%f, %f);" % (
            color, size, float(x1), float(y1), float(x2), float(y2)
        )

    def draw_polyline(self, list_of_vertices, color=None, size=None):
        r"""
        Return the TIKZ code for a polyline.

        INPUT:

        - ``list_of_vertices`` -- list of points

        - ``color`` -- string (default: ``None``); the color of the line.
          If set to ``None``, the color is chosen according the
          drawing option given by ``_drawing_tool``.

        - ``size`` -- integer (default: ``None``); the size of the line.
          If set to ``None``, the size is chosen according the
          drawing option given by ``_drawing_tool``.

        OUTPUT: the code of a polyline in TIKZ

        EXAMPLES::

            sage: from sage.combinat.parallelogram_polyomino import (
            ....:     _drawing_tool, ParallelogramPolyominoesOptions
            ....: )
            sage: opt = ParallelogramPolyominoesOptions['tikz_options']
            sage: dt = _drawing_tool(opt)
            sage: dt.draw_polyline([[1, 1], [-1, -1], [0,0]])
            '\n  \\draw[color=black, line width=1] (1.000000, 1.000000) --
            (-1.000000, -1.000000);\n  \\draw[color=black, line width=1]
            (-1.000000, -1.000000) -- (0.000000, 0.000000);'
        """
        res = ""
        for i in range(len(list_of_vertices)-1):
            res += self.draw_line(
                list_of_vertices[i], list_of_vertices[i+1], color, size)
        return res

    def draw_point(self, p1, color=None, size=None):
        r"""
        Return the TIKZ code for a point.

        INPUT:

        - ``p1`` -- a point

        - ``color`` -- string (default: ``None``); the color of the line.
          If set to ``None``, the color is chosen according the
          drawing option given by ``_drawing_tool``.

        - ``size`` -- integer (default: ``None``); the size of the line.
          If set to ``None``, the size is chosen according the
          drawing option given by ``_drawing_tool``.

        OUTPUT: the code of a point in TIKZ

        EXAMPLES::

            sage: from sage.combinat.parallelogram_polyomino import (
            ....:     _drawing_tool, ParallelogramPolyominoesOptions
            ....: )
            sage: opt = ParallelogramPolyominoesOptions['tikz_options']
            sage: dt = _drawing_tool(opt)
            sage: dt.draw_point([1, 1])
            '\n  \\filldraw[color=black] (1.000000, 1.000000) circle (3.5pt);'
        """
        if color is None:
            color = self._color_point
        if size is None:
            size = self._point_size
        [x1, y1] = self.XY(p1)
        return "\n  \\filldraw[color=%s] (%f, %f) circle (%spt);" % (
            color, float(x1), float(y1), size
        )


class ParallelogramPolyomino(ClonableList,
        metaclass=InheritComparisonClasscallMetaclass):
    r"""
    Parallelogram Polyominoes.

    A parallelogram polyomino is a finite connected union of cells
    whose boundary can be decomposed in two paths, the upper and the lower
    paths, which are comprised of north and east unit steps and meet only at
    their starting and final points.

    Parallelogram Polyominoes can be defined with those two paths.

    EXAMPLES::

        sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
        sage: pp
        [[0, 1], [1, 0]]
    """
    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        r"""
        Ensure that parallelogram polyominoes created by the enumerated sets
        are instances of :class:`ParallelogramPolyomino` and have the same
        parent.

        TESTS::

            sage: issubclass(
            ....:     ParallelogramPolyominoes().element_class,
            ....:     ParallelogramPolyomino
            ....: )
            True
            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp.parent()
            Parallelogram polyominoes
            sage: str(type(pp)) == (
            ....:      "<class 'sage.combinat.parallelogram_polyomino" +
            ....:      ".ParallelogramPolyominoes_all_with_category" +
            ....:      ".element_class'>"
            ....: )
            True

            sage: pp1 = ParallelogramPolyominoes()([[0, 1], [1, 0]])
            sage: pp1.parent() is pp.parent()
            True
            sage: type(pp1) is type(pp)
            True

            sage: pp1 = ParallelogramPolyominoes(2)([[0, 1], [1, 0]])
            sage: pp1.parent() is pp.parent()
            True
            sage: type(pp1) is type(pp)
            True
        """
        return cls._auto_parent._element_constructor_(*args, **opts)

    @lazy_class_attribute
    def _auto_parent(cls):
        r"""
        The automatic parent of the elements of this class.

        When calling the constructor of an element of this class, one needs a
        parent. This class attribute specifies which parent is used.

        EXAMPLES::

            sage: ParallelogramPolyomino._auto_parent
            Parallelogram polyominoes
            sage: ParallelogramPolyominoes()([[0, 1], [1, 0]]).parent()
            Parallelogram polyominoes
        """
        return ParallelogramPolyominoes()

    def _ascii_art_(self):
        """
        TESTS::

            sage: ascii_art(ParallelogramPolyomino([[0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1], [1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0]]))
            ***
            ****
             ***
              ***
               **
                *
            sage: ascii_art(ParallelogramPolyomino([[0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1], [1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0]]))
            **
            ***
            ***
            ***
              **
               **
        """
        from sage.typeset.ascii_art import AsciiArt

        data = zip(self.lower_widths(), self.upper_widths())
        txt = []
        for x,y in data:
            txt += [' ' * x + '*' * (y - x)]

        return AsciiArt(txt)

    def _unicode_art_(self):
        """
        TESTS::

            sage: unicode_art(ParallelogramPolyomino([[0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1], [1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0]]))
            ┌┬┬┐
            ├┼┼┼┐
            └┼┼┼┤
             └┼┼┼┐
              └┼┼┤
               └┼┤
                └┘
            sage: unicode_art(ParallelogramPolyomino([[0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1], [1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0]]))
            ┌┬┐
            ├┼┼┐
            ├┼┼┤
            ├┼┼┤
            └┴┼┼┐
              └┼┼┐
               └┴┘
        """
        from sage.typeset.unicode_art import UnicodeArt

        data = list(zip(self.lower_widths(), self.upper_widths()))

        txt = ['┌' + '┬' * (data[0][1] - 1) + '┐']
        for i in range(1, len(data)):
            x1, y1 = data[i-1]
            x2, y2 = data[i]
            line = [' ' * x1]
            if x1 == x2:
                line += ['├']
            else:
                line += ['└' + '┴' * (x2 - x1 - 1) + '┼']
            line += ['┼' * (y1 - x2 - 1)]
            if y1 == y2:
                line += ['┤']
            else:
                line += ['┼' + '┬' * (y2 - y1 - 1) + '┐']
            txt += [''.join(line)]
        txt += [' ' * data[-1][0] + '└' + '┴' * (data[-1][1] - data[-1][0] - 1) + '┘']

        return UnicodeArt(txt, baseline=0)

    def check(self):
        r"""
        This method raises an error if the internal data of the class does not
        represent a parallelogram polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 0, 1, 1, 0, 0, 1, 0, 0]
            ....:     ]
            ....: )
            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp = ParallelogramPolyomino([[1], [1]])

            sage: pp = ParallelogramPolyomino(                # indirect doctest
            ....:     [[1, 0], [0, 1]]
            ....: )
            Traceback (most recent call last):
            ...
            ValueError: the lower and upper paths are crossing

            sage: pp = ParallelogramPolyomino([[1], [0, 1]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the lower and upper paths have different sizes (2 != 1)

            sage: pp = ParallelogramPolyomino([[1], [0]])     # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the two paths have distinct ends

            sage: pp = ParallelogramPolyomino([[0], [1]])     # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the two paths have distinct ends

            sage: pp = ParallelogramPolyomino([[0], [0]])     # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the lower or the upper path can...t be equal to [0]

            sage: pp = ParallelogramPolyomino([[], [0]])      # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the lower or the upper path can...t be equal to []

            sage: pp = ParallelogramPolyomino([[0], []])      # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the lower or the upper path can...t be equal to []

            sage: pp = ParallelogramPolyomino([[], []])       # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the lower or the upper path can...t be equal to []
        """
        lower_path = self.lower_path()
        upper_path = self.upper_path()
        if lower_path == [0] and upper_path == [0]:
            raise ValueError(
                "the lower or the upper path can't be equal to [0]"
            )
        if lower_path == [] or upper_path == []:
            raise ValueError(
                "the lower or the upper path can't be equal to []"
            )
        if len(upper_path) != len(lower_path):
            raise ValueError(
                "the lower and upper paths have different sizes (%s != %s)" % (
                    len(upper_path), len(lower_path)
                )
            )
        p_up = [0, 0]
        p_down = [0, 0]
        for i in range(len(upper_path)-1):
            p_up[1-upper_path[i]] += 1
            p_down[1-lower_path[i]] += 1
            if (p_up[0] <= p_down[0] or p_down[1] <= p_up[1]):
                raise ValueError("the lower and upper paths are crossing")
        p_up[1 - upper_path[-1]] += 1
        p_down[1 - lower_path[-1]] += 1
        if (p_up[0] != p_down[0] or p_up[1] != p_down[1]):
            raise ValueError("the two paths have distinct ends")

    def __hash__(self):
        r"""
        Return the hash code of the parallelogram polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 0, 1, 1, 0, 0, 1, 0, 0]
            ....:     ]
            ....: )
            sage: hash(pp) == hash((
            ....:     (0, 0, 0, 1, 0, 1, 0, 1, 1), (1, 0, 1, 1, 0, 0, 1, 0, 0)
            ....: ))
            True

            sage: PPS = ParallelogramPolyominoes(8)
            sage: D = { PPS[0]: True, PPS[1]: True }
            sage: D[PPS[0]] = False
            sage: import pprint
            sage: pp = pprint.PrettyPrinter()
            sage: pp.pprint(D)
            {[[0, 0, 0, 0, 0, 0, 0, 1], [1, 0, 0, 0, 0, 0, 0, 0]]: False,
             [[0, 0, 0, 0, 0, 0, 1, 1], [1, 0, 0, 0, 0, 0, 1, 0]]: True}
        """
        return hash(tuple(map(tuple, list(self))))

    def __copy__(self):
        r"""
        Copy a parallelogram Polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 0, 1, 1, 0, 0, 1, 0, 0]
            ....:     ]
            ....: )
            sage: pp1 = copy(pp)
            sage: pp1 is pp
            False
            sage: pp1 == pp
            True
            sage: pp1
            [[0, 0, 0, 1, 0, 1, 0, 1, 1], [1, 0, 1, 1, 0, 0, 1, 0, 0]]
        """
        return ParallelogramPolyomino([self.lower_path(), self.upper_path()])

    def __init__(self, parent, value, check=True):
        r"""
        Construct a parallelogram polyomino.

        The input is a pair of lower path and upper path.

        The lower and upper paths of the empty parallelogram polyomino are
        [1] and [1].

        EXAMPLES::

            sage: lower_path = [0, 0, 1, 0, 1, 1]
            sage: upper_path = [1, 1, 0, 1, 0, 0]
            sage: pp = ParallelogramPolyomino([lower_path, upper_path])
            sage: pp
            [[0, 0, 1, 0, 1, 1], [1, 1, 0, 1, 0, 0]]

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp
            [[0, 1], [1, 0]]

            sage: pp = ParallelogramPolyomino([[1], [1]])
            sage: pp
            [[1], [1]]
        """
        ClonableList.__init__(self, parent, value)
        if check:
            if not isinstance(value, (list, tuple)):
                raise ValueError(
                    "value %s must be a list or a tuple" % value)
            self.check()
        self._options = None

    def reflect(self) -> ParallelogramPolyomino:
        r"""
        Return the parallelogram polyomino obtained by switching rows and
        columns.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino([[0,0,0,0,1,1,0,1,0,1], [1,0,1,0,0,1,1,0,0,0]])
            sage: pp.heights(), pp.upper_heights()
            ([4, 3, 2, 3], [0, 1, 3, 3])
            sage: pp = pp.reflect()
            sage: pp.widths(), pp.lower_widths()
            ([4, 3, 2, 3], [0, 1, 3, 3])

            sage: pp = ParallelogramPolyomino([[0,0,0,1,1], [1,0,0,1,0]])
            sage: ascii_art(pp)
            *
            *
            **
            sage: ascii_art(pp.reflect())
            ***
              *

        TESTS::

           sage: pp = ParallelogramPolyomino([[1], [1]])
           sage: pp.reflect()
           [[1], [1]]
        """
        if self.size() == 1:
            return self
        a, b = self
        return ParallelogramPolyomino([[1 - v for v in b],
                                       [1 - v for v in a]])

    def rotate(self) -> ParallelogramPolyomino:
        r"""
        Return the parallelogram polyomino obtained by rotation of 180 degrees.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino([[0,0,0,1,1], [1,0,0,1,0]])
            sage: ascii_art(pp)
            *
            *
            **
            sage: ascii_art(pp.rotate())
            **
             *
             *
        """
        a, b = self
        return ParallelogramPolyomino([b[::-1], a[::-1]])

    def _to_dyck_delest_viennot(self):
        r"""
        Convert to a Dyck word using the Delest-Viennot bijection.

        This bijection is described on page 179 and page 180 Figure 6
        in the article [DeVi1984]_, where it is called the classical
        bijection `\gamma`.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino([[0, 1, 0, 0, 1, 1], [1, 1, 1, 0, 0, 0]])
            sage: pp._to_dyck_delest_viennot()
            [1, 1, 0, 1, 1, 0, 1, 0, 0, 0]

        TESTS::

            sage: pp = ParallelogramPolyomino([[1], [1]])
            sage: pp._to_dyck_delest_viennot()
            []
        """
        from sage.combinat.dyck_word import DyckWord
        dyck = []
        dick_size = self.size() - 1
        if not dick_size:
            return DyckWord([])
        upper_path = self.upper_path()
        lower_path = self.lower_path()
        dyck.append(1 - lower_path[0])
        for i in range(1, dick_size):
            dyck.append(upper_path[i])
            dyck.append(1 - lower_path[i])
        dyck.append(upper_path[dick_size])
        return DyckWord(dyck)

    def _to_dyck_delest_viennot_peaks_valleys(self):
        r"""
        Convert to a Dyck word using the Delest-Viennot bijection `\beta`.

        This bijection is described on page 182 and Figure 8 in the
        article [DeVi1984]_.  It returns the unique Dyck path whose
        peak heights are the column heights and whose valley heights
        are the overlaps between adjacent columns.

        EXAMPLES:

        This is the example in Figure 8 of [DeVi1984]_::

            sage: pp = ParallelogramPolyomino([[0,0,0,0,1,1,0,1,0,1], [1,0,1,0,0,1,1,0,0,0]])
            sage: pp._to_dyck_delest_viennot_peaks_valleys()
            [1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0]

        TESTS::

            sage: pp = ParallelogramPolyomino([[1], [1]])
            sage: pp._to_dyck_delest_viennot_peaks_valleys()
            []
        """
        from sage.combinat.dyck_word import DyckWord
        a = self.heights()
        u = self.upper_heights()
        b = [0] + [a[i]-u[i+1]+u[i]-1 for i in range(len(a)-1)] + [0]
        dyck = []
        for i in range(len(a)):
            dyck.extend([1] * (a[i] - b[i]))
            dyck.extend([0] * (a[i] - b[i + 1]))
        return DyckWord(dyck)

    @combinatorial_map(name="To Dyck word")
    def to_dyck_word(self, bijection=None):
        r"""
        Convert to a Dyck word.

        INPUT:

        - ``bijection`` -- string or ``None`` (default: ``None``); the name of
          the bijection. If it is set to ``None`` then the ``'Delest-Viennot'``
          bijection is used.
          Expected values are ``None``, ``'Delest-Viennot'``, or ``'Delest-Viennot-beta'``.

        OUTPUT: a Dyck word

        EXAMPLES::

            sage: pp = ParallelogramPolyomino([[0, 1, 0, 0, 1, 1], [1, 1, 1, 0, 0, 0]])
            sage: pp.to_dyck_word()
            [1, 1, 0, 1, 1, 0, 1, 0, 0, 0]
            sage: pp.to_dyck_word(bijection='Delest-Viennot')
            [1, 1, 0, 1, 1, 0, 1, 0, 0, 0]

            sage: pp.to_dyck_word(bijection='Delest-Viennot-beta')
            [1, 0, 1, 1, 1, 0, 1, 0, 0, 0]
        """
        if bijection is None or bijection == 'Delest-Viennot':
            return self._to_dyck_delest_viennot()
        if bijection == 'Delest-Viennot-beta':
            return self._to_dyck_delest_viennot_peaks_valleys()
        raise ValueError("the given bijection is not valid")

    @staticmethod
    def _from_dyck_word_delest_viennot(dyck):
        r"""
        Convert a Dyck word to a parallelogram polyomino using the Delest
        Viennot bijection.

        This bijection is described on page 179 and page 180 Figure 6 in
        the article [DeVi1984]_, where it is called the classical
        bijection `\gamma`.

        INPUT:

        - ``dyck`` -- a Dyck word

        OUTPUT: a parallelogram polyomino

        EXAMPLES::

            sage: dyck = DyckWord([1, 1, 0, 1, 1, 0, 1, 0, 0, 0])
            sage: ParallelogramPolyomino._from_dyck_word_delest_viennot(dyck)
            [[0, 1, 0, 0, 1, 1], [1, 1, 1, 0, 0, 0]]

        TESTS::

            sage: gamma = ParallelogramPolyomino._to_dyck_delest_viennot
            sage: gamma_inv = ParallelogramPolyomino._from_dyck_word_delest_viennot
            sage: all(all(D == gamma(gamma_inv(D)) for D in DyckWords(n)) for n in range(7))
            True
        """
        l = [1] + list(dyck) + [0]
        word_up = []
        word_down = []
        for i in range(0, len(l), 2):
            word_up.append(l[i])
            word_down.append(1 - l[i + 1])
        return ParallelogramPolyomino([word_down, word_up])

    @staticmethod
    def _from_dyck_word_delest_viennot_peaks_valleys(dyck):
        r"""
        Convert a Dyck word to a parallelogram polyomino using the Delest
        Viennot bijection `\beta`.

        This bijection is described on page 182 and Figure 8 in
        the article [DeVi1984]_.

        INPUT:

        - ``dyck`` -- a Dyck word

        OUTPUT: a parallelogram polyomino

        EXAMPLES::

            sage: dyck = DyckWord([1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0])
            sage: ParallelogramPolyomino._from_dyck_word_delest_viennot_peaks_valleys(dyck)
            [[0, 0, 0, 0, 1, 1, 0, 1, 0, 1], [1, 0, 1, 0, 0, 1, 1, 0, 0, 0]]

            sage: dyck = DyckWord([1,1,0,1,1,1,1,1,0,0,1,0,0,0,0,0,1,1,1,0,0,1,0,0])
            sage: ParallelogramPolyomino._from_dyck_word_delest_viennot_peaks_valleys(dyck)
            [[0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1], [1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0]]

        TESTS::

            sage: beta = ParallelogramPolyomino._to_dyck_delest_viennot_peaks_valleys
            sage: beta_inv = ParallelogramPolyomino._from_dyck_word_delest_viennot_peaks_valleys
            sage: all(all(D == beta(beta_inv(D)) for D in DyckWords(n)) for n in range(7))
            True
        """
        if not dyck:
            return ParallelogramPolyomino([[1], [1]])
        a = []
        b = [0]
        h = 0
        for i in range(len(dyck)-1):
            if dyck[i] == 1:
                h += 1
                if dyck[i+1] == 0:
                    a.append(h)
            else:
                if dyck[i+1] == 1:
                    b.append(h)
                h -= 1
        b.append(0)
        word_down = []
        word_up = []
        for i in range(len(a)):
            word_down.extend([0]*(a[i]-b[i]) + [1])
            word_up.extend([1]+[0]*(a[i]-b[i+1]))
        return ParallelogramPolyomino([word_down, word_up])

    @staticmethod
    def from_dyck_word(dyck, bijection=None):
        r"""
        Convert a Dyck word to parallelogram polyomino.

        INPUT:

        - ``dyck`` -- a Dyck word

        - ``bijection`` -- string or ``None`` (default: ``None``); the
          bijection to use. See :meth:`to_dyck_word` for more details.

        OUTPUT: a parallelogram polyomino

        EXAMPLES::

            sage: dyck = DyckWord([1, 1, 0, 1, 1, 0, 1, 0, 0, 0])
            sage: ParallelogramPolyomino.from_dyck_word(dyck)
            [[0, 1, 0, 0, 1, 1], [1, 1, 1, 0, 0, 0]]
            sage: ParallelogramPolyomino.from_dyck_word(dyck, bijection='Delest-Viennot')
            [[0, 1, 0, 0, 1, 1], [1, 1, 1, 0, 0, 0]]
            sage: ParallelogramPolyomino.from_dyck_word(dyck, bijection='Delest-Viennot-beta')
            [[0, 0, 1, 0, 1, 1], [1, 1, 1, 0, 0, 0]]
        """
        if bijection is None or bijection == 'Delest-Viennot':
            return ParallelogramPolyomino._from_dyck_word_delest_viennot(dyck)
        if bijection == 'Delest-Viennot-beta':
            return ParallelogramPolyomino._from_dyck_word_delest_viennot_peaks_valleys(dyck)
        raise ValueError("the given bijection is not valid")

    def _to_binary_tree_Aval_Boussicault(self, position=None):
        r"""
        Convert to a binary tree using the Aval-Boussicault algorithm.

        You can use the parameter ``position`` to use the bijection on
        a new parallelogram polyomino (PP). This PP is obtained by cutting the
        PP in such a way the cell at position ``position`` becomes the
        top-left most corner of the PP.

        Reference: [ABBS2013]_

        INPUT:

        - ``position`` -- the cell position; this is a recursive parameter
          It should not be used directly

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 1, 0, 1, 1, 0, 0, 0, 1, 0]
            ....:     ]
            ....: )
            sage: pp._to_binary_tree_Aval_Boussicault()
            [[., [[., .], [[., [., .]], .]]], [[., .], .]]

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp._to_binary_tree_Aval_Boussicault()
            [., .]

            sage: pp = ParallelogramPolyomino([[1], [1]])
            sage: pp._to_binary_tree_Aval_Boussicault()
            .
        """
        from sage.combinat.binary_tree import BinaryTree
        if position is None:
            position = [0, 0]
        if self.size() == 1:
            return BinaryTree()
        result = [BinaryTree(), BinaryTree()]
        right_son = list(position)
        left_son = list(position)
        w = right_son[1] + 1
        h = left_son[0] + 1
        while w < self.width():
            if self[right_son[0]][w] == 1:
                if self[right_son[0]-1][w] == 0:
                    right_son[1] = w
                    result[1] = self._to_binary_tree_Aval_Boussicault(
                        right_son
                    )
                    break
            w += 1
        while h < self.height():
            if self[h][left_son[1]] == 1:
                if self[h][left_son[1]-1] == 0:
                    left_son[0] = h
                    result[0] = self._to_binary_tree_Aval_Boussicault(left_son)
                    break
            h += 1
        return BinaryTree(result)

    @combinatorial_map(name="To binary tree")
    def to_binary_tree(self, bijection=None):
        r"""
        Convert to a binary tree.

        INPUT:

        - ``bijection`` -- string or ``None`` (default: ``None``); the name of
          bijection to use for the conversion. The possible values are ``None``
          or ``'Aval-Boussicault'``. The ``None`` value is equivalent to
          ``'Aval-Boussicault'``.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 1, 0, 1, 1, 0, 0, 0, 1, 0]
            ....:     ]
            ....: )
            sage: pp.to_binary_tree()
            [[., [[., .], [[., [., .]], .]]], [[., .], .]]

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp.to_binary_tree()
            [., .]

            sage: pp = ParallelogramPolyomino([[1], [1]])
            sage: pp.to_binary_tree()
            .
        """
        if bijection is None or bijection == 'Aval-Boussicault':
            return self._to_binary_tree_Aval_Boussicault([0, 0])

    def _to_ordered_tree_via_dyck(self):
        r"""
        Convert the parallelogram polyominoe (PP) by using first the
        Delest-Viennot bijection between PP and Dyck paths, and then
        by using the classical bijection between Dyck paths and
        ordered trees.

        This last bijection is described in [DerZak1980]_ (see page 12 and
        Figure 3.1 of page 13).

        See :meth:`_to_dyck_delest_viennot` for the exact references.
        See also :meth:`to_ordered_tree()`.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp._to_ordered_tree_via_dyck()
            [[]]

            sage: pp = ParallelogramPolyomino([[0, 1, 1], [1, 1, 0]])
            sage: pp._to_ordered_tree_via_dyck()
            [[[]]]

            sage: pp = ParallelogramPolyomino([[0, 0, 1], [1, 0, 0]])
            sage: pp._to_ordered_tree_via_dyck()
            [[], []]

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 1, 0, 1, 1, 0, 0, 0, 1, 0]
            ....:     ]
            ....: )
            sage: pp._to_ordered_tree_via_dyck()
            [[[[]], [[[]], []]], [[]]]
        """
        return self._to_dyck_delest_viennot().to_ordered_tree()

    def _to_ordered_tree_Bou_Socci(self):
        r"""
        Return the ordered tree using the Boussicault-Socci bijection.

        This bijection is described in the article [BRS2015]_.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 1, 0, 1, 1, 0, 0, 0, 1, 0]
            ....:     ]
            ....: )
            sage: pp._to_ordered_tree_Bou_Socci()
            [[[[[]], [[[]]]]], [[]]]
            sage: pp.to_ordered_tree(bijection='Boussicault-Socci')
            [[[[[]], [[[]]]]], [[]]]

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp._to_ordered_tree_Bou_Socci()
            [[]]

            sage: pp = ParallelogramPolyomino([[1], [1]])
            sage: pp._to_ordered_tree_Bou_Socci()
            []
        """
        from sage.combinat.ordered_tree import OrderedTree
        from sage.combinat.binary_tree import BinaryTree

        def make_tree(b_tree, d):
            r"""
            This is a technical function that converts binary tree to ordered
            tree with the following construction.

            Add a virtual root v such that the root become:

            - the left son of v if ``d`` is equal to 0;

            - the right son of v if ``d`` is equal to 1;

            Then now the vertices of the ordered tree are the vertices of
            the binary tree and the virtual root.

            The edges are defined as follow:

            - if v1 is a left (resp. right) son of v2 and v2 is a right
              (resp. left) son of v3, then, in the ordered tree, v2 is the
              father of v1;
            - if v1 is a left (resp. right) son of v2 and v2 is a left (resp.
              right) son of v3, then, in the ordered tree, v2 and v3 are
              brother and v2 are on the left of v3.

            For example, if d is equal to 1

            ::

                      1
                     / \
                    2   3
                   /
                  7
                 / \
                8   4
                     \
                      5
                       \
                        6

            becomes

            ::

                    _____v_____
                   |           |
                ___1___ ____   3
                |      |    |
                2   ___7__  8
                    |  |  |
                    4  5  6

            and if d = 0, it becomes

            ::

                _________v________
                |  |      |      |
                1  2  ____7___   8
                |     |   |   |
                3     4   5   6

            INPUT:

            - ``b_tree`` -- a binary tree

            - ``d`` -- 0 or 1

            OUTPUT:

            An ordered tree.
            """
            if b_tree == BinaryTree():
                return OrderedTree([])
            res = []
            res.append(make_tree(b_tree[1 - d], 1 - d))
            res += make_tree(b_tree[d], d)
            return OrderedTree(res)
        return make_tree(
            self.to_binary_tree(bijection='Aval-Boussicault'), 1)

    @combinatorial_map(name="To ordered tree")
    def to_ordered_tree(self, bijection=None):
        r"""
        Return an ordered tree from the parallelogram polyomino.

        Different bijections can be specified.

        The bijection 'via dyck and Delest-Viennot' is the composition of
        :meth:`_to_dyck_delest_viennot` and the classical bijection between
        dyck paths and ordered trees.

        The bijection between Dyck Word and ordered trees is described
        in [DerZak1980]_ (See page 12 and 13 and Figure 3.1).

        The bijection 'Boussicault-Socci' is described in [BRS2015]_.

        INPUT:

        - ``bijection`` -- string or ``None`` (default: ``None``); the name of
          bijection to use for the conversion. The possible value are ``None``,
          ``'Boussicault-Socci'`` or ``'via dyck and Delest-Viennot'``.
          The ``None`` value is equivalent to the ``'Boussicault-Socci'``
          value.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 1, 0, 1, 1, 0, 0, 0, 1, 0]
            ....:     ]
            ....: )
            sage: pp.to_ordered_tree()
            [[[[[]], [[[]]]]], [[]]]

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp.to_ordered_tree()
            [[]]

            sage: pp = ParallelogramPolyomino([[1], [1]])
            sage: pp.to_ordered_tree()
            []

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 1, 0, 1, 1, 0, 0, 0, 1, 0]
            ....:     ]
            ....: )
            sage: pp.to_ordered_tree('via dyck and Delest-Viennot')
            [[[[]], [[[]], []]], [[]]]
        """
        if bijection is None or bijection == 'Boussicault-Socci':
            return self._to_ordered_tree_Bou_Socci()
        if bijection == 'via dyck and Delest-Viennot':
            return self._to_ordered_tree_via_dyck()

    def get_options(self):
        r"""
        Return all the options of the object.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp.get_options()
            Current options for ParallelogramPolyominoes_size
              - display:            'list'
              - drawing_components: {'bounce_0': False,
             'bounce_1': False,
             'bounce_values': False,
             'diagram': True,
             'tree': False}
              - latex:              'drawing'
              - tikz_options:       {'color_bounce_0': 'red',
             'color_bounce_1': 'blue',
             'color_line': 'black',
             'color_point': 'black',
             'line_size': 1,
             'mirror': None,
             'point_size': 3.5,
             'rotation': 0,
             'scale': 1,
             'translation': [0, 0]}
        """
        if self._options is None:
            return self.parent().get_options()
        return self._options

    def set_options(self, *get_value, **set_value):
        r"""
        Set new options to the object.
        See :class:`LocalOptions` for more info.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 0, 0, 1, 0, 1, 0, 1],
            ....:         [1, 0, 0, 0, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: pp
            [[0, 0, 0, 0, 1, 0, 1, 0, 1], [1, 0, 0, 0, 1, 1, 0, 0, 0]]
            sage: pp.set_options(display='drawing')
            sage: pp
            [1 0 0]
            [1 0 0]
            [1 0 0]
            [1 1 1]
            [0 1 1]
            [0 0 1]

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: view(PP) # not tested
            sage: pp.set_options(
            ....:     drawing_components=dict(
            ....:         diagram = True,
            ....:         bounce_0 = True,
            ....:         bounce_1 = True,
            ....:     )
            ....: )
            sage: view(PP) # not tested
        """
        if self._options is None:
            self._options = deepcopy(self.get_options())
        self._options(*get_value, **set_value)

    def upper_path(self) -> list:
        r"""
        Get the upper path of the parallelogram polyomino.

        EXAMPLES::

            sage: lower_path = [0, 0, 1, 0, 1, 1]
            sage: upper_path = [1, 1, 0, 1, 0, 0]
            sage: pp = ParallelogramPolyomino([lower_path, upper_path])
            sage: pp.upper_path()
            [1, 1, 0, 1, 0, 0]
        """
        return list(ClonableList.__getitem__(self, 1))

    def lower_path(self) -> list:
        r"""
        Get the lower path of the parallelogram polyomino.

        EXAMPLES::

            sage: lower_path = [0, 0, 1, 0, 1, 1]
            sage: upper_path = [1, 1, 0, 1, 0, 0]
            sage: pp = ParallelogramPolyomino([lower_path, upper_path])
            sage: pp.lower_path()
            [0, 0, 1, 0, 1, 1]
        """
        return list(ClonableList.__getitem__(self, 0))

    @staticmethod
    def _prefix_lengths(word, up):
        r"""
        Convert a word to a list of lengths using the following algorithm:

        1) convert each 1-``up`` letter of the word by the number of ``up``
           located on the left in the word;
        2) remove all the ``up`` letters and return the resulting list of
           integers.

        INPUT:

        - ``word`` -- a word of 0 and 1

        - ``up`` -- 0 or 1 (a letter of the word)

        OUTPUT: list of integers

        EXAMPLES::

            sage: ParallelogramPolyomino._prefix_lengths([], 1)
            []
            sage: ParallelogramPolyomino._prefix_lengths([], 0)
            []
            sage: ParallelogramPolyomino._prefix_lengths([1,1,0,1,0,0,1], 1)
            [2, 3, 3]
            sage: ParallelogramPolyomino._prefix_lengths([1,1,0,1,0,0,1], 0)
            [0, 0, 1, 3]
        """
        res = []
        h = 0
        for e in word:
            if e == up:
                h += 1
            else:
                res.append(h)
        return res

    def upper_heights(self):
        r"""
        Return the list of heights associated to each vertical step of the
        parallelogram polyomino's upper path.

        OUTPUT: list of integers

        EXAMPLES::

            sage: ParallelogramPolyomino([[0, 1], [1, 0]]).upper_heights()
            [0]
            sage: ParallelogramPolyomino(
            ....:     [[0, 0, 1, 1, 0, 1, 1, 1], [1, 0, 1, 1, 0, 1, 1, 0]]
            ....: ).upper_heights()
            [0, 1, 1, 2, 2]
        """
        return ParallelogramPolyomino._prefix_lengths(self.upper_path(), 0)

    def lower_heights(self):
        r"""
        Return the list of heights associated to each vertical step of the
        parallelogram polyomino's lower path.

        OUTPUT: list of integers

        EXAMPLES::

            sage: ParallelogramPolyomino([[0, 1], [1, 0]]).lower_heights()
            [1]
            sage: ParallelogramPolyomino(
            ....:     [[0, 0, 1, 1, 0, 1, 1, 1], [1, 0, 1, 1, 0, 1, 1, 0]]
            ....: ).lower_heights()
            [2, 2, 3, 3, 3]
        """
        return ParallelogramPolyomino._prefix_lengths(self.lower_path(), 0)

    def upper_widths(self):
        r"""
        Return the list of widths associated to each horizontal step of the
        parallelogram polyomino's upper path.

        OUTPUT: list of integers

        EXAMPLES::

            sage: ParallelogramPolyomino([[0, 1], [1, 0]]).upper_widths()
            [1]
            sage: ParallelogramPolyomino(
            ....:     [[0, 0, 1, 1, 0, 1, 1, 1], [1, 0, 1, 1, 0, 1, 1, 0]]
            ....: ).upper_widths()
            [1, 3, 5]
        """
        return ParallelogramPolyomino._prefix_lengths(self.upper_path(), 1)

    def lower_widths(self):
        r"""
        Return the list of widths associated to each horizontal step of the
        parallelogram polyomino's lower path.

        OUTPUT: list of integers

        EXAMPLES::

            sage: ParallelogramPolyomino([[0, 1], [1, 0]]).lower_widths()
            [0]
            sage: ParallelogramPolyomino(
            ....:     [[0, 0, 1, 1, 0, 1, 1, 1], [1, 0, 1, 1, 0, 1, 1, 0]]
            ....: ).lower_widths()
            [0, 0, 2]
        """
        return ParallelogramPolyomino._prefix_lengths(self.lower_path(), 1)

    def widths(self) -> list:
        r"""
        Return a list of the widths of the parallelogram polyomino.

        Namely, the parallelogram polyomino is split row by row and the
        method returns the list containing the sizes of the rows.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 0, 1, 1, 0, 0, 1, 0, 0]
            ....:     ]
            ....: )
            sage: pp.widths()
            [1, 3, 3, 3, 2]

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp.widths()
            [1]

            sage: pp = ParallelogramPolyomino([[1], [1]])
            sage: pp.widths()
            []
        """
        widths = []
        uw = self.upper_widths()
        lw = self.lower_widths()
        for i in range(len(lw)):
            widths.append(uw[i] - lw[i])
        return widths

    def degree_convexity(self) -> int:
        r"""
        Return the degree convexity of a parallelogram polyomino.

        A convex polyomino is said to be k-convex if every pair of its cells
        can be connected by a monotone path (path with south and east steps)
        with at most k changes of direction.
        The degree of convexity of a convex polyomino P is the smallest integer
        k such that P is k-convex.

        If the parallelogram polyomino is empty, the function return -1.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 0, 1, 1, 0, 0, 1, 0, 0]
            ....:     ]
            ....: )
            sage: pp.degree_convexity()
            3

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp.degree_convexity()
            0

            sage: pp = ParallelogramPolyomino([[1], [1]])
            sage: pp.degree_convexity()
            -1
        """
        l0 = len(self.bounce_path(direction=0))
        l1 = len(self.bounce_path(direction=1))
        return min(l0, l1) - 1

    def is_flat(self) -> bool:
        r"""
        Return whether the two bounce paths join together in the rightmost cell
        of the bottom row of P.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 0, 1, 1, 0, 0, 1, 0, 0]
            ....:     ]
            ....: )
            sage: pp.is_flat()
            False

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp.is_flat()
            True

            sage: pp = ParallelogramPolyomino([[1], [1]])
            sage: pp.is_flat()
            True
        """
        l0 = len(self.bounce_path(direction=0))
        l1 = len(self.bounce_path(direction=1))
        return l0 == l1

    def is_k_directed(self, k) -> bool:
        r"""
        Return whether the Polyomino Parallelogram is k-directed.

        A convex polyomino is said to be k-convex if every pair of its cells
        can be connected by a monotone path (path with south and east steps)
        with at most k changes of direction.

        The degree of convexity of a convex polyomino P is the smallest integer
        k such that P is k-convex.

        INPUT:

        - ``k`` -- nonnegative integer

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 0, 1, 1, 0, 0, 1, 0, 0]
            ....:     ]
            ....: )
            sage: pp.is_k_directed(3)
            True
            sage: pp.is_k_directed(4)
            True
            sage: pp.is_k_directed(5)
            True
            sage: pp.is_k_directed(0)
            False
            sage: pp.is_k_directed(1)
            False
            sage: pp.is_k_directed(2)
            False

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp.is_k_directed(0)
            True
            sage: pp.is_k_directed(1)
            True

            sage: pp = ParallelogramPolyomino([[1], [1]])
            sage: pp.is_k_directed(0)
            True
            sage: pp.is_k_directed(1)
            True
        """
        return self.degree_convexity() <= k

    def heights(self) -> list:
        r"""
        Return a list of heights of the parallelogram polyomino.

        Namely, the parallelogram polyomino is split column by column and
        the method returns the list containing the sizes of the columns.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 0, 1, 1, 0, 0, 1, 0, 0]
            ....:     ]
            ....: )
            sage: pp.heights()
            [3, 3, 4, 2]

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp.heights()
            [1]

            sage: pp = ParallelogramPolyomino([[1], [1]])
            sage: pp.heights()
            [0]
        """
        uh = self.upper_heights()
        lh = self.lower_heights()
        return [a - b for a, b in zip(lh, uh)]

    def width(self):
        r"""
        Return the width of the parallelogram polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 1, 0, 0, 1, 1, 0, 1, 1, 1],
            ....:         [1, 1, 1, 0, 1, 0, 0, 1, 1, 0]
            ....:     ]
            ....: )
            sage: pp.width()
            6

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp.width()
            1

            sage: pp = ParallelogramPolyomino([[1], [1]])
            sage: pp.width()
            1
        """
        if self.size() == 1:
            return 1
        return self.upper_widths()[-1]

    def height(self):
        r"""
        Return the height of the parallelogram polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 1, 0, 0, 1, 1, 0, 1, 1, 1],
            ....:         [1, 1, 1, 0, 1, 0, 0, 1, 1, 0]
            ....:     ]
            ....: )
            sage: pp.height()
            4

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp.height()
            1

            sage: pp = ParallelogramPolyomino([[1], [1]])
            sage: pp.height()
            0
        """
        if self.size() == 1:
            return 0
        return self.lower_heights()[-1]

    def cell_is_inside(self, w, h):
        r"""
        Determine whether the cell at a given position
        is inside the parallelogram polyomino.

        INPUT:

        - ``w`` -- the x coordinate of the box position

        - ``h`` -- the y coordinate of the box position

        OUTPUT:

        Return 0 if there is no cell at the given position,
        return 1 if there is a cell.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 1, 0, 0, 1, 1, 0, 1, 1, 1],
            ....:         [1, 1, 1, 0, 1, 0, 0, 1, 1, 0]
            ....:     ]
            ....: )
            sage: pp.cell_is_inside(0, 0)
            1
            sage: pp.cell_is_inside(1, 0)
            1
            sage: pp.cell_is_inside(0, 1)
            0
            sage: pp.cell_is_inside(3, 0)
            0
            sage: pp.cell_is_inside(pp.width()-1,pp.height()-1)
            1
            sage: pp.cell_is_inside(pp.width(),pp.height()-1)
            0
            sage: pp.cell_is_inside(pp.width()-1,pp.height())
            0
        """
        lower_widths = self.lower_widths()
        widths = self.widths()

        if h >= len(widths) or h < 0:
            return 0
        if lower_widths[h] <= w and w < lower_widths[h] + widths[h]:
            return 1
        return 0

    @cached_method
    def get_array(self):
        r"""
        Return an array of 0s and 1s such that the 1s represent the boxes of
        the parallelogram polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 0, 0, 1, 0, 1, 0, 1],
            ....:         [1, 0, 0, 0, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: matrix(pp.get_array())
            [1 0 0]
            [1 0 0]
            [1 0 0]
            [1 1 1]
            [0 1 1]
            [0 0 1]

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp.get_array()
            [[1]]

            sage: pp = ParallelogramPolyomino([[1], [1]])
            sage: pp.get_array()
            []
        """
        width = self.width()
        height = self.height()
        return [
            [self.cell_is_inside(w, h) for w in range(width)]
            for h in range(height)
        ]

    class _polyomino_row:
        r"""
        This is an internal class representing a single row of
        a parallelogram polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 0, 0, 1, 0, 1, 0, 1],
            ....:         [1, 0, 0, 0, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: row = ParallelogramPolyomino._polyomino_row(pp, 4)
            sage: row
            [0, 1, 1]
        """

        def __init__(self, polyomino, row):
            r"""
            The constructor of the class.

            EXAMPLES::

                sage: pp = ParallelogramPolyomino(
                ....:     [
                ....:         [0, 0, 0, 0, 1, 0, 1, 0, 1],
                ....:         [1, 0, 0, 0, 1, 1, 0, 0, 0]
                ....:     ]
                ....: )
                sage: row = ParallelogramPolyomino._polyomino_row(pp, 4)
            """
            self.polyomino = polyomino
            self.row = row

        def __getitem__(self, column):
            r"""
            Return 0 or 1 if the is a cell inside the specific column inside the
            row.

            EXAMPLES::

                sage: pp = ParallelogramPolyomino(
                ....:     [
                ....:         [0, 0, 0, 0, 1, 0, 1, 0, 1],
                ....:         [1, 0, 0, 0, 1, 1, 0, 0, 0]
                ....:     ]
                ....: )
                sage: matrix(pp.get_array())
                [1 0 0]
                [1 0 0]
                [1 0 0]
                [1 1 1]
                [0 1 1]
                [0 0 1]

                sage: row = ParallelogramPolyomino._polyomino_row(pp, 4)
                sage: [row[-1], row[0], row[1], row[2], row[3]]
                [0, 0, 1, 1, 0]
            """
            if (self.is_inside() and
                    0 <= column and column < self.polyomino.width()):
                return self.polyomino.get_array()[self.row][column]
            return 0

        def is_inside(self) -> bool:
            r"""
            Return ``True`` if the row is inside the parallelogram polyomino,
            return ``False`` otherwise.

            EXAMPLES::

                sage: PP = ParallelogramPolyomino
                sage: pp = PP(
                ....:     [
                ....:         [0, 0, 0, 0, 1, 0, 1, 0, 1],
                ....:         [1, 0, 0, 0, 1, 1, 0, 0, 0]
                ....:     ]
                ....: )
                sage: matrix(pp.get_array())
                [1 0 0]
                [1 0 0]
                [1 0 0]
                [1 1 1]
                [0 1 1]
                [0 0 1]

                sage: [
                ....:     PP._polyomino_row(pp, i).is_inside()
                ....:     for i in [-1,0,3,5,6]
                ....: ]
                [False, True, True, True, False]
            """
            return 0 <= self.row and self.row < self.polyomino.height()

        def is_outside(self) -> bool:
            r"""
            Return ``True`` if the row is outside the parallelogram polyomino,
            return ``False`` otherwise.

            EXAMPLES::

                sage: PP = ParallelogramPolyomino
                sage: pp = PP(
                ....:     [
                ....:         [0, 0, 0, 0, 1, 0, 1, 0, 1],
                ....:         [1, 0, 0, 0, 1, 1, 0, 0, 0]
                ....:     ]
                ....: )
                sage: matrix(pp.get_array())
                [1 0 0]
                [1 0 0]
                [1 0 0]
                [1 1 1]
                [0 1 1]
                [0 0 1]

                sage: [
                ....:     PP._polyomino_row(pp, i).is_outside()
                ....:     for i in [-1,0,3,5,6]
                ....: ]
                [True, False, False, False, True]
            """
            return not self.is_inside()

        def __repr__(self) -> str:
            r"""
            Return a string representation of ``self``.

            EXAMPLES::

                sage: PP = ParallelogramPolyomino
                sage: pp = PP(
                ....:     [
                ....:         [0, 0, 0, 0, 1, 0, 1, 0, 1],
                ....:         [1, 0, 0, 0, 1, 1, 0, 0, 0]
                ....:     ]
                ....: )
                sage: pp[-1]
                The (outside) row -1 of the parallelogram
                sage: pp[0]
                [1, 0, 0]
                sage: pp[5]
                [0, 0, 1]
                sage: pp[6]
                The (outside) row 6 of the parallelogram
            """
            if self.is_outside():
                return "The (outside) row %s of the parallelogram" % (self.row)
            else:
                return str(self.polyomino.get_array()[self.row])

    def __getitem__(self, row):
        r"""
        Return the row of the parallelogram polyomino.

        The index of the row can be out of range of the height of
        the parallelogram polyomino.
        In that case, the row returned is outside the parallelogram
        polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]
            ....: )
            sage: pp[0].is_inside()
            True
            sage: pp[3].is_outside()
            True
            sage: pp[-1].is_outside()
            True
            sage: pp[0][0]
            1
            sage: pp[[0, 1]]
            1
            sage: pp[2][0]
            0
            sage: pp[-1][0]
            0
            sage: pp[[4, 4]]
            0
        """
        if isinstance(row, list):
            return self._polyomino_row(self, row[0])[row[1]]
        return self._polyomino_row(self, row)

    def bounce_path(self, direction=1):
        r"""
        Return the bounce path of the parallelogram polyomino.

        The bounce path is a path with two steps (1, 0) and (0, 1).

        If 'direction' is 1 (resp. 0), the bounce path is the path
        starting at position (h=1, w=0) (resp. (h=0, w=1)) with
        initial direction, the vector (0, 1) (resp. (1, 0)), and turning
        each time the path crosses the perimeter of the parallelogram
        polyomino.

        The path is coded by a list of integers. Each integer represents
        the size of the path between two turnings.

        You can visualize the two bounce paths by using the following
        commands.

        INPUT:

        - ``direction`` -- the initial direction of the bounce path (see above
          for the definition)

        EXAMPLES::

            sage: PP = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]
            ....: )
            sage: PP.bounce_path(direction=1)
            [2, 2, 1]
            sage: PP.bounce_path(direction=0)
            [2, 1, 1, 1]

            sage: PP = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 1, 1, 0, 0, 1, 1],
            ....:         [1, 1, 1, 0, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: PP.bounce_path(direction=1)
            [3, 1, 2, 2]
            sage: PP.bounce_path(direction=0)
            [2, 4, 2]

            sage: PP = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]
            ....: )
            sage: PP.set_options(
            ....:     drawing_components=dict(
            ....:         diagram = True
            ....:         , bounce_0 = True
            ....:         , bounce_1 = True
            ....:     )
            ....: )
            sage: view(PP) # not tested

            sage: PP = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: PP.bounce_path(direction=1)
            [1]
            sage: PP.bounce_path(direction=0)
            [1]

            sage: PP = ParallelogramPolyomino([[1], [1]])
            sage: PP.bounce_path(direction=1)
            []
            sage: PP.bounce_path(direction=0)
            []

        TESTS::

            sage: PP = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]
            ....: )
            sage: PP.bounce_path(direction=1)
            [2, 2, 1]
            sage: PP.bounce_path(direction=0)
            [2, 1, 1, 1]
        """
        result = []
        pos = [0, 0]
        pos[direction] -= 1
        old = list(pos)
        ne = list(pos)
        ne[direction] += 1
        while self[ne] == 1:
            pos[direction] += 1
            while self[pos] == 1:
                pos[direction] += 1
            pos[direction] -= 1
            result.append(pos[direction]-old[direction])
            direction = 1 - direction
            old[0], old[1] = pos
            ne[0], ne[1] = pos
            ne[direction] += 1
        return result

    def bounce(self, direction=1):
        r"""
        Return the bounce of the parallelogram polyomino.

        Let ``p`` be the bounce path of the parallelogram
        polyomino (:meth:`bounce_path`). The bounce is defined by:

            ``sum([(1+ floor(i/2))*p[i] for i in range(len(p))])``

        INPUT:

        - ``direction`` -- the initial direction of the bounce path
          (see :meth:`bounce_path` for the definition)

        EXAMPLES::

            sage: PP = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]
            ....: )
            sage: PP.bounce(direction=1)
            6
            sage: PP.bounce(direction=0)
            7

            sage: PP = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 1, 1, 0, 0, 1, 1],
            ....:         [1, 1, 1, 0, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: PP.bounce(direction=1)
            12
            sage: PP.bounce(direction=0)
            10

            sage: PP = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: PP.bounce(direction=1)
            1
            sage: PP.bounce(direction=0)
            1

            sage: PP = ParallelogramPolyomino([[1], [1]])
            sage: PP.bounce(direction=1)
            0
            sage: PP.bounce(direction=0)
            0
        """
        return sum((1 + i//2) * pi
                   for i, pi in enumerate(self.bounce_path(direction)))

    def area(self):
        r"""
        Return the area of the parallelogram polyomino. The area of a
        parallelogram polyomino is the number of cells of the parallelogram
        polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1],
            ....:         [1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0]
            ....:     ]
            ....: )
            sage: pp.area()
            13

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp.area()
            1

            sage: pp = ParallelogramPolyomino([[1], [1]])
            sage: pp.area()
            0
        """
        return sum(self.heights())

    def _repr_(self) -> str:
        r"""
        Return a string representation of the parallelogram polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]
            ....: )
            sage: pp # indirect doctest
            [[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]
            sage: pp.set_options(display='drawing')
            sage: pp # indirect doctest
            [1 1 0]
            [1 1 0]
            [0 1 1]
        """
        return self.get_options()._dispatch(self, '_repr_', 'display')

    def _repr_list(self) -> str:
        r"""
        Return a string representation with list style.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]
            ....: )
            sage: pp._repr_() == pp._repr_list()
            True
            sage: pp._repr_list()
            '[[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]'
        """
        return ClonableList._repr_(self)

    def _repr_drawing(self) -> str:
        r"""
        Return a string representing a drawing of the parallelogram polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]
            ....: )
            sage: pp._repr_() == pp._repr_drawing()
            False
            sage: pp._repr_drawing()
            '[1 1 0]\n[1 1 0]\n[0 1 1]'
        """
        return str(matrix(self.get_array()))

    def get_tikz_options(self):
        r"""
        Return all the tikz options permitting to draw the parallelogram
        polyomino.

        See :class:`LocalOption` to have more informations about the
        modification of those options.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp.get_tikz_options()
            {'color_bounce_0': 'red',
             'color_bounce_1': 'blue',
             'color_line': 'black',
             'color_point': 'black',
             'line_size': 1,
             'mirror': None,
             'point_size': 3.5,
             'rotation': 0,
             'scale': 1,
             'translation': [0, 0]}
        """
        return self.get_options()['tikz_options']

    def _to_tikz_diagram(self):
        r"""
        Return the tikz code of the diagram representing ``self``.

        TESTS::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 0, 0, 1, 1], [1, 0, 1, 1, 0, 0, 0, 0]]
            ....: )
            sage: pp.to_tikz() == pp._to_tikz_diagram()
            True
            sage: print(pp.to_tikz())
            <BLANKLINE>
              \draw[color=black, line width=1] (0.000000, 5.000000) --
            (0.000000, 3.000000);
              \draw[color=black, line width=1] (3.000000, 4.000000) --
            (3.000000, 0.000000);
              \draw[color=black, line width=1] (0.000000, 5.000000) --
            (1.000000, 5.000000);
              \draw[color=black, line width=1] (1.000000, 0.000000) --
            (3.000000, 0.000000);
              \draw[color=black, line width=1] (1.000000, 5.000000) --
            (1.000000, 0.000000);
              \draw[color=black, line width=1] (2.000000, 4.000000) --
            (2.000000, 0.000000);
              \draw[color=black, line width=1] (0.000000, 4.000000) --
            (3.000000, 4.000000);
              \draw[color=black, line width=1] (0.000000, 3.000000) --
            (3.000000, 3.000000);
              \draw[color=black, line width=1] (1.000000, 2.000000) --
            (3.000000, 2.000000);
              \draw[color=black, line width=1] (1.000000, 1.000000) --
            (3.000000, 1.000000);
        """
        tikz_options = self.get_tikz_options()
        grid_width = self.width() + 1
        grid_height = self.height() + 1
        drawing_tool = _drawing_tool(
            tikz_options,
            XY=lambda v: [v[0], grid_height-1-v[1]]
        )
        res = ""
        if self.size() == 1:
            res += drawing_tool.draw_line([0, 0], [1, 0])
            return res
        res += drawing_tool.draw_line([0, 0], [0, self.lower_heights()[0]])
        res += drawing_tool.draw_line(
            [grid_width-1, self.upper_heights()[grid_width-2]],
            [grid_width-1, self.lower_heights()[grid_width-2]]
        )
        res += drawing_tool.draw_line([0, 0], [self.upper_widths()[0], 0])
        res += drawing_tool.draw_line(
            [self.lower_widths()[grid_height-2], grid_height-1],
            [self.upper_widths()[grid_height-2], grid_height-1]
        )
        for w in range(1, grid_width-1):
            h1 = self.upper_heights()[w-1]
            h2 = self.lower_heights()[w]
            res += drawing_tool.draw_line([w, h1], [w, h2])
        for h in range(1, grid_height-1):
            w1 = self.lower_widths()[h-1]
            w2 = self.upper_widths()[h]
            res += drawing_tool.draw_line([w1, h], [w2, h])
        return res

    def _to_tikz_bounce(self, directions=None):
        r"""
        Return the tikz code to display one or both bounces of ``self``.

        See :meth:`ParallelogramPolyomino.bounce_path` for more information
        about the bounce.

        TESTS::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 0, 1, 1, 0, 1, 1], [1, 0, 1, 1, 0, 1, 0, 0]]
            ....: )
            sage: pp.to_tikz() == pp._to_tikz_bounce()
            False
            sage: pp.set_options(drawing_components=dict(
            ....:     diagram=False, bounce_0=True)
            ....: )
            sage: pp.to_tikz() == pp._to_tikz_bounce([0])
            True
            sage: pp.set_options(
            ....:     drawing_components=dict(diagram=False, bounce_1=True)
            ....: )
            sage: pp.to_tikz() == pp._to_tikz_bounce([1])
            True
            sage: pp.set_options(
            ....:     drawing_components=dict(
            ....:         diagram=False, bounce_0= True, bounce_1=True
            ....:     )
            ....: )
            sage: pp.to_tikz() == pp._to_tikz_bounce([0,1])
            True
            sage: pp.to_tikz() == pp._to_tikz_bounce()
            True
            sage: pp.set_options(
            ....:     drawing_components=dict(diagram=True, bounce_0=True)
            ....: )
            sage: print(pp.to_tikz()) # indirect doctest
            <BLANKLINE>
              \draw[color=black, line width=1] (0.000000, 4.000000) --
            (0.000000, 1.000000);
              \draw[color=black, line width=1] (4.000000, 2.000000) --
            (4.000000, 0.000000);
              \draw[color=black, line width=1] (0.000000, 4.000000) --
            (1.000000, 4.000000);
              \draw[color=black, line width=1] (2.000000, 0.000000) --
            (4.000000, 0.000000);
              \draw[color=black, line width=1] (1.000000, 4.000000) --
            (1.000000, 1.000000);
              \draw[color=black, line width=1] (2.000000, 3.000000) --
            (2.000000, 0.000000);
              \draw[color=black, line width=1] (3.000000, 3.000000) --
            (3.000000, 0.000000);
              \draw[color=black, line width=1] (0.000000, 3.000000) --
            (3.000000, 3.000000);
              \draw[color=black, line width=1] (0.000000, 2.000000) --
            (4.000000, 2.000000);
              \draw[color=black, line width=1] (0.000000, 1.000000) --
            (4.000000, 1.000000);
              \draw[color=red, line width=2] (1.000000, 4.000000) --
            (1.000000, 1.000000);
              \draw[color=red, line width=2] (1.000000, 1.000000) --
            (4.000000, 1.000000);
              \draw[color=red, line width=2] (4.000000, 1.000000) --
            (4.000000, 0.000000);
        """
        if directions is None:
            directions = [0, 1]
        res = ""
        tikz_options = self.get_tikz_options()
        grid_height = self.height() + 1
        drawing_tool = _drawing_tool(
            tikz_options,
            XY=lambda v: [v[0], grid_height-1-v[1]]
        )

        def draw_bounce(direction, color):
            r"""
            Return the TIKZ code of the bounce path of ``self``.

            See :meth:`ParallelogramPolyomino.bounce_path` for more information
            about the bounce.
            """
            if (len(self.bounce_path(direction)) >
                    len(self.bounce_path(1 - direction))):
                increase_size_line = 1
            else:
                increase_size_line = 0
            res = ""
            bp = self.bounce_path(direction)
            pos = [0, 0]
            pos[1-direction] += 1
            old = list(pos)
            for e in bp:
                pos[direction] += e
                res += drawing_tool.draw_line(
                    [old[1], old[0]], [pos[1], pos[0]],
                    color=color,
                    size=2*tikz_options['line_size'] + increase_size_line,
                )
                old[0], old[1] = pos
                direction = 1-direction
            return res
        if len(self.bounce_path(0)) > len(self.bounce_path(1)):
            if 0 in directions:
                res += draw_bounce(0, tikz_options['color_bounce_0'])
            if 1 in directions:
                res += draw_bounce(1, tikz_options['color_bounce_1'])
        else:
            if 1 in directions:
                res += draw_bounce(1, tikz_options['color_bounce_1'])
            if 0 in directions:
                res += draw_bounce(0, tikz_options['color_bounce_0'])
        return res

    def _to_tikz_tree(self):
        r"""
        Return the tikz code to display a node inside the boxes which are
        nodes.
        See :meth:`ParallelogramPolyomino.box_is_node` for more information.

        TESTS::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 0, 1, 1, 0, 1, 1], [1, 0, 1, 1, 0, 1, 0, 0]]
            ....: )
            sage: pp.to_tikz() == pp._to_tikz_tree()
            False
            sage: pp.set_options(
            ....:     drawing_components=dict(diagram=False, tree=True)
            ....: )
            sage: pp.to_tikz() == pp._to_tikz_tree()
            True
            sage: pp.set_options(
            ....:     drawing_components=dict(diagram=True, tree=True)
            ....: )
            sage: print(pp.to_tikz())  # indirect doctest
            <BLANKLINE>
              \draw[color=black, line width=1] (0.000000, 4.000000) --
            (0.000000, 1.000000);
              \draw[color=black, line width=1] (4.000000, 2.000000) --
              (4.000000, 0.000000);
              \draw[color=black, line width=1] (0.000000, 4.000000) --
              (1.000000, 4.000000);
              \draw[color=black, line width=1] (2.000000, 0.000000) --
              (4.000000, 0.000000);
              \draw[color=black, line width=1] (1.000000, 4.000000) --
              (1.000000, 1.000000);
              \draw[color=black, line width=1] (2.000000, 3.000000) --
              (2.000000, 0.000000);
              \draw[color=black, line width=1] (3.000000, 3.000000) --
              (3.000000, 0.000000);
              \draw[color=black, line width=1] (0.000000, 3.000000) --
              (3.000000, 3.000000);
              \draw[color=black, line width=1] (0.000000, 2.000000) --
              (4.000000, 2.000000);
              \draw[color=black, line width=1] (0.000000, 1.000000) --
              (4.000000, 1.000000);
              \filldraw[color=black] (0.500000, 2.500000) circle (3.5pt);
              \filldraw[color=black] (0.500000, 1.500000) circle (3.5pt);
              \filldraw[color=black] (2.500000, 0.500000) circle (3.5pt);
              \filldraw[color=black] (1.500000, 2.500000) circle (3.5pt);
              \filldraw[color=black] (2.500000, 2.500000) circle (3.5pt);
              \filldraw[color=black] (3.500000, 1.500000) circle (3.5pt);
              \filldraw[color=black] (0.500000, 3.500000) circle (3.5pt);
        """
        res = ""
        tikz_options = self.get_tikz_options()
        if self.size() == 1:
            return res
        grid_height = self.height() + 1
        drawing_tool = _drawing_tool(
            tikz_options,
            XY=lambda v: [v[0] + .5, grid_height-1-v[1] - .5]
        )
        for node in self.get_BS_nodes():
            res += drawing_tool.draw_point([node[1], node[0]])
        res += drawing_tool.draw_point([0, 0])
        return res

    def _get_node_position_at_row(self, row):
        r"""
        Return the position of the leftmost cell in the row indexed by ``row``
        of the array obtained with ``get_array``.

        INPUT:

        - ``row`` -- the index of the row

        OUTPUT: a [row,column] position of the cell

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 0, 0, 1, 0, 1, 0, 1],
            ....:         [1, 0, 0, 0, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: matrix(pp.get_array())
            [1 0 0]
            [1 0 0]
            [1 0 0]
            [1 1 1]
            [0 1 1]
            [0 0 1]
            sage: pp._get_node_position_at_row(0)
            [0, 0]
            sage: pp._get_node_position_at_row(1)
            [1, 0]
            sage: pp._get_node_position_at_row(2)
            [2, 0]
            sage: pp._get_node_position_at_row(3)
            [3, 0]
            sage: pp._get_node_position_at_row(4)
            [4, 1]
            sage: pp._get_node_position_at_row(5)
            [5, 2]
        """
        h = row
        for w in range(self.width()):
            if self[h][w] == 1:
                return [h, w]
        return None

    def _get_node_position_at_column(self, column):
        r"""
        Return the position of the topmost cell in the column indexed by
        ``column`` of the array obtained with ``get_array``.

        INPUT:

        - ``column`` -- the index of the column

        OUTPUT: a [row,column] position of the cell

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 0, 0, 1, 1], [1, 0, 1, 1, 0, 0, 0, 0]]
            ....: )
            sage: matrix(pp.get_array())
            [1 0 0]
            [1 1 1]
            [0 1 1]
            [0 1 1]
            [0 1 1]
            sage: pp._get_node_position_at_column(0)
            [0, 0]
            sage: pp._get_node_position_at_column(1)
            [1, 1]
            sage: pp._get_node_position_at_column(2)
            [1, 2]
        """
        w = column
        for h in range(self.height()):
            if self[h][w] == 1:
                return [h, w]
        return None

    def get_node_position_from_box(self, box_position, direction, nb_crossed_nodes=None):
        r"""
        This function starts from a cell inside a parallelogram polyomino and
        a direction.

        If ``direction`` is equal to 0, the function selects the column
        associated with the y-coordinate of ``box_position`` and then returns
        the topmost cell of the column that is on the top of ``box_position``
        (the cell of ``box_position`` is included).

        If ``direction`` is equal to 1, the function selects the row
        associated with the x-coordinate of ``box_position`` and then returns
        the leftmost cell of the row that is on the left of ``box_position``.
        (the cell of ``box_position`` is included).

        This function updates the entry of ``nb_crossed_nodes``. The function
        increases the entry of ``nb_crossed_nodes`` by the number of boxes that
        is a node (see ``box_is_node``) located on the top if ``direction``
        is 0 (resp. on the left if ``direction`` is 1) of ``box_position``
        (cell at ``box_position`` is excluded).

        INPUT:

        - ``box_position`` -- the position of the statring cell

        - ``direction`` -- the direction (0 or 1)

        - ``nb_crossed_nodes`` -- ``[0]`` (default) a list containing just one
          integer

        OUTPUT: a [row,column] position of the cell

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 0, 0, 1, 1], [1, 0, 1, 1, 0, 0, 0, 0]]
            ....: )
            sage: matrix(pp.get_array())
            [1 0 0]
            [1 1 1]
            [0 1 1]
            [0 1 1]
            [0 1 1]
            sage: l = [0]
            sage: pp.get_node_position_from_box([3, 2], 0, l)
            [1, 2]
            sage: l
            [1]
            sage: l = [0]
            sage: pp.get_node_position_from_box([3, 2], 1, l)
            [3, 1]
            sage: l
            [1]
            sage: l = [0]
            sage: pp.get_node_position_from_box([1, 2], 0, l)
            [1, 2]
            sage: l
            [0]
            sage: l = [0]
            sage: pp.get_node_position_from_box([1, 2], 1, l)
            [1, 0]
            sage: l
            [2]
            sage: l = [0]
            sage: pp.get_node_position_from_box([3, 1], 0, l)
            [1, 1]
            sage: l
            [2]
            sage: l = [0]
            sage: pp.get_node_position_from_box([3, 1], 1, l)
            [3, 1]
            sage: l
            [0]
        """
        if nb_crossed_nodes is None:
            nb_crossed_nodes = [0]
        pos = list(box_position)
        if self[pos[0]][pos[1]] == 0:
            return None
        while self[pos[0]][pos[1]] != 0:
            pos[direction] -= 1
            if self.box_is_node(pos):
                nb_crossed_nodes[0] += 1
        pos[direction] += 1
        return pos

    def box_is_node(self, pos) -> bool:
        r"""
        Return ``True`` if the box contains a node in the context of the
        Aval-Boussicault bijection between parallelogram polyomino and binary
        tree.

        A box is a node if there is no cell on the top of the box in the
        same column or on the left of the box.in the same row.

        INPUT:

        - ``pos`` -- the [x,y] coordinate of the box

        OUTPUT: boolean

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 0, 0, 1, 1], [1, 1, 0, 1, 0, 0, 0, 0]]
            ....: )
            sage: pp.set_options(display='drawing')
            sage: pp
            [1 1 0]
            [1 1 1]
            [0 1 1]
            [0 1 1]
            [0 1 1]
            sage: pp.box_is_node([2,1])
            True
            sage: pp.box_is_node([2,0])
            False
            sage: pp.box_is_node([1,1])
            False
        """
        if self[pos[0]][pos[1]] == 0:
            return False
        if self[pos[0] - 1][pos[1]] == 0:
            return True
        if self[pos[0]][pos[1] - 1] == 0:
            return True
        return False

    def box_is_root(self, box) -> bool:
        r"""
        Return ``True`` if the box contains the root of the tree : it
        is the top-left box of the parallelogram polyomino.

        INPUT:

        - ``box`` -- the x,y coordinate of the cell

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 0, 0, 1, 1], [1, 1, 0, 1, 0, 0, 0, 0]]
            ....: )
            sage: pp.box_is_root([0, 0])
            True
            sage: pp.box_is_root([0, 1])
            False
        """
        return box[0] == 0 and box[1] == 0

    def _get_number_of_nodes_in_the_bounding_path(self, box, direction):
        r"""
        When we draw the bounding path from ``box`` to the top-left cell of
        ``self``, the path is crossing some cells containing some nodes
        defined by the Boussicault-Socci bijection
        (see :meth:`_to_ordered_tree_Bou_Socci`).

        This function returns a list of numbers that represent the number of
        nodes minus 1 that the path is crossing between each bounding.
        The starting box is excluded from the count of nodes.

        This function is a specialized tool for
        :meth:`_get_path_in_pair_of_tree_from_row()` and
        :meth:`_get_path_in_pair_of_tree_from_column()`
        each number is reduced by one to compute the path in the ordered tree
        of those functions.

        INPUT:

        - ``box`` -- the x,y coordinate of the starting point of the bounding
          path
        - ``direction`` -- the initial direction of the bounding path (1 or 0,
          1 for left and 0 for top)

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 0, 0, 1, 1], [1, 0, 1, 1, 0, 0, 0, 0]]
            ....: )
            sage: pp._get_number_of_nodes_in_the_bounding_path([4, 2], 1)
            [0, 0, 2, 0]
            sage: pp._get_number_of_nodes_in_the_bounding_path([3, 2], 1)
            [0, 0, 1, 0]
            sage: pp._get_number_of_nodes_in_the_bounding_path([2, 2], 1)
            [0, 0, 0, 0]
            sage: pp._get_number_of_nodes_in_the_bounding_path([1, 2], 1)
            [0, 1]
            sage: pp._get_number_of_nodes_in_the_bounding_path([4, 2], 0)
            [0, 1, 0]
            sage: pp._get_number_of_nodes_in_the_bounding_path([3, 2], 0)
            [0, 1, 0]
            sage: pp._get_number_of_nodes_in_the_bounding_path([2, 2], 0)
            [0, 1, 0]
            sage: pp._get_number_of_nodes_in_the_bounding_path([1, 2], 0)
            [0, 1, -1]

            sage: pp._get_number_of_nodes_in_the_bounding_path([4, 1], 1)
            [0, 0, 2, -1]
            sage: pp._get_number_of_nodes_in_the_bounding_path([3, 1], 1)
            [0, 0, 1, -1]
            sage: pp._get_number_of_nodes_in_the_bounding_path([2, 1], 1)
            [0, 0, 0, -1]
            sage: pp._get_number_of_nodes_in_the_bounding_path([1, 1], 1)
            [0, 0]
            sage: pp._get_number_of_nodes_in_the_bounding_path([4, 1], 0)
            [0, 0, 2]
            sage: pp._get_number_of_nodes_in_the_bounding_path([3, 1], 0)
            [0, 0, 1]
            sage: pp._get_number_of_nodes_in_the_bounding_path([2, 1], 0)
            [0, 0, 0]
            sage: pp._get_number_of_nodes_in_the_bounding_path([1, 1], 0)
            [0, 0, -1]

            sage: pp._get_number_of_nodes_in_the_bounding_path([1, 0], 1)
            [0, -1]
            sage: pp._get_number_of_nodes_in_the_bounding_path([0, 0], 1)
            []
            sage: pp._get_number_of_nodes_in_the_bounding_path([1, 0], 0)
            [0]
            sage: pp._get_number_of_nodes_in_the_bounding_path([0, 0], 0)
            []
        """
        path = []
        while not self.box_is_root(box):
            nb_sons = [0]
            box = self.get_node_position_from_box(box, direction, nb_sons)
            direction = 1 - direction
            path.append(nb_sons[0] - 1)
        path.reverse()
        return path

    def _get_path_in_pair_of_tree_from_row(self, line):
        r"""
        When we draw the bounding path from the left-most cell of ``line`` to
        the top-left cell of ``self``, the path is bounding in some cells that
        are nodes in the ordered tree of the Boussicault-Socci bijection.
        This function returns the path of the bounding path inside the ordered
        tree.

        The path in the ordered tree is encoded as a list of integers.
        The first integer represents the son choice between the sons of the
        root.
        Recursively, an integer represent the son choice between the sons of
        the current father.

        The bijection is described in the paper [BRS2015]_
        at page 7, the first (resp. second) ordered tree is obtained by
        gluing all roots of the ordered forest F_e (resp. F_s) to a virtual
        root. An example can be read, page 8, Figure 6.

        INPUT:

        - ``line`` -- the x coordinate of the line

        OUTPUT: list of integers

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 0, 0, 1, 1], [1, 0, 1, 1, 0, 0, 0, 0]]
            ....: )
            sage: pp._get_path_in_pair_of_tree_from_row(4)
            [0, 0, 2]
            sage: pp._get_path_in_pair_of_tree_from_row(3)
            [0, 0, 1]
            sage: pp._get_path_in_pair_of_tree_from_row(2)
            [0, 0, 0]
            sage: pp._get_path_in_pair_of_tree_from_row(1)
            [0]
            sage: pp._get_path_in_pair_of_tree_from_row(0)
            []

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 0, 0, 1, 1], [1, 1, 0, 1, 0, 0, 0, 0]]
            ....: )
            sage: pp._get_path_in_pair_of_tree_from_row(4)
            [0, 2]
            sage: pp._get_path_in_pair_of_tree_from_row(3)
            [0, 1]
            sage: pp._get_path_in_pair_of_tree_from_row(2)
            [0, 0]
            sage: pp._get_path_in_pair_of_tree_from_row(1)
            [0]
            sage: pp._get_path_in_pair_of_tree_from_row(0)
            []
        """
        pos = self._get_node_position_at_row(line)
        return self._get_number_of_nodes_in_the_bounding_path(pos, 0)

    def _get_path_in_pair_of_tree_from_column(self, column):
        r"""
        When we draw the bounding path from the top-most cell of ``column``
        to the top-left cell of ``self``, the path is bounding in some cells
        that are nodes in the ordered tree of the Boussicault-Socci bijection.
        This function returns the path of the bounding path inside the ordered
        tree.

        The path in the ordered tree is encoded as a list of integers.
        The first integer represents the son choice between the sons of the
        root.
        Recursively, an integer represent the son choice between the sons of
        the current father.

        The bijection is described in the paper [BRS2015]_
        at page 7, the first (resp. second) ordered tree is obtained by
        gluing all roots of the ordered forest F_e (resp. F_s) to a virtual
        root. An example can be read, page 8, Figure 6.

        INPUT:

        - ``column`` -- the y coordinate of the column

        OUTPUT: list of integers

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 0, 0, 1, 1], [1, 0, 1, 1, 0, 0, 0, 0]]
            ....: )
            sage: pp._get_path_in_pair_of_tree_from_column(2)
            [0, 1]
            sage: pp._get_path_in_pair_of_tree_from_column(1)
            [0, 0]
            sage: pp._get_path_in_pair_of_tree_from_column(0)
            []

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 0, 0, 1, 1], [1, 1, 0, 1, 0, 0, 0, 0]]
            ....: )
            sage: pp._get_path_in_pair_of_tree_from_column(2)
            [0, 0]
            sage: pp._get_path_in_pair_of_tree_from_row(1)
            [0]
            sage: pp._get_path_in_pair_of_tree_from_row(0)
            []
        """
        pos = self._get_node_position_at_column(column)
        return self._get_number_of_nodes_in_the_bounding_path(pos, 1)

    def get_BS_nodes(self):
        r"""
        Return the list of cells containing node of the left and right planar
        tree in the Boussicault-Socci bijection.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 0, 0, 1, 1], [1, 1, 0, 1, 0, 0, 0, 0]]
            ....: )
            sage: pp.set_options(display='drawing')
            sage: pp
            [1 1 0]
            [1 1 1]
            [0 1 1]
            [0 1 1]
            [0 1 1]
            sage: sorted(pp.get_BS_nodes())
            [[0, 1], [1, 0], [1, 2], [2, 1], [3, 1], [4, 1]]

        You can draw the point inside the parallelogram polyomino by typing
        (the left nodes are in blue, and the right node are in red) ::

            sage: pp.set_options(drawing_components=dict(tree=True))
            sage: view(pp) # not tested
        """
        result = []
        for h in range(1, self.height()):
            result.append(self._get_node_position_at_row(h))
        for w in range(1, self.width()):
            result.append(self._get_node_position_at_column(w))
        return result

    def get_right_BS_nodes(self):
        r"""
        Return the list of cells containing node of the right planar tree in
        the Boussicault-Socci bijection between parallelogram polyominoes
        and pair of ordered trees.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 0, 0, 1, 1], [1, 1, 0, 1, 0, 0, 0, 0]]
            ....: )
            sage: pp.set_options(display='drawing')
            sage: pp
            [1 1 0]
            [1 1 1]
            [0 1 1]
            [0 1 1]
            [0 1 1]
            sage: sorted(pp.get_right_BS_nodes())
            [[1, 0], [1, 2]]

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 0, 0, 1, 1], [1, 0, 1, 1, 0, 0, 0, 0]]
            ....: )
            sage: pp.set_options(display='drawing')
            sage: pp
            [1 0 0]
            [1 1 1]
            [0 1 1]
            [0 1 1]
            [0 1 1]
            sage: sorted(pp.get_right_BS_nodes())
            [[1, 0], [1, 1], [1, 2], [2, 1], [3, 1], [4, 1]]

        You can draw the point inside the parallelogram polyomino by typing,
        (the left nodes are in blue, and the right node are in red) ::

            sage: pp.set_options(drawing_components=dict(tree=True))
            sage: view(pp) # not tested
        """
        result = []
        for h in range(1, self.height()):
            path2 = self._get_path_in_pair_of_tree_from_row(h)
            if len(path2) % 2 == 1:
                result.append(self._get_node_position_at_row(h))
        for w in range(1, self.width()):
            path2 = self._get_path_in_pair_of_tree_from_column(w)
            if len(path2) % 2 == 0:
                result.append(self._get_node_position_at_column(w))
        return result

    def get_left_BS_nodes(self):
        r"""
        Return the list of cells containing node of the left planar tree in
        the Boussicault-Socci bijection between parallelogram polyominoes
        and pair of ordered trees.

        OUTPUT: list of [row,column] position of cells

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 0, 0, 1, 1], [1, 1, 0, 1, 0, 0, 0, 0]]
            ....: )
            sage: pp.set_options(display='drawing')
            sage: pp
            [1 1 0]
            [1 1 1]
            [0 1 1]
            [0 1 1]
            [0 1 1]
            sage: sorted(pp.get_left_BS_nodes())
            [[0, 1], [2, 1], [3, 1], [4, 1]]

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 0, 0, 1, 1], [1, 0, 1, 1, 0, 0, 0, 0]]
            ....: )
            sage: pp.set_options(display='drawing')
            sage: pp
            [1 0 0]
            [1 1 1]
            [0 1 1]
            [0 1 1]
            [0 1 1]
            sage: sorted(pp.get_left_BS_nodes())
            []

        You can draw the point inside the parallelogram polyomino by typing
        (the left nodes are in blue, and the right node are in red) ::

            sage: pp.set_options(drawing_components=dict(tree=True))
            sage: view(pp) # not tested
        """
        result = []
        for h in range(1, self.height()):
            path2 = self._get_path_in_pair_of_tree_from_row(h)
            if len(path2) % 2 == 0:
                result.append(self._get_node_position_at_row(h))
        for w in range(1, self.width()):
            path2 = self._get_path_in_pair_of_tree_from_column(w)
            if len(path2) % 2 == 1:
                result.append(self._get_node_position_at_column(w))
        return result

    def to_tikz(self):
        r"""
        Return the tikz code of the parallelogram polyomino.

        This code is the code present inside a tikz latex environment.

        We can modify the output with the options.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0,0,0,1,1,0,1,0,0,1,1,1],[1,1,1,0,0,1,1,0,0,1,0,0]]
            ....: )
            sage: print(pp.to_tikz())
            <BLANKLINE>
              \draw[color=black, line width=1] (0.000000, 6.000000) --
            (0.000000, 3.000000);
              \draw[color=black, line width=1] (6.000000, 2.000000) --
            (6.000000, 0.000000);
              \draw[color=black, line width=1] (0.000000, 6.000000) --
            (3.000000, 6.000000);
              \draw[color=black, line width=1] (3.000000, 0.000000) --
            (6.000000, 0.000000);
              \draw[color=black, line width=1] (1.000000, 6.000000) --
            (1.000000, 3.000000);
              \draw[color=black, line width=1] (2.000000, 6.000000) --
            (2.000000, 2.000000);
              \draw[color=black, line width=1] (3.000000, 6.000000) --
            (3.000000, 0.000000);
              \draw[color=black, line width=1] (4.000000, 4.000000) --
            (4.000000, 0.000000);
              \draw[color=black, line width=1] (5.000000, 4.000000) --
            (5.000000, 0.000000);
              \draw[color=black, line width=1] (0.000000, 5.000000) --
            (3.000000, 5.000000);
              \draw[color=black, line width=1] (0.000000, 4.000000) --
            (5.000000, 4.000000);
              \draw[color=black, line width=1] (0.000000, 3.000000) --
            (5.000000, 3.000000);
              \draw[color=black, line width=1] (2.000000, 2.000000) --
            (6.000000, 2.000000);
              \draw[color=black, line width=1] (3.000000, 1.000000) --
            (6.000000, 1.000000);
            sage: pp.set_options(
            ....:     drawing_components=dict(
            ....:         diagram=True,
            ....:         tree=True,
            ....:         bounce_0=True,
            ....:         bounce_1=True
            ....:     )
            ....: )
            sage: print(pp.to_tikz())
            <BLANKLINE>
              \draw[color=black, line width=1] (0.000000, 6.000000) --
            (0.000000, 3.000000);
              \draw[color=black, line width=1] (6.000000, 2.000000) --
            (6.000000, 0.000000);
              \draw[color=black, line width=1] (0.000000, 6.000000) --
            (3.000000, 6.000000);
              \draw[color=black, line width=1] (3.000000, 0.000000) --
            (6.000000, 0.000000);
              \draw[color=black, line width=1] (1.000000, 6.000000) --
            (1.000000, 3.000000);
              \draw[color=black, line width=1] (2.000000, 6.000000) --
            (2.000000, 2.000000);
              \draw[color=black, line width=1] (3.000000, 6.000000) --
            (3.000000, 0.000000);
              \draw[color=black, line width=1] (4.000000, 4.000000) --
            (4.000000, 0.000000);
              \draw[color=black, line width=1] (5.000000, 4.000000) --
            (5.000000, 0.000000);
              \draw[color=black, line width=1] (0.000000, 5.000000) --
            (3.000000, 5.000000);
              \draw[color=black, line width=1] (0.000000, 4.000000) --
            (5.000000, 4.000000);
              \draw[color=black, line width=1] (0.000000, 3.000000) --
            (5.000000, 3.000000);
              \draw[color=black, line width=1] (2.000000, 2.000000) --
            (6.000000, 2.000000);
              \draw[color=black, line width=1] (3.000000, 1.000000) --
            (6.000000, 1.000000);
              \draw[color=blue, line width=3] (0.000000, 5.000000) --
            (3.000000, 5.000000);
              \draw[color=blue, line width=3] (3.000000, 5.000000) --
            (3.000000, 2.000000);
              \draw[color=blue, line width=3] (3.000000, 2.000000) --
            (5.000000, 2.000000);
              \draw[color=blue, line width=3] (5.000000, 2.000000) --
            (5.000000, 0.000000);
              \draw[color=blue, line width=3] (5.000000, 0.000000) --
            (6.000000, 0.000000);
              \draw[color=red, line width=2] (1.000000, 6.000000) --
            (1.000000, 3.000000);
              \draw[color=red, line width=2] (1.000000, 3.000000) --
            (5.000000, 3.000000);
              \draw[color=red, line width=2] (5.000000, 3.000000) --
            (5.000000, 0.000000);
              \draw[color=red, line width=2] (5.000000, 0.000000) --
            (6.000000, 0.000000);
              \filldraw[color=black] (0.500000, 4.500000) circle (3.5pt);
              \filldraw[color=black] (0.500000, 3.500000) circle (3.5pt);
              \filldraw[color=black] (2.500000, 2.500000) circle (3.5pt);
              \filldraw[color=black] (3.500000, 1.500000) circle (3.5pt);
              \filldraw[color=black] (3.500000, 0.500000) circle (3.5pt);
              \filldraw[color=black] (1.500000, 5.500000) circle (3.5pt);
              \filldraw[color=black] (2.500000, 5.500000) circle (3.5pt);
              \filldraw[color=black] (3.500000, 3.500000) circle (3.5pt);
              \filldraw[color=black] (4.500000, 3.500000) circle (3.5pt);
              \filldraw[color=black] (5.500000, 1.500000) circle (3.5pt);
              \filldraw[color=black] (0.500000, 5.500000) circle (3.5pt);
        """
        res = ""
        drawing_components = self.get_options()['drawing_components']
        if 'diagram' in drawing_components and drawing_components["diagram"]:
            res += self._to_tikz_diagram()
        directions = []
        if 'bounce_0' in drawing_components and drawing_components["bounce_0"]:
            directions.append(0)
        if 'bounce_1' in drawing_components and drawing_components["bounce_1"]:
            directions.append(1)
        if len(directions) != 0:
            res += self._to_tikz_bounce(directions)
        if 'tree' in drawing_components and drawing_components["tree"]:
            res += self._to_tikz_tree()
        return res

    def geometry(self) -> list:
        r"""
        Return a pair [h, w] containing the height and the width of the
        parallelogram polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 1, 1, 1, 1], [1, 1, 1, 1, 0]]
            ....: )
            sage: pp.geometry()
            [1, 4]

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp.geometry()
            [1, 1]

            sage: pp = ParallelogramPolyomino([[1], [1]])
            sage: pp.geometry()
            [0, 1]
        """
        return [self.height(), self.width()]

    def _plot_diagram(self):
        r"""
        Return a plot of the diagram representing ``self``.

        TESTS::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 1, 1, 1, 1], [1, 1, 1, 1, 0]]
            ....: )
            sage: pp._plot_diagram()                                                    # needs sage.plot
            Graphics object consisting of 7 graphics primitives

            sage: pp = ParallelogramPolyomino([
            ....:     [0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
            ....:     [1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0]
            ....: ])
            sage: pp._plot_diagram()                                                    # needs sage.plot
            Graphics object consisting of 25 graphics primitives
        """
        G = Graphics()

        # Draw the inner grid
        for i,u,v in zip(range(self.height()-1), self.upper_widths()[1:], self.lower_widths()):
            G += line([(u,-i-1),(v,-i-1)],rgbcolor=(0,0,0))
        for i,u,v in zip(range(self.width()-1), self.upper_heights()[1:], self.lower_heights()):
            G += line([(i+1,-u),(i+1,-v)],rgbcolor=(0,0,0))

        # Draw the outer border
        lower_heights = [0] + self.lower_heights()
        for i in range(self.width()):
            if lower_heights[i] != lower_heights[i+1]:
                G += line([(i,-lower_heights[i]),(i,-lower_heights[i+1])],rgbcolor=(0,0,0),thickness=2)
        upper_heights = self.upper_heights() + [self.height()]
        for i in range(self.width()):
            if upper_heights[i] != upper_heights[i+1]:
                G += line([(i+1,-upper_heights[i]),(i+1,-upper_heights[i+1])],rgbcolor=(0,0,0),thickness=2)

        lower_widths = self.lower_widths() + [self.width()]
        for i in range(self.height()):
            if lower_widths[i] != lower_widths[i+1]:
                G += line([(lower_widths[i],-i-1),(lower_widths[i+1],-i-1)],rgbcolor=(0,0,0),thickness=2)
        upper_widths = [0] + self.upper_widths()
        for i in range(self.height()):
            if upper_widths[i] != upper_widths[i+1]:
                G += line([(upper_widths[i],-i),(upper_widths[i+1],-i)],rgbcolor=(0,0,0),thickness=2)

        return G

    def _plot_bounce(self, directions=None):
        r"""
        Return a plot of the bounce paths of ``self``.

        INPUT:

        - ``directions`` -- direction(s) `0` and/or `1` of the bounce paths

        TESTS::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 1, 1, 1, 1], [1, 1, 1, 1, 0]]
            ....: )
            sage: pp._plot_bounce(directions=[1])                                       # needs sage.plot
            Graphics object consisting of 1 graphics primitive

            sage: pp = ParallelogramPolyomino([
            ....:     [0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
            ....:     [1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0]
            ....: ])
            sage: pp._plot_bounce(directions=[0,1])                                     # needs sage.plot
            Graphics object consisting of 9 graphics primitives
        """
        if directions is None:
            directions = [0, 1]
        G = Graphics()
        if 0 in directions:
            a,b = (1,0)
            for bounce,u in enumerate(self.bounce_path(direction=0)):
                if bounce & 1:
                    u,v = a+u,b
                else:
                    u,v = a,b+u
                G += line([(a-.1,-b),(u-.1,-v)], rgbcolor=(1,0,0), thickness=1.5)
                a,b = u,v
        if 1 in directions:
            a,b = (0,1)
            for bounce,u in enumerate(self.bounce_path(direction=1)):
                if bounce & 1:
                    u,v = a,b+u
                else:
                    u,v = a+u,b
                G += line([(a,-b+.1),(u,-v+.1)], rgbcolor=(0,0,1), thickness=1.5)
                a,b = u,v
        return G

    def _plot_bounce_values(self, bounce=0):
        r"""
        Return a plot containing the value of bounce along the specified bounce path.

        TESTS::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 1, 1, 1, 1], [1, 1, 1, 1, 0]]
            ....: )
            sage: pp._plot_bounce_values()                                              # needs sage.plot
            Graphics object consisting of 4 graphics primitives

            sage: pp = ParallelogramPolyomino([
            ....:     [0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
            ....:     [1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0]
            ....: ])
            sage: pp._plot_bounce_values(bounce=1)                                      # needs sage.plot
            Graphics object consisting of 10 graphics primitives
        """
        G = Graphics()

        # Bounce path from the top
        if bounce == 0:
            a,b = (0,-1)
            for bounce,u in enumerate(self.bounce_path(direction=0)):
                if bounce & 1:
                    u,v = a+u,b
                else:
                    u,v = a,b+u
                for i in range(a,u+1):
                    for j in range(b,v+1):
                        if (i,j) != (a,b):
                            G += text(str(bounce//2 + 1), (i+.5,-j-.5),rgbcolor=(0,0,0))
                a,b = u,v
        #Bounce path from the left
        else:
            a,b = (-1,0)
            for bounce,u in enumerate(self.bounce_path(direction=1)):
                if bounce & 1:
                    u,v = a,b+u
                else:
                    u,v = a+u,b
                for i in range(a,u+1):
                    for j in range(b,v+1):
                        if (i,j) != (a,b):
                            G += text(str(bounce//2 + 1), (i+.5,-j-.5),rgbcolor=(0,0,0))
                a,b = u,v
        return G

    def _plot_tree(self):
        r"""
        Return a plot of the nodes of the tree.

        TESTS::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 1, 1, 1, 1], [1, 1, 1, 1, 0]]
            ....: )
            sage: pp._plot_tree()                                                       # needs sage.plot
            Graphics object consisting of 2 graphics primitives

            sage: pp = ParallelogramPolyomino([
            ....:     [0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
            ....:     [1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0]
            ....: ])
            sage: pp._plot_tree()                                                       # needs sage.plot
            Graphics object consisting of 2 graphics primitives
        """
        G = Graphics()
        G += point(points=((v+.5,-u-.5) for u,v in self.get_BS_nodes()),size=20)
        G += point([.5, -.5],size=20)
        return G

    def plot(self):
        r"""
        Return a plot of ``self``.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino([[0,1],[1,0]])
            sage: pp.plot()                                                             # needs sage.plot
            Graphics object consisting of 4 graphics primitives
            sage: pp.set_options(
            ....:     drawing_components=dict(
            ....:         diagram=True,
            ....:         bounce_0=True,
            ....:         bounce_1=True,
            ....:         bounce_values=0,
            ....:     )
            ....: )
            sage: pp.plot()                                                             # needs sage.plot
            Graphics object consisting of 7 graphics primitives
        """
        G = Graphics()

        drawing_components = self.get_options()['drawing_components']
        if 'diagram' in drawing_components and drawing_components["diagram"]:
            G += self._plot_diagram()
        directions = []
        if 'bounce_0' in drawing_components and drawing_components["bounce_0"]:
            directions.append(0)
        if 'bounce_1' in drawing_components and drawing_components["bounce_1"]:
            directions.append(1)
        if len(directions) != 0:
            G += self._plot_bounce(directions)
        if 'bounce_values' in drawing_components and drawing_components["bounce_values"] is not False:
            G += self._plot_bounce_values()
        if 'tree' in drawing_components and drawing_components["tree"]:
            G += self._plot_tree()

        G.set_aspect_ratio(1)
        G.axes(False)
        return G

    def size(self) -> int:
        r"""
        Return the size of the parallelogram polyomino.

        The size of a parallelogram polyomino is its half-perimeter.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 0, 0, 1, 0, 1, 1], [1, 0, 0, 0, 1, 1, 0, 0]]
            ....: )
            sage: pp.size()
            8

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp.size()
            2

            sage: pp = ParallelogramPolyomino([[1], [1]])
            sage: pp.size()
            1
        """
        return len(self.upper_path())

    def _latex_(self):
        r"""
        Return a LaTeX version of ``self``.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino([[0,1],[1,0]])
            sage: latex(pp)
            <BLANKLINE>
            \begin{tikzpicture}[scale=1]
            ...
            \end{tikzpicture}

        For more on the latex options, see
        :meth:`ParallelogramPolyominoes.options`.
        """
        return self.get_options()._dispatch(self, '_latex_', 'latex')

    def _latex_drawing(self):
        r"""
        Return a LaTeX version of ``self`` in a drawing style.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino([[0,1],[1,0]])
            sage: print(pp._latex_drawing())
            <BLANKLINE>
            \begin{tikzpicture}[scale=1]
            ...
            \end{tikzpicture}
        """
        latex.add_package_to_preamble_if_available("tikz")
        tikz_options = self.get_tikz_options()
        res = "\n\\begin{tikzpicture}[scale=%s]" % (tikz_options['scale'])
        res += self.to_tikz()
        res += "\n\\end{tikzpicture}"
        return res

    def _latex_list(self):
        r"""
        Return a LaTeX version of ``self`` in a list style.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino([[0,1],[1,0]])
            sage: pp._latex_list()
            '\\[[[0, 1], [1, 0]]\\]'
        """
        return "\\[%s\\]" % self._repr_list()


class ParallelogramPolyominoesFactory(SetFactory):
    r"""
    The parallelogram polyominoes factory.

    EXAMPLES::

        sage: PPS = ParallelogramPolyominoes(size=4)
        sage: PPS
        Parallelogram polyominoes of size 4

        sage: sorted(PPS)
        [[[0, 0, 0, 1], [1, 0, 0, 0]],
         [[0, 0, 1, 1], [1, 0, 1, 0]],
         [[0, 0, 1, 1], [1, 1, 0, 0]],
         [[0, 1, 0, 1], [1, 1, 0, 0]],
         [[0, 1, 1, 1], [1, 1, 1, 0]]]

        sage: PPS = ParallelogramPolyominoes()
        sage: PPS
        Parallelogram polyominoes
        sage: PPS.cardinality()
        +Infinity
    """

    def __call__(self, size=None, policy=None):
        r"""
        Return a family of parallelogram polyominoes enumerated with the
        parameter constraints.

        INPUT:

        - ``size`` -- integer (default: ``None``); the size of the parallelogram
                      polyominoes contained in the family.
                      If set to ``None``, the family returned contains all
                      the parallelogram polyominoes.

        EXAMPLES::

            sage: PPS = ParallelogramPolyominoes(size=4)
            sage: PPS
            Parallelogram polyominoes of size 4
            sage: sorted(PPS)
            [[[0, 0, 0, 1], [1, 0, 0, 0]],
             [[0, 0, 1, 1], [1, 0, 1, 0]],
             [[0, 0, 1, 1], [1, 1, 0, 0]],
             [[0, 1, 0, 1], [1, 1, 0, 0]],
             [[0, 1, 1, 1], [1, 1, 1, 0]]]

            sage: PPS = ParallelogramPolyominoes()
            sage: PPS
            Parallelogram polyominoes
            sage: PPS.cardinality()
            +Infinity

            sage: PPS = ParallelogramPolyominoes(size=None)
            sage: PPS
            Parallelogram polyominoes
            sage: PPS.cardinality()
            +Infinity
        """
        if policy is None:
            policy = self._default_policy

        if isinstance(size, (Integer, int)):
            return ParallelogramPolyominoes_size(size, policy)
        if size is None:
            return ParallelogramPolyominoes_all(policy)
        raise ValueError("invalid argument for Parallelogram Polyominoes "
                         "Factory")

    @lazy_attribute
    def _default_policy(self):
        r"""
        Return a default policy.

        EXAMPLES::

            sage: str(ParallelogramPolyominoes._default_policy) == (
            ....:    "Set factory policy for " +
            ....:    "<class 'sage.combinat.parallelogram_polyomino" +
            ....:    ".ParallelogramPolyomino'> " +
            ....:    "with parent Parallelogram polyominoes" +
            ....:    "[=Factory for parallelogram polyominoes(())]"
            ....: )
            True
        """
        return TopMostParentPolicy(self, (), ParallelogramPolyomino)

    def _repr_(self) -> str:
        r"""
        Return the string representation of the parallelogram polyominoes
        factory.

        EXAMPLES::

            sage: ParallelogramPolyominoes
            Factory for parallelogram polyominoes
        """
        return "Factory for parallelogram polyominoes"


ParallelogramPolyominoes = ParallelogramPolyominoesFactory()
ParallelogramPolyominoes.__doc__ = \
    ParallelogramPolyominoesFactory.__call__.__doc__


class ParallelogramPolyominoes_size(
    ParentWithSetFactory, UniqueRepresentation
):
    r"""
    The parallelogram polyominoes of size `n`.

    EXAMPLES::

        sage: PPS = ParallelogramPolyominoes(4)
        sage: PPS
        Parallelogram polyominoes of size 4
        sage: sorted(PPS)
        [[[0, 0, 0, 1], [1, 0, 0, 0]],
         [[0, 0, 1, 1], [1, 0, 1, 0]],
         [[0, 0, 1, 1], [1, 1, 0, 0]],
         [[0, 1, 0, 1], [1, 1, 0, 0]],
         [[0, 1, 1, 1], [1, 1, 1, 0]]]
    """

    def __init__(self, size, policy):
        r"""
        Construct a set of Parallelogram Polyominoes of a given size.

        EXAMPLES::

            sage: ParallelogramPolyominoes(4)
            Parallelogram polyominoes of size 4
        """
        self._size = size
        ParentWithSetFactory.__init__(
            self, (size, ), policy, category=FiniteEnumeratedSets()
        )

    def _repr_(self) -> str:
        r"""
        Return the string representation of the set of
        parallelogram polyominoes

        EXAMPLES::

            sage: ParallelogramPolyominoes(4)
            Parallelogram polyominoes of size 4
        """
        return "Parallelogram polyominoes of size %s" % (self._size)

    def an_element(self):
        r"""
        Return an element of a parallelogram polyomino of a given size.

        EXAMPLES::

            sage: PPS = ParallelogramPolyominoes(4)
            sage: PPS.an_element() in PPS
            True
        """
        return next(self.__iter__())

    def check_element(self, el, check):
        r"""
        Check is a given element `el` is in the set of parallelogram
        polyominoes of a fixed size.

        EXAMPLES::

            sage: PPS = ParallelogramPolyominoes(3)
            sage: ParallelogramPolyomino(                     # indirect doctest
            ....:     [[0, 1, 1], [1, 1, 0]]
            ....: ) in PPS
            True
        """
        if el.size() != self.size():
            raise ValueError(
                "the parallelogram polyomino has a wrong size: %s" % el.size())

    def cardinality(self):
        r"""
        Return the number of parallelogram polyominoes.

        The number of parallelogram polyominoes of size n is given by
        the Catalan number `c_{n-1}`.

        EXAMPLES::

            sage: ParallelogramPolyominoes(1).cardinality()
            1
            sage: ParallelogramPolyominoes(2).cardinality()
            1
            sage: ParallelogramPolyominoes(3).cardinality()
            2
            sage: ParallelogramPolyominoes(4).cardinality()
            5

            sage: all(
            ....:     ParallelogramPolyominoes(i).cardinality()
            ....:     == len(list(ParallelogramPolyominoes(i)))
            ....:     for i in range(1,7)
            ....: )
            True
        """
        return catalan_number(self.size() - 1)

    def __iter__(self):
        r"""
        Return a parallelogram polyomino generator.

        EXAMPLES::

            sage: len(list(ParallelogramPolyominoes(4))) == 5
            True
            sage: all(
            ....:     pp in ParallelogramPolyominoes()
            ....:     for pp in ParallelogramPolyominoes(4)
            ....: )
            True
        """
        from sage.combinat.dyck_word import DyckWords
        for dyck in DyckWords(self.size() - 1):
            yield ParallelogramPolyomino.from_dyck_word(dyck)

    def get_options(self):
        r"""
        Return all the options associated with all the elements of
        the set of parallelogram polyominoes with a fixed size.

        EXAMPLES::

            sage: pps = ParallelogramPolyominoes(5)
            sage: pps.get_options()
            Current options for ParallelogramPolyominoes_size
              - display:            'list'
            ...
        """
        return self.options

    def size(self):
        r"""
        Return the size of the parallelogram polyominoes generated by this
        parent.

        EXAMPLES::

            sage: ParallelogramPolyominoes(0).size()
            0
            sage: ParallelogramPolyominoes(1).size()
            1
            sage: ParallelogramPolyominoes(5).size()
            5
        """
        return self._size

    def set_options(self, *get_value, **set_value):
        r"""
        Set new options to the object.

        EXAMPLES::

            sage: PPS = ParallelogramPolyominoes(3)
            sage: PPS.set_options(
            ....:     drawing_components=dict(
            ....:         diagram = True,
            ....:         bounce_0 = True,
            ....:         bounce_1 = True,
            ....:     )
            ....: )
            sage: pp = PPS[0]
            sage: view(pp) # not tested
        """
        self.options(*get_value, **set_value)

    options = ParallelogramPolyominoesOptions
    r"""
    The options for ParallelogramPolyominoes.
    """


class ParallelogramPolyominoes_all(
    ParentWithSetFactory, DisjointUnionEnumeratedSets
):
    r"""
    This class enumerates all the parallelogram polyominoes.

    EXAMPLES::

        sage: PPS = ParallelogramPolyominoes()
        sage: PPS
        Parallelogram polyominoes
    """

    def __init__(self, policy):
        r"""
        Construct the set of all parallelogram polyominoes.

        EXAMPLES::

            sage: PPS = ParallelogramPolyominoes()
            sage: PPS
            Parallelogram polyominoes

            sage: ParallelogramPolyomino([[0, 1, 1], [1, 1, 0]])  in PPS
            True

            sage: PPS = ParallelogramPolyominoes()
            sage: next(PPS.__iter__()) in PPS
            True
        """
        ParentWithSetFactory.__init__(
            self, (), policy, category=FiniteEnumeratedSets()
        )
        DisjointUnionEnumeratedSets.__init__(
            self, Family(
                PositiveIntegers(),
                lambda n: ParallelogramPolyominoes_size(
                    n, policy=self.facade_policy()
                )
            ),
            facade=True, keepkey=False, category=self.category()
        )

    def _repr_(self) -> str:
        r"""
        Return a string representation of the set of parallelogram polyominoes.

        EXAMPLES::

            sage: PPS = ParallelogramPolyominoes()
            sage: PPS
            Parallelogram polyominoes
        """
        return "Parallelogram polyominoes"

    def check_element(self, el, check):
        r"""
        Check is a given element `el` is in the set of parallelogram
        polyominoes.

        EXAMPLES::

            sage: PPS = ParallelogramPolyominoes()
            sage: ParallelogramPolyomino(                     # indirect doctest
            ....:     [[0, 1, 1], [1, 1, 0]]
            ....: ) in PPS
            True
        """
        pass

    def get_options(self):
        r"""
        Return all the options associated with the set of
        parallelogram polyominoes.

        EXAMPLES::

            sage: PPS = ParallelogramPolyominoes()
            sage: options = PPS.get_options()
            sage: options
            Current options for ParallelogramPolyominoes_size
              - display:            'list'
            ...
        """
        return self.options

    def set_options(self, *get_value, **set_value):
        r"""
        Set new options to the object.

        EXAMPLES::

            sage: PPS = ParallelogramPolyominoes()
            sage: PPS.set_options(
            ....:     drawing_components=dict(
            ....:         diagram = True,
            ....:         bounce_0 = True,
            ....:         bounce_1 = True,
            ....:     )
            ....: )
            sage: pp = next(iter(PPS))
            sage: view(pp) # not tested
        """
        self.options(*get_value, **set_value)

    options = ParallelogramPolyominoesOptions
    r"""
    The options for ParallelogramPolyominoes.
    """
