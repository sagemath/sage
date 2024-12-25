r"""
Triangulations of a point configuration

A point configuration is a finite set of points in Euclidean space or,
more generally, in projective space. A triangulation is a simplicial
decomposition of the convex hull of a given point configuration such
that all vertices of the simplices end up lying on points of the
configuration. That is, there are no new vertices apart from the
initial points.

Note that points that are not vertices of the convex hull need not be
used in the triangulation. A triangulation that does make use of all
points of the configuration is called fine, and you can restrict
yourself to such triangulations if you want. See
:class:`PointConfiguration` and
:meth:`~PointConfiguration.restrict_to_fine_triangulations` for
more details.

Finding a single triangulation and listing all connected
triangulations is implemented natively in this package. However, for
more advanced options [TOPCOM]_ needs to be installed; see :ref:`spkg_topcom`.

.. NOTE::

    TOPCOM and the internal algorithms tend to enumerate
    triangulations in a different order. This is why we always
    explicitly specify the engine as ``engine='topcom'`` or
    ``engine='internal'`` in the doctests. In your own applications,
    you do not need to specify the engine. By default, TOPCOM is used
    if it is available and the internal algorithms are used otherwise.

EXAMPLES:

First, we select the internal implementation for enumerating
triangulations::

    sage: PointConfiguration.set_engine('internal')   # to make doctests independent of TOPCOM

A 2-dimensional point configuration::

    sage: p = PointConfiguration([[0,0], [0,1], [1,0], [1,1], [-1,-1]]); p
    A point configuration in affine 2-space over Integer Ring consisting
    of 5 points. The triangulations of this point configuration are
    assumed to be connected, not necessarily fine, not necessarily regular.

.. PLOT::
    :width: 300 px

    p = PointConfiguration([[-1,-1], [1,1], [1,0], [0,1], [0,0]])
    sphinx_plot(p.plot(axes=False))

A triangulation of it::

    sage: t = p.triangulate(); t  # a single triangulation
    (<1,3,4>, <2,3,4>)
    sage: len(t)
    2
    sage: t[0]
    (1, 3, 4)
    sage: t[1]
    (2, 3, 4)
    sage: list(t)
    [(1, 3, 4), (2, 3, 4)]
    sage: t.plot(axes=False)                                                            # needs sage.plot
    Graphics object consisting of 12 graphics primitives

.. PLOT::
    :width: 300 px

    p = PointConfiguration([[0,0], [0,1], [1,0], [1,1], [-1,-1]])
    t = p.triangulate()
    sphinx_plot(t.plot(axes=False))

List triangulations of it::

    sage: list(p.triangulations())
    [(<1,3,4>, <2,3,4>),
     (<0,1,3>, <0,1,4>, <0,2,3>, <0,2,4>),
     (<1,2,3>, <1,2,4>),
     (<0,1,2>, <0,1,4>, <0,2,4>, <1,2,3>)]
    sage: p_fine = p.restrict_to_fine_triangulations(); p_fine
    A point configuration in affine 2-space over Integer Ring consisting
    of 5 points. The triangulations of this point configuration are
    assumed to be connected, fine, not necessarily regular.
    sage: list(p_fine.triangulations())
    [(<0,1,3>, <0,1,4>, <0,2,3>, <0,2,4>),
     (<0,1,2>, <0,1,4>, <0,2,4>, <1,2,3>)]

A 3-dimensional point configuration::

    sage: p = [[0,-1,-1], [0,0,1], [0,1,0], [1,-1,-1], [1,0,1], [1,1,0]]
    sage: points = PointConfiguration(p)
    sage: triang = points.triangulate()
    sage: triang.plot(axes=False)                                                       # needs sage.plot
    Graphics3d Object

.. PLOT::
    :width: 300 px

    p = [[0,-1,-1], [0,0,1], [0,1,0], [1,-1,-1], [1,0,1], [1,1,0]]
    points = PointConfiguration(p)
    triang = points.triangulate()
    sphinx_plot(triang.plot(axes=False))

The standard example of a non-regular triangulation (requires TOPCOM)::

    sage: # optional - topcom
    sage: PointConfiguration.set_engine('topcom')
    sage: p = PointConfiguration([[-1,-5/9], [0,10/9], [1,-5/9],
    ....:                         [-2,-10/9], [0,20/9], [2,-10/9]])
    sage: p_regular = p.restrict_to_regular_triangulations(True)
    sage: regular = p_regular.triangulations_list()
    sage: p_nonregular = p.restrict_to_regular_triangulations(False)
    sage: nonregular = p_nonregular.triangulations_list()
    sage: len(regular)
    16
    sage: len(nonregular)
    2
    sage: nonregular[0].plot(aspect_ratio=1, axes=False)                                # needs sage.plot
    Graphics object consisting of 25 graphics primitives
    sage: PointConfiguration.set_engine('internal')   # to make doctests independent of TOPCOM

Note that the points need not be in general position. That is, the
points may lie in a hyperplane and the linear dependencies will be
removed before passing the data to TOPCOM which cannot handle it::

    sage: points = [[0,0,0,1], [0,3,0,1], [3,0,0,1], [0,0,1,1],
    ....:           [0,3,1,1], [3,0,1,1], [1,1,2,1]]
    sage: points = [p + [1,2,3] for p in points]
    sage: pc = PointConfiguration(points)
    sage: pc.ambient_dim()
    7
    sage: pc.dim()
    3
    sage: pc.triangulate()
    (<0,1,2,6>, <0,1,3,6>, <0,2,3,6>, <1,2,4,6>, <1,3,4,6>, <2,3,5,6>, <2,4,5,6>)
    sage: _ in pc.triangulations()
    True
    sage: len(pc.triangulations_list())
    26

AUTHORS:

- Volker Braun: initial version, 2010

- Josh Whitney: added functionality for computing
  volumes and secondary polytopes of PointConfigurations

- Marshall Hampton: improved documentation and doctest coverage

- Volker Braun: rewrite using Parent/Element and categories. Added
  a Point class. More doctests. Less zombies.

- Volker Braun: Cythonized parts of it, added a C++ implementation
  of the bistellar flip algorithm to enumerate all connected
  triangulations.

- Volker Braun 2011: switched the triangulate() method to the
  placing triangulation (faster).
"""

########################################################################
# Note: The doctests that make use of TOPCOM are
#       marked # optional - topcom
#       If you have it installed, run doctests as
#
#   sage -tp 4 --long --optional=sage,topcom sage/geometry/triangulation/
########################################################################


########################################################################
#       Copyright (C) 2010 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2010 Josh Whitney <josh.r.whitney@gmail.com>
#       Copyright (C) 2010 Marshall Hampton <hamptonio@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
########################################################################

import itertools

from sage.features import FeatureNotPresentError
from sage.features.topcom import TOPCOMExecutable
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector

from copy import copy
import sys
import pexpect


from sage.geometry.triangulation.base import \
    PointConfiguration_base, Point, ConnectedTriangulationsIterator

from sage.geometry.triangulation.element import Triangulation


########################################################################
class PointConfiguration(UniqueRepresentation, PointConfiguration_base):
    """
    A collection of points in Euclidean (or projective) space.

    This is the parent class for the triangulations of the point
    configuration. There are a few options to specifically select what
    kind of triangulations are admissible.

    INPUT:

    The constructor accepts the following arguments:

    - ``points`` -- the points; technically, any iterable of iterables
      will do. In particular, a :class:`PointConfiguration` can be passed.

    - ``projective`` -- boolean (default: ``False``); whether the
      point coordinates should be interpreted as projective (``True``)
      or affine (``False``) coordinates. If necessary, points are
      projectivized by setting the last homogeneous coordinate to one
      and/or affine patches are chosen internally.

    - ``connected`` -- boolean (default: ``True``); whether the
      triangulations should be connected to the regular triangulations
      via bistellar flips. These are much easier to compute than all
      triangulations.

    - ``fine`` -- boolean (default: ``False``); whether the
      triangulations must be fine, that is, make use of all points of
      the configuration

    - ``regular`` -- boolean or ``None`` (default: ``None``); whether
      the triangulations must be regular. A regular triangulation is
      one that is induced by a piecewise-linear convex support
      function. In other words, the shadows of the faces of a
      polyhedron in one higher dimension.

      * ``True``: Only regular triangulations.

      * ``False``: Only non-regular triangulations.

      * ``None`` (default): Both kinds of triangulation.

    - ``star`` -- either ``None`` or a point; whether the
      triangulations must be star. A triangulation is star if all
      maximal simplices contain a common point. The central point can
      be specified by its index (an integer) in the given points or by
      its coordinates (anything iterable.)

    EXAMPLES::

        sage: p = PointConfiguration([[0,0], [0,1], [1,0], [1,1], [-1,-1]]); p
        A point configuration in affine 2-space over Integer Ring
        consisting of 5 points. The triangulations of this point
        configuration are assumed to be connected, not necessarily fine,
        not necessarily regular.
        sage: p.triangulate()  # a single triangulation
        (<1,3,4>, <2,3,4>)
    """

    # we cache the output of _have_TOPCOM() in this class variable
    _have_TOPCOM_cached = None

    # whether to use TOPCOM. Will be set to True or False during
    # initialization. All implementations should check this boolean
    # variable to decide whether to call TOPCOM or not
    _use_TOPCOM = None

    @classmethod
    def _have_TOPCOM(cls):
        r"""
        Return whether TOPCOM is installed.

        EXAMPLES::

            sage: PointConfiguration._have_TOPCOM()    # optional - topcom
            True
        """
        if PointConfiguration._have_TOPCOM_cached is not None:
            return PointConfiguration._have_TOPCOM_cached

        try:
            out = next(PointConfiguration._TOPCOM_exec('points2placingtriang',
                                                  '[[0,1],[1,1]]', verbose=False))
            PointConfiguration._have_TOPCOM_cached = True
            assert out == '{{0,1}}',\
                'TOPCOM ran but did not produce the correct output!'
        except (FeatureNotPresentError, pexpect.ExceptionPexpect):
            PointConfiguration._have_TOPCOM_cached = False

        PointConfiguration.set_engine('auto')
        return PointConfiguration._have_TOPCOM_cached

    @staticmethod
    def __classcall__(cls, points, projective=False, connected=True, fine=False, regular=None, star=None):
        r"""
        Normalize the constructor arguments to be unique keys.

        EXAMPLES::

            sage: pc1 = PointConfiguration([[1,2], [2,3], [3,4]], connected=True)
            sage: pc2 = PointConfiguration(((1,2), (2,3), (3,4)), regular=None)
            sage: pc1 is pc2   # indirect doctest
            True
        """
        if isinstance(points, PointConfiguration_base):
            pc = points
            points = tuple( p.projective() for p in points )
            projective = True
            defined_affine = pc.is_affine()
        elif projective:
            points = tuple( tuple(p) for p in points )
            defined_affine = False
        else:
            points = tuple( tuple(p)+(1,) for p in points )
            defined_affine = True
        if star is not None and star not in ZZ:
            star_point = tuple(star)
            if len(star_point) < len(points[0]):
                star_point = tuple(star)+(1,)
            star = points.index(star_point)
        return super().__classcall__(cls, points, connected, fine,
                                     regular, star, defined_affine)

    def __init__(self, points, connected, fine, regular, star, defined_affine):
        """
        Initialize a :class:`PointConfiguration` object.

        EXAMPLES::

            sage: p = PointConfiguration([[0,4], [2,3], [3,2], [4,0], [3,-2], [2,-3],
            ....:                         [0,-4], [-2,-3], [-3,-2], [-4,0], [-3,2], [-2,3]])
            sage: len(p.triangulations_list())  # long time (26s on sage.math, 2012)
            16796

        TESTS::

            sage: TestSuite(p).run()
        """
        # first, test if we have TOPCOM and set up class variables accordingly
        PointConfiguration._have_TOPCOM()

        assert connected in [True, False], 'Unknown value: connected=' + str(connected)
        self._connected = connected
        if not connected and not PointConfiguration._have_TOPCOM():
            raise ValueError('You must install TOPCOM to find non-connected triangulations.')

        assert fine in [True, False], 'Unknown value: fine=' + str(fine)
        self._fine = fine

        assert regular in [True, False, None], 'Unknown value: regular=' + str(regular)
        self._regular = regular
        if regular is not None and not PointConfiguration._have_TOPCOM():
            raise ValueError('You must install TOPCOM to test for regularity.')

        assert star is None or star in ZZ, 'Unknown value: fine=' + str(star)
        self._star = star

        PointConfiguration_base.__init__(self, points, defined_affine)

    @classmethod
    def set_engine(cls, engine='auto'):
        r"""
        Set the engine used to compute triangulations.

        INPUT:

        - ``engine`` -- either ``'auto'`` (default), ``'internal'``, or
          ``'topcom'``. The latter two instruct this package to always use
          its own triangulation algorithms or TOPCOM's algorithms,
          respectively. By default (``'auto'``), internal routines are used.

        EXAMPLES::

            sage: # optional - topcom
            sage: p = PointConfiguration([[0,0], [0,1], [1,0], [1,1], [-1,-1]])
            sage: p.set_engine('internal')   # to make doctests independent of TOPCOM
            sage: p.triangulate()
            (<1,3,4>, <2,3,4>)
            sage: p.set_engine('topcom')
            sage: p.triangulate()
            (<0,1,2>, <0,1,4>, <0,2,4>, <1,2,3>)
            sage: p.set_engine('internal')
        """
        engine = engine.lower()
        if engine not in ['auto', 'topcom', 'internal']:
            raise ValueError('Unknown value for "engine": '+str(engine))

        PointConfiguration._use_TOPCOM = (engine == 'topcom')

    def star_center(self):
        r"""
        Return the center used for star triangulations.

        .. SEEALSO:: :meth:`restrict_to_star_triangulations`.

        OUTPUT:

        A :class:`~sage.geometry.triangulation.base.Point` if a
        distinguished star central point has been fixed.
        :exc:`ValueError` exception is raised otherwise.

        EXAMPLES::

            sage: pc = PointConfiguration([(1,0), (-1,0), (0,1), (0,2)], star=(0,1)); pc
            A point configuration in affine 2-space over Integer Ring
            consisting of 4 points. The triangulations of this point
            configuration are assumed to be connected, not necessarily
            fine, not necessarily regular, and star with center P(0, 1).
            sage: pc.star_center()
            P(0, 1)

            sage: pc_nostar = pc.restrict_to_star_triangulations(None); pc_nostar
            A point configuration in affine 2-space over Integer Ring
            consisting of 4 points. The triangulations of this point
            configuration are assumed to be connected, not necessarily
            fine, not necessarily regular.
            sage: pc_nostar.star_center()
            Traceback (most recent call last):
            ...
            ValueError: The point configuration has no star center defined.
        """
        if self._star is None:
            raise ValueError('The point configuration has no star center defined.')
        else:
            return self[self._star]

    def __reduce__(self):
        r"""
        Override __reduce__ to correctly pickle/unpickle.

        TESTS::

            sage: p = PointConfiguration([[0, 1], [0, 0], [1, 0], [1,1]])
            sage: loads(p.dumps()) is p
            True

            sage: p = PointConfiguration([[0, 1, 1], [0, 0, 1], [1, 0, 1], [1,1, 1]],
            ....:                        projective=True)
            sage: loads(p.dumps()) is p
            True
        """
        if self.is_affine():
            points = tuple( p.affine() for p in self )
            return (PointConfiguration, (points, False,
                                         self._connected, self._fine, self._regular, self._star))
        else:
            points = tuple( p.projective() for p in self )
            return (PointConfiguration, (points, True,
                                         self._connected, self._fine, self._regular, self._star))

    def an_element(self):
        """
        Synonymous for :meth:`triangulate`.

        TESTS::

            sage: p = PointConfiguration([[0, 1], [0, 0], [1, 0], [1,1]])
            sage: p.an_element()
            (<0,1,3>, <1,2,3>)
        """
        return self.triangulate()

    def _element_constructor_(self, e):
        """
        Construct a triangulation.

        TESTS::

            sage: p = PointConfiguration([[0, 1], [0, 0], [1, 0], [1,1]])
            sage: p._element_constructor_([(0,1,2), (2,3,0)])
            (<0,1,2>, <0,2,3>)
        """
        return self.element_class(e, parent=self)

    Element = Triangulation

    def __iter__(self):
        """
        Iterate through the points of the point configuration.

        OUTPUT:

        Returns projective coordinates of the points. See also the
        ``PointConfiguration.points()`` method, which returns affine
        coordinates.

        EXAMPLES::

            sage: p = PointConfiguration([[1,1], [2,2], [3,3]])
            sage: list(p)     # indirect doctest
            [P(1, 1), P(2, 2), P(3, 3)]
            sage: [p[i] for i in range(p.n_points())]
            [P(1, 1), P(2, 2), P(3, 3)]
            sage: list(p.points())
            [P(1, 1), P(2, 2), P(3, 3)]
            sage: [p.point(i) for i in range(p.n_points())]
            [P(1, 1), P(2, 2), P(3, 3)]
        """
        yield from self.points()

    def _repr_(self):
        r"""
        Return a string representation.

        TESTS::

            sage: p = PointConfiguration([[1,1,1], [-1,1,1], [1,-1,1], [-1,-1,1], [1,1,-1],
            ....:                         [-1,1,-1], [1,-1,-1], [-1,-1,-1], [0,0,0]])
            sage: p._repr_()
            'A point configuration in affine 3-space over Integer Ring
            consisting of 9 points. The triangulations of this point
            configuration are assumed to be connected, not necessarily
            fine, not necessarily regular.'

            sage: PointConfiguration([[1, 1, 1], [-1, 1, 1], [1, -1, 1], [-1, -1, 1]],
            ....:                    projective=True)
            A point configuration in projective 2-space over Integer
            Ring consisting of 4 points. The triangulations of this
            point configuration are assumed to be connected,
            not necessarily fine, not necessarily regular.
        """
        s = 'A point configuration in'
        if self.is_affine():
            s += ' affine'
        else:
            s += ' projective'
        s += " %s-space over %s" % (self.ambient_dim(),self.base_ring())
        if len(self) == 1:
            s += ' consisting of '+str(len(self))+' point. '
        else:
            s += ' consisting of '+str(len(self))+' points. '

        s += 'The triangulations of this point configuration are assumed to be'

        if self._connected:
            s += ' connected,'
        else:
            s += ' not necessarily connected,'

        if self._fine:
            s += ' fine,'
        else:
            s += ' not necessarily fine,'

        if self._regular:
            s += ' regular'
        elif self._regular is False: # may be False or None, with different meanings
            s += ' irregular'
        else:
            s += ' not necessarily regular'

        if self._star is None:
            s += '.'
        else:
            s += ', and star with center '+str(self.star_center())+'.'
        if self.n_points() == 0:
            s = 'The pointless empty configuration'
        return s

    def _TOPCOM_points(self):
        r"""
        Convert the list of input points to a string that can be fed
        to TOPCOM.

        TESTS::

            sage: p = PointConfiguration([[1,1,1], [-1,1,1], [1,-1,1], [-1,-1,1], [1,1,-1]])
            sage: p._TOPCOM_points()
            '[[0,0,0,1],[-2,0,0,1],[0,-2,0,1],[-2,-2,0,1],[0,0,-2,1]]'
        """
        s = '['
        s += ','.join([
                '[' + ','.join(map(str,p.reduced_projective())) + ']'
                for p in self ])
        s += ']'
        return s

    @classmethod
    def _TOPCOM_exec(cls, executable, input_string, verbose=True):
        r"""
        Run TOPCOM.

        INPUT:

        - ``executable`` -- string; the name of the executable

        - ``input_string`` -- string; will be piped into the running
          executable's stdin

        - ``verbose`` -- boolean; whether to print out the TOPCOM
          interaction

        TESTS::

            sage: p = PointConfiguration([[1,1,1], [-1,1,1], [1,-1,1], [-1,-1,1], [1,1,-1]])
            sage: out = p._TOPCOM_exec('points2placingtriang',
            ....:                      '[[0,0,0,1],[-2,0,0,1],[0,-2,0,1],[-2,-2,0,1],[0,0,-2,1]]',
            ....:                      verbose=True)
            sage: list(out)       # optional - topcom
            #### TOPCOM input ####
            # points2placingtriang
            # [[0,0,0,1],[-2,0,0,1],[0,-2,0,1],[-2,-2,0,1],[0,0,-2,1]]
            #### TOPCOM output ####
            # {{0,1,2,4},{1,2,3,4}}
            #######################
            ['{{0,1,2,4},{1,2,3,4}}']
        """
        timeout = 600
        executable_name, *args = executable.split()
        executable_absname = TOPCOMExecutable(executable_name).absolute_filename()
        proc = pexpect.spawn(executable_absname, args, timeout=timeout)
        proc.expect(r'Evaluating Commandline Options \.\.\.')
        proc.expect(r'\.\.\. done\.')
        proc.setecho(0)
        assert proc.readline().strip() == b''

        if verbose:
            print("#### TOPCOM input ####")
            print("# " + executable)
            print("# " + input_string)
            sys.stdout.flush()

        proc.send(input_string)
        proc.send('X\nX\n')

        if verbose:
            print("#### TOPCOM output ####")
            sys.stdout.flush()

        while True:
            try:
                line = proc.readline().strip()
                if not isinstance(line, str):
                    line = line.decode()
            except pexpect.TIMEOUT:
                if verbose:
                    print('# Still running ' + str(executable))
                continue
            if len(line) == 0:  # EOF
                break
            if verbose:
                print("# " + line)
                sys.stdout.flush()

            try:
                yield line.strip()
            except GeneratorExit:
                proc.close(force=True)
                return

        if verbose:
            print("#######################")
            sys.stdout.flush()

    def _TOPCOM_communicate(self, executable, verbose=True):
        r"""
        Execute TOPCOM and parse the output into a
        :class:`~sage.geometry.triangulation.element.Triangulation`.

        TESTS::

            sage: p = PointConfiguration([[1,1,1], [-1,1,1], [1,-1,1], [-1,-1,1], [1,1,-1]])
            sage: out = p._TOPCOM_communicate('points2placingtriang', verbose=True)
            sage: list(out)       # optional - topcom
            #### TOPCOM input ####
            # points2placingtriang
            # [[0,0,0,1],[-2,0,0,1],[0,-2,0,1],[-2,-2,0,1],[0,0,-2,1]]
            #### TOPCOM output ####
            # {{0,1,2,4},{1,2,3,4}}
            #######################
            [(<0,1,2,4>, <1,2,3,4>)]
        """
        for line in self._TOPCOM_exec(executable,
                                      self._TOPCOM_points(), verbose):
            triangulation = line[ line.find('{{')+2 : line.rfind('}}') ]
            triangulation = triangulation.split('},{')
            triangulation = [ [ QQ(t) for t in triangle.split(',') ]
                              for triangle in triangulation ]

            if self._star is not None:
                o = self._star
                if not all( t.count(o) > 0 for t in triangulation):
                    continue

            yield self(triangulation)

    def _TOPCOM_triangulations(self, verbose=True):
        r"""
        Return all triangulations satisfying the restrictions imposed.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0], [0,1], [1,0], [1,1], [-1,-1]])
            sage: iter = p._TOPCOM_triangulations(verbose=True)
            sage: next(iter)     # optional - topcom
            #### TOPCOM input ####
            # points2triangs
            # [[0,0,1],[0,1,1],[1,0,1],[1,1,1],[-1,-1,1]]
            #### TOPCOM output ####
            # T[0] := {{0,1,2},{0,1,4},{0,2,4},{1,2,3}};
            (<0,1,2>, <0,1,4>, <0,2,4>, <1,2,3>)
        """
        command = 'points2'

        if not self._connected:
            command += 'all'

        if self._fine:
            command += 'fine'

        command += 'triangs'

        if self._regular:
            command += ' --regular'
        if self._regular is False:
            command += ' --nonregular'

        yield from self._TOPCOM_communicate(command, verbose)

    def _TOPCOM_triangulate(self, verbose=True):
        r"""
        Return one (in no particular order) triangulation subject
        to all restrictions imposed previously.

        INPUT:

        - ``verbose`` -- boolean; whether to print out the TOPCOM
          interaction

        OUTPUT:

        A :class:`~sage.geometry.triangulation.element.Triangulation`
        satisfying all restrictions imposed. This raises a :exc:`ValueError`
        if no such triangulation exists.

        EXAMPLES::

            sage: # optional - topcom
            sage: p = PointConfiguration([[0,0], [0,1], [1,0], [1,1], [-1,-1]])
            sage: p.set_engine('topcom')
            sage: p._TOPCOM_triangulate(verbose=False)
            (<0,1,2>, <0,1,4>, <0,2,4>, <1,2,3>)
            sage: list( p.triangulate() )
            [(0, 1, 2), (0, 1, 4), (0, 2, 4), (1, 2, 3)]
            sage: p.set_engine('internal')
        """
        assert self._regular is not False, \
            'When asked for a single triangulation TOPCOM ' + \
            'always returns a regular triangulation.'

        command = "points2"
        if self._fine:
            command += "finetriang"
        else:
            command += "placingtriang"

        return next(self._TOPCOM_communicate(command, verbose))

    def restrict_to_regular_triangulations(self, regular=True):
        """
        Restrict to regular triangulations.

        NOTE:

        Regularity testing requires the optional TOPCOM package.

        INPUT:

        - ``regular`` -- ``True``, ``False``, or ``None``; whether to
          restrict to regular triangulations, irregular
          triangulations, or lift any restrictions on regularity

        OUTPUT:

        A new :class:`PointConfiguration` with the same points, but
        whose triangulations will all be regular as specified. See
        :class:`PointConfiguration` for details.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0], [0,1], [1,0], [1,1], [-1,-1]]); p
            A point configuration in affine 2-space over Integer Ring
            consisting of 5 points. The triangulations of this point
            configuration are assumed to be connected, not necessarily
            fine, not necessarily regular.
            sage: len(p.triangulations_list())
            4
            sage: PointConfiguration.set_engine('topcom')
            sage: p_regular = p.restrict_to_regular_triangulations() # optional - topcom
            sage: len(p_regular.triangulations_list())               # optional - topcom
            4
            sage: p == p_regular.restrict_to_regular_triangulations(regular=None) # optional - topcom
            True
            sage: PointConfiguration.set_engine('internal')
        """
        return PointConfiguration(self,
                                  connected=self._connected,
                                  fine=self._fine,
                                  regular=regular,
                                  star=self._star)

    def restrict_to_connected_triangulations(self, connected=True):
        """
        Restrict to connected triangulations.

        NOTE:

        Finding non-connected triangulations requires the optional
        TOPCOM package.

        INPUT:

        - ``connected`` -- boolean; whether to restrict to
          triangulations that are connected by bistellar flips to the
          regular triangulations

        OUTPUT:

        A new :class:`PointConfiguration` with the same points, but
        whose triangulations will all be in the connected
        component. See :class:`PointConfiguration` for details.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0], [0,1], [1,0], [1,1], [-1,-1]]); p
            A point configuration in affine 2-space over Integer Ring
            consisting of 5 points. The triangulations of this point
            configuration are assumed to be connected, not necessarily
            fine, not necessarily regular.
            sage: len(p.triangulations_list())
            4
            sage: PointConfiguration.set_engine('topcom')
            sage: p_all = p.restrict_to_connected_triangulations(connected=False)  # optional - topcom
            sage: len(p_all.triangulations_list())                                 # optional - topcom
            4
            sage: p == p_all.restrict_to_connected_triangulations(connected=True)  # optional - topcom
            True
            sage: PointConfiguration.set_engine('internal')
        """
        return PointConfiguration(self,
                                  connected=connected,
                                  fine=self._fine,
                                  regular=self._regular,
                                  star=self._star)

    def restrict_to_fine_triangulations(self, fine=True):
        """
        Restrict to fine triangulations.

        INPUT:

        - ``fine`` -- boolean; whether to restrict to fine triangulations

        OUTPUT:

        A new :class:`PointConfiguration` with the same points, but
        whose triangulations will all be fine. See
        :class:`PointConfiguration` for details.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0], [0,1], [1,0], [1,1], [-1,-1]])
            sage: p
            A point configuration in affine 2-space over Integer Ring
            consisting of 5 points. The triangulations of this point
            configuration are assumed to be connected, not necessarily
            fine, not necessarily regular.

            sage: len(p.triangulations_list())
            4
            sage: p_fine = p.restrict_to_fine_triangulations()
            sage: len(p.triangulations_list())
            4
            sage: p == p_fine.restrict_to_fine_triangulations(fine=False)
            True
        """
        return PointConfiguration(self,
                                  connected=self._connected,
                                  fine=fine,
                                  regular=self._regular,
                                  star=self._star)

    def restrict_to_star_triangulations(self, star):
        """
        Restrict to star triangulations with the given point as the
        center.

        INPUT:

        - ``origin`` -- ``None`` or an integer or the coordinates of a
          point. An integer denotes the index of the central point. If
          ``None`` is passed, any restriction on the starshape will be
          removed.

        OUTPUT:

        A new :class:`PointConfiguration` with the same points, but
        whose triangulations will all be star. See
        :class:`PointConfiguration` for details.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0], [0,1], [1,0], [1,1], [-1,-1]])
            sage: len(list(p.triangulations()))
            4
            sage: p_star =  p.restrict_to_star_triangulations(0)
            sage: p_star is p.restrict_to_star_triangulations((0,0))
            True
            sage: p_star.triangulations_list()
            [(<0,1,3>, <0,1,4>, <0,2,3>, <0,2,4>)]
            sage: p_newstar = p_star.restrict_to_star_triangulations(1)  # pick different origin
            sage: p_newstar.triangulations_list()
            [(<1,2,3>, <1,2,4>)]
            sage: p == p_star.restrict_to_star_triangulations(star=None)
            True
        """
        return PointConfiguration(self,
                                  connected=self._connected,
                                  fine=self._fine,
                                  regular=self._regular,
                                  star=star)

    def triangulations(self, verbose=False):
        r"""
        Return all triangulations.

        - ``verbose`` -- boolean (default: ``False``); whether to
          print out the TOPCOM interaction, if any

        OUTPUT:

        A generator for the triangulations satisfying all the
        restrictions imposed. Each triangulation is returned as a
        :class:`~sage.geometry.triangulation.element.Triangulation` object.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0], [0,1], [1,0], [1,1], [-1,-1]])
            sage: iter = p.triangulations()
            sage: next(iter)
            (<1,3,4>, <2,3,4>)
            sage: next(iter)
            (<0,1,3>, <0,1,4>, <0,2,3>, <0,2,4>)
            sage: next(iter)
            (<1,2,3>, <1,2,4>)
            sage: next(iter)
            (<0,1,2>, <0,1,4>, <0,2,4>, <1,2,3>)
            sage: p.triangulations_list()
            [(<1,3,4>, <2,3,4>),
             (<0,1,3>, <0,1,4>, <0,2,3>, <0,2,4>),
             (<1,2,3>, <1,2,4>),
             (<0,1,2>, <0,1,4>, <0,2,4>, <1,2,3>)]
            sage: p_fine = p.restrict_to_fine_triangulations()
            sage: p_fine.triangulations_list()
            [(<0,1,3>, <0,1,4>, <0,2,3>, <0,2,4>),
             (<0,1,2>, <0,1,4>, <0,2,4>, <1,2,3>)]

         Note that we explicitly asked the internal algorithm to
         compute the triangulations. Using TOPCOM, we obtain the same
         triangulations but in a different order::

            sage: # optional - topcom
            sage: p.set_engine('topcom')
            sage: iter = p.triangulations()
            sage: next(iter)
            (<0,1,2>, <0,1,4>, <0,2,4>, <1,2,3>)
            sage: next(iter)
            (<0,1,3>, <0,1,4>, <0,2,3>, <0,2,4>)
            sage: next(iter)
            (<1,2,3>, <1,2,4>)
            sage: next(iter)
            (<1,3,4>, <2,3,4>)
            sage: p.triangulations_list()
            [(<0,1,2>, <0,1,4>, <0,2,4>, <1,2,3>),
             (<0,1,3>, <0,1,4>, <0,2,3>, <0,2,4>),
             (<1,2,3>, <1,2,4>),
             (<1,3,4>, <2,3,4>)]
            sage: p_fine = p.restrict_to_fine_triangulations()
            sage: p_fine.set_engine('topcom')
            sage: p_fine.triangulations_list()
            [(<0,1,2>, <0,1,4>, <0,2,4>, <1,2,3>),
             (<0,1,3>, <0,1,4>, <0,2,3>, <0,2,4>)]
            sage: p.set_engine('internal')
        """
        if self._use_TOPCOM:
            yield from self._TOPCOM_triangulations(verbose)
        else:
            if not self._connected:
                raise ValueError('Need TOPCOM to find disconnected triangulations.')
            if (self._regular is not None):
                raise ValueError('Need TOPCOM to test for regularity.')
            ci = ConnectedTriangulationsIterator(self, star=self._star, fine=self._fine)
            for encoded_triangulation in ci:
                yield self(encoded_triangulation)

    def triangulations_list(self, verbose=False):
        r"""
        Return all triangulations.

        INPUT:

        - ``verbose`` -- boolean; whether to print out the TOPCOM
          interaction, if any

        OUTPUT:

        A list of triangulations (see
        :class:`~sage.geometry.triangulation.element.Triangulation`)
        satisfying all restrictions imposed previously.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0], [0,1], [1,0], [1,1]])
            sage: p.triangulations_list()
            [(<0,1,2>, <1,2,3>), (<0,1,3>, <0,2,3>)]
            sage: list(map(list, p.triangulations_list()))
            [[(0, 1, 2), (1, 2, 3)], [(0, 1, 3), (0, 2, 3)]]
            sage: p.set_engine('topcom')
            sage: p.triangulations_list()      # optional - topcom
            [(<0,1,2>, <1,2,3>), (<0,1,3>, <0,2,3>)]
            sage: p.set_engine('internal')
        """
        return list(self.triangulations(verbose))

    def triangulate(self, verbose=False):
        r"""
        Return one (in no particular order) triangulation.

        INPUT:

        - ``verbose`` -- boolean; whether to print out the TOPCOM
          interaction, if any

        OUTPUT:

        A :class:`~sage.geometry.triangulation.element.Triangulation`
        satisfying all restrictions imposed. This raises a :exc:`ValueError`
        if no such triangulation exists.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0], [0,1], [1,0], [1,1], [-1,-1]])
            sage: p.triangulate()
            (<1,3,4>, <2,3,4>)
            sage: list( p.triangulate() )
            [(1, 3, 4), (2, 3, 4)]

        Using TOPCOM yields a different, but equally good, triangulation::

            sage: # optional - topcom
            sage: p.set_engine('topcom')
            sage: p.triangulate()
            (<0,1,2>, <0,1,4>, <0,2,4>, <1,2,3>)
            sage: list(p.triangulate())
            [(0, 1, 2), (0, 1, 4), (0, 2, 4), (1, 2, 3)]
            sage: p.set_engine('internal')
        """
        if self._use_TOPCOM and self._regular is not False:
            try:
                return self._TOPCOM_triangulate(verbose)
            except StopIteration:
                # either topcom did not return a triangulation or we filtered it out
                pass

        if self._connected and not self._fine and self._regular is not False and self._star is None:
            return self.placing_triangulation()

        try:
            return next(self.triangulations(verbose))
        except StopIteration:
            # there is no triangulation
            pass
        raise ValueError('No triangulation with the required properties.')

    def convex_hull(self):
        """
        Return the convex hull of the point configuration.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0], [0,1], [1,0], [1,1], [-1,-1]])
            sage: p.convex_hull()
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 4 vertices
        """
        try:
            return self._polyhedron
        except AttributeError:
            pass

        from sage.geometry.polyhedron.constructor import Polyhedron
        pts = [p.reduced_affine() for p in self.points()]
        self._polyhedron = Polyhedron(vertices=pts)
        return self._polyhedron

    @cached_method
    def restricted_automorphism_group(self):
        r"""
        Return the restricted automorphism group.

        First, let the linear automorphism group be the subgroup of
        the affine group `AGL(d,\RR) = GL(d,\RR) \ltimes \RR^d`
        preserving the `d`-dimensional point configuration. The
        affine group acts in the usual way `\vec{x}\mapsto
        A\vec{x}+b` on the ambient space.

        The restricted automorphism group is the subgroup of the
        linear automorphism group generated by permutations of
        points. See [BSS2009]_ for more details and a description of the
        algorithm.

        OUTPUT:

        A
        :class:`PermutationGroup<sage.groups.perm_gps.permgroup.PermutationGroup_generic>`
        that is isomorphic to the restricted automorphism group is
        returned.

        Note that in Sage, permutation groups always act on positive
        integers while lists etc. are indexed by nonnegative
        integers. The indexing of the permutation group is chosen to
        be shifted by ``+1``. That is, the transposition ``(i,j)`` in
        the permutation group corresponds to exchange of ``self[i-1]``
        and ``self[j-1]``.

        EXAMPLES::

            sage: pyramid = PointConfiguration([[1,0,0], [0,1,1], [0,1,-1],
            ....:                               [0,-1,-1], [0,-1,1]])
            sage: G = pyramid.restricted_automorphism_group()                           # needs sage.graphs sage.groups
            sage: G == PermutationGroup([[(3,5)], [(2,3),(4,5)], [(2,4)]])              # needs sage.graphs sage.groups
            True
            sage: DihedralGroup(4).is_isomorphic(G)                                     # needs sage.graphs sage.groups
            True

        The square with an off-center point in the middle. Note that
        the middle point breaks the restricted automorphism group
        `D_4` of the convex hull::

            sage: square = PointConfiguration([(3/4,3/4), (1,1), (1,-1), (-1,-1), (-1,1)])
            sage: square.restricted_automorphism_group()                                # needs sage.graphs sage.groups
            Permutation Group with generators [(3,5)]
            sage: DihedralGroup(1).is_isomorphic(_)                                     # needs sage.graphs sage.groups
            True
        """
        v_list = [ vector(p.projective()) for p in self ]
        Qinv = sum( v.column() * v.row() for v in v_list ).inverse()

        # construct the graph
        from sage.graphs.graph import Graph

        # Was set to sparse = False, but there is a problem with Graph
        # backends. It should probably be set back to sparse = False as soon as
        # the backends are fixed.
        G = Graph(sparse=True)
        for i in range(len(v_list)):
            for j in range(i+1,len(v_list)):
                v_i = v_list[i]
                v_j = v_list[j]
                G.add_edge(i+1,j+1, v_i * Qinv * v_j)

        return G.automorphism_group(edge_labels=True)

    def face_codimension(self, point):
        r"""
        Return the smallest `d\in\ZZ` such that ``point`` is
        contained in the interior of a codimension-`d` face.

        EXAMPLES::

            sage: triangle = PointConfiguration([[0,0], [1,-1], [1,0], [1,1]])
            sage: triangle.point(2)
            P(1, 0)
            sage: triangle.face_codimension(2)
            1
            sage: triangle.face_codimension([1,0])
            1

        This also works for degenerate cases like the tip of the
        pyramid over a square (which saturates four inequalities)::

            sage: pyramid = PointConfiguration([[1,0,0], [0,1,1], [0,1,-1],
            ....:                               [0,-1,-1], [0,-1,1]])
            sage: pyramid.face_codimension(0)
            3
        """
        try:
            p = vector(self.point(point).reduced_affine())
        except TypeError:
            p = vector(point)

        inequalities = []
        for ieq in self.convex_hull().inequality_generator():
            if (ieq.A()*p + ieq.b() == 0):
                inequalities += [ ieq.vector() ]
        return matrix(inequalities).rank()

    def face_interior(self, dim=None, codim=None):
        """
        Return points by the codimension of the containing face in the convex hull.

        EXAMPLES::

            sage: triangle = PointConfiguration([[-1,0], [0,0], [1,-1], [1,0], [1,1]])
            sage: triangle.face_interior()
            ((1,), (3,), (0, 2, 4))
            sage: triangle.face_interior(dim=0)    # the vertices of the convex hull
            (0, 2, 4)
            sage: triangle.face_interior(codim=1)  # interior of facets
            (3,)
        """
        assert not (dim is not None and codim is not None), "You cannot specify both dim and codim."

        if (dim is not None):
            return self.face_interior()[self.convex_hull().dim()-dim]
        if (codim is not None):
            return self.face_interior()[codim]

        try:
            return self._face_interior
        except AttributeError:
            pass

        d = [ self.face_codimension(i) for i in range(self.n_points()) ]

        return tuple( tuple(i for i in range(self.n_points()) if d[i] == codim )
                      for codim in range(self.dim()+1) )

    def exclude_points(self, point_idx_list):
        """
        Return a new point configuration with the given points
        removed.

        INPUT:

        - ``point_idx_list`` -- list of integers; the indices of
          points to exclude

        OUTPUT:

        A new :class:`PointConfiguration` with the given points
        removed.

        EXAMPLES::

            sage: p = PointConfiguration([[-1,0], [0,0], [1,-1], [1,0], [1,1]])
            sage: list(p)
            [P(-1, 0), P(0, 0), P(1, -1), P(1, 0), P(1, 1)]
            sage: q = p.exclude_points([3])
            sage: list(q)
            [P(-1, 0), P(0, 0), P(1, -1), P(1, 1)]
            sage: p.exclude_points(p.face_interior(codim=1)).points()
            (P(-1, 0), P(0, 0), P(1, -1), P(1, 1))
        """
        points = [self.point(i) for i in range(self.n_points())
                  if i not in point_idx_list]
        return PointConfiguration(points,
                                  projective=False,
                                  connected=self._connected,
                                  fine=self._fine,
                                  regular=self._regular,
                                  star=self._star)

    def volume(self, simplex=None):
        """
        Find `n!` times the `n`-volume of a simplex of dimension `n`.

        INPUT:

        - ``simplex`` -- (optional argument) a simplex from a
          triangulation T specified as a list of point indices

        OUTPUT:

        * If a simplex was passed as an argument: `n!` * (volume of ``simplex``).

        * Without argument: `n!` * (the total volume of the convex hull).

        EXAMPLES:

        The volume of the standard simplex should always be 1::

            sage: p = PointConfiguration([[0,0], [1,0], [0,1], [1,1]])
            sage: p.volume([0,1,2])
            1
            sage: simplex = p.triangulate()[0]  # first simplex of triangulation
            sage: p.volume(simplex)
            1

        The square can be triangulated into two minimal simplices, so
        in the "integral" normalization its volume equals two::

            sage: p.volume()
            2

        .. NOTE::

            We return `n!` * (metric volume of the simplex) to ensure that
            the volume is an integer.  Essentially, this normalizes
            things so that the volume of the standard `n`-simplex is 1.
            See [GKZ1994]_ page 182.
        """
        if (simplex is None):
            return sum([ self.volume(s) for s in self.triangulate() ])

        #Form a matrix whose columns are the points of simplex
        #with the first point of simplex shifted to the origin.
        v = [ self.point(i).reduced_affine_vector() for i in simplex ]
        m = matrix([ v_i - v[0] for v_i in v[1:] ])
        return abs(m.det())

    def secondary_polytope(self):
        r"""
        Calculate the secondary polytope of the point configuration.

        For a definition of the secondary polytope, see [GKZ1994]_ page 220
        Definition 1.6.

        Note that if you restricted the admissible triangulations of
        the point configuration then the output will be the
        corresponding face of the whole secondary polytope.

        OUTPUT:

        The secondary polytope of the point configuration as an
        instance of
        :class:`~sage.geometry.polyhedron.base.Polyhedron_base`.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0], [1,0], [2,1], [1,2], [0,1]])
            sage: poly = p.secondary_polytope()
            sage: poly.vertices_matrix()
            [1 1 3 3 5]
            [3 5 1 4 1]
            [4 2 5 2 4]
            [2 4 2 5 4]
            [5 3 4 1 1]
            sage: poly.Vrepresentation()
            (A vertex at (1, 3, 4, 2, 5),
             A vertex at (1, 5, 2, 4, 3),
             A vertex at (3, 1, 5, 2, 4),
             A vertex at (3, 4, 2, 5, 1),
             A vertex at (5, 1, 4, 4, 1))
            sage: poly.Hrepresentation()
            (An equation (0, 0, 1, 2, 1) x - 13 == 0,
             An equation (1, 0, 0, 2, 2) x - 15 == 0,
             An equation (0, 1, 0, -3, -2) x + 13 == 0,
             An inequality (0, 0, 0, -1, -1) x + 7 >= 0,
             An inequality (0, 0, 0, 1, 0) x - 2 >= 0,
             An inequality (0, 0, 0, -2, -1) x + 11 >= 0,
             An inequality (0, 0, 0, 0, 1) x - 1 >= 0,
             An inequality (0, 0, 0, 3, 2) x - 14 >= 0)
        """
        from sage.geometry.polyhedron.constructor import Polyhedron
        #TODO: once restriction to regular triangulations is fixed,
        #change the next line to only take the regular triangulations,
        #since they are the vertices of the secondary polytope anyway.
        l = self.triangulations_list()
        return Polyhedron(vertices=[x.gkz_phi() for x in l])

    def circuits_support(self):
        r"""
        A generator for the supports of the circuits of the point configuration.

        See :meth:`circuits` for details.

        OUTPUT:

        A generator for the supports `C_-\cup C_+` (returned as a
        Python tuple) for all circuits of the point configuration.

        EXAMPLES::

            sage: p = PointConfiguration([(0,0), (+1,0), (-1,0), (0,+1), (0,-1)])
            sage: sorted(p.circuits_support())
            [(0, 1, 2), (0, 3, 4), (1, 2, 3, 4)]
        """
        n = len(self)
        U = [self[i].reduced_projective() for i in range(n)]

        # the index set of U
        I = set(range(n))
        # The (indices of) known independent elements of U
        independent_k = [(i,) for i in range(n)]
        supports_k = []

        for k in range(2, self.dim() + 3):

            # possibly linear dependent subsets
            supports_knext = set()
            possible_dependency = set()
            for indep in independent_k:
                indep_plus_one = [ tuple(sorted(indep+(i,))) for i in (I-set(indep)) ]
                possible_dependency.update(indep_plus_one)
            for supp in supports_k:
                supp_plus_one = [ tuple(sorted(supp+(i,))) for i in (I-set(supp)) ]
                possible_dependency.difference_update(supp_plus_one)
                supports_knext.update(supp_plus_one)

            # remember supports and independents for the next k-iteration
            supports_k = list(supports_knext)
            independent_k = []
            for idx in possible_dependency:
                rk = matrix([ U[i] for i in idx ]).rank()
                if rk == k:
                    independent_k.append(idx)
                else:
                    supports_k.append(idx)
                    yield idx
        assert independent_k == []  # there are no independent (self.dim()+3)-tuples

    def circuits(self):
        r"""
        Return the circuits of the point configuration.

        Roughly, a circuit is a minimal linearly dependent subset of
        the points. That is, a circuit is a partition

        .. MATH::

            \{ 0, 1, \dots, n-1 \} = C_+ \cup C_0 \cup C_-

        such that there is an (unique up to an overall normalization) affine
        relation

        .. MATH::

            \sum_{i\in C_+}  \alpha_i \vec{p}_i =
            \sum_{j\in C_-}  \alpha_j \vec{p}_j

        with all positive (or all negative) coefficients, where
        `\vec{p}_i=(p_1,\dots,p_k,1)` are the projective coordinates
        of the `i`-th point.

        OUTPUT:

        The list of (unsigned) circuits as triples `(C_+, C_0,
        C_-)`. The swapped circuit `(C_-, C_0, C_+)` is not returned
        separately.

        EXAMPLES::

            sage: p = PointConfiguration([(0,0), (+1,0), (-1,0), (0,+1), (0,-1)])
            sage: sorted(p.circuits())
            [((0,), (1, 2), (3, 4)), ((0,), (3, 4), (1, 2)), ((1, 2), (0,), (3, 4))]


        TESTS::

            sage: U=matrix([
            ....:    [ 0, 0, 0, 0, 0, 2, 4,-1, 1, 1, 0, 0, 1, 0],
            ....:    [ 0, 0, 0, 1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0],
            ....:    [ 0, 2, 0, 0, 0, 0,-1, 0, 1, 0, 1, 0, 0, 1],
            ....:    [ 0, 1, 1, 0, 0, 1, 0,-2, 1, 0, 0,-1, 1, 1],
            ....:    [ 0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0]
            ....: ])
            sage: p = PointConfiguration(U.columns())
            sage: len(p.circuits())    # long time
            218
        """
        try:
            return self._circuits
        except AttributeError:
            pass

        n = len(self)
        U = [self[i].reduced_projective() for i in range(n)]

        Circuits = ()
        for support in self.circuits_support():
            m = matrix([ U[i] for i in support ]).transpose()
            ker = m.right_kernel().basis()[0]
            assert len(ker) == len(support)
            Cplus = [ support[i] for i in range(len(support)) if ker[i] > 0 ]
            Cminus = [ support[i] for i in range(len(support)) if ker[i] < 0 ]
            Czero = set( range(n) ).difference(support)
            Circuits += ( (tuple(Cplus), tuple(Czero), tuple(Cminus)), )
        self._circuits = Circuits
        return Circuits

    def positive_circuits(self, *negative):
        r"""
        Return the positive part of circuits with fixed negative part.

        A circuit is a pair `(C_+, C_-)`, each consisting of a subset
        (actually, an ordered tuple) of point indices.

        INPUT:

        - ``*negative`` -- integer; the indices of points

        OUTPUT: a tuple of all circuits with `C_-` = ``negative``

        EXAMPLES::

            sage: p = PointConfiguration([(1,0,0), (0,1,0), (0,0,1), (-2,0,-1), (-2,-1,0),
            ....:                         (-3,-1,-1), (1,1,1), (-1,0,0), (0,0,0)])
            sage: sorted(p.positive_circuits(8))
            [(0, 1, 2, 5), (0, 1, 4), (0, 2, 3), (0, 3, 4, 6), (0, 5, 6), (0, 7)]
            sage: p.positive_circuits(0,5,6)
            ((8,),)
        """
        pos = ()
        negative = tuple(sorted(negative))
        for circuit in self.circuits():
            Cpos = circuit[0]
            Cneg = circuit[2]
            if Cpos == negative:
                pos += ( Cneg, )
            elif Cneg == negative:
                pos += ( Cpos, )
        return pos

    def bistellar_flips(self):
        r"""
        Return the bistellar flips.

        OUTPUT:

        The bistellar flips as a tuple. Each flip is a pair
        `(T_+,T_-)` where `T_+` and `T_-` are partial triangulations
        of the point configuration.

        EXAMPLES::

            sage: pc = PointConfiguration([(0,0),(1,0),(0,1),(1,1)])
            sage: pc.bistellar_flips()
            (((<0,1,3>, <0,2,3>), (<0,1,2>, <1,2,3>)),)
            sage: Tpos, Tneg = pc.bistellar_flips()[0]
            sage: Tpos.plot(axes=False)                                                 # needs sage.plot
            Graphics object consisting of 11 graphics primitives
            sage: Tneg.plot(axes=False)                                                 # needs sage.plot
            Graphics object consisting of 11 graphics primitives

        The 3d analog::

            sage: pc = PointConfiguration([(0,0,0),(0,2,0),(0,0,2),(-1,0,0),(1,1,1)])
            sage: pc.bistellar_flips()
            (((<0,1,2,3>, <0,1,2,4>), (<0,1,3,4>, <0,2,3,4>, <1,2,3,4>)),)

        A 2d flip on the base of the pyramid over a square::

            sage: pc = PointConfiguration([(0,0,0),(0,2,0),(0,0,2),(0,2,2),(1,1,1)])
            sage: pc.bistellar_flips()
            (((<0,1,3>, <0,2,3>), (<0,1,2>, <1,2,3>)),)
            sage: Tpos, Tneg = pc.bistellar_flips()[0]
            sage: Tpos.plot(axes=False)                                                 # needs sage.plot
            Graphics3d Object
        """
        flips = []
        for C in self.circuits():
            Cpos = list(C[0])
            Cneg = list(C[2])
            Tpos = [Cpos + Cneg[0:i] + Cneg[i+1:len(Cneg)]
                    for i in range(len(Cneg))]
            Tneg = [Cneg + Cpos[0:i] + Cpos[i+1:len(Cpos)]
                    for i in range(len(Cpos))]
            flips.append((self.element_class(Tpos, parent=self, check=False),
                          self.element_class(Tneg, parent=self, check=False)))
        return tuple(flips)

    def lexicographic_triangulation(self):
        r"""
        Return the lexicographic triangulation.

        The algorithm was taken from [PUNTOS]_.

        EXAMPLES::

            sage: p = PointConfiguration([(0,0), (+1,0), (-1,0), (0,+1), (0,-1)])
            sage: p.lexicographic_triangulation()
            (<1,3,4>, <2,3,4>)

        TESTS::

            sage: U = matrix([
            ....:    [ 0, 0, 0, 0, 0, 2, 4,-1, 1, 1, 0, 0, 1, 0],
            ....:    [ 0, 0, 0, 1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0],
            ....:    [ 0, 2, 0, 0, 0, 0,-1, 0, 1, 0, 1, 0, 0, 1],
            ....:    [ 0, 1, 1, 0, 0, 1, 0,-2, 1, 0, 0,-1, 1, 1],
            ....:    [ 0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0]
            ....: ])
            sage: pc = PointConfiguration(U.columns())
            sage: pc.lexicographic_triangulation()
            (<1,3,4,7,10,13>, <1,3,4,8,10,13>, <1,3,6,7,10,13>, <1,3,6,8,10,13>,
             <1,4,6,7,10,13>, <1,4,6,8,10,13>, <2,3,4,6,7,12>, <2,3,4,7,12,13>,
             <2,3,6,7,12,13>, <2,4,6,7,12,13>, <3,4,5,6,9,12>, <3,4,5,8,9,12>,
             <3,4,6,7,11,12>, <3,4,6,9,11,12>, <3,4,7,10,11,13>, <3,4,7,11,12,13>,
             <3,4,8,9,10,12>, <3,4,8,10,12,13>, <3,4,9,10,11,12>, <3,4,10,11,12,13>,
             <3,5,6,8,9,12>, <3,6,7,10,11,13>, <3,6,7,11,12,13>, <3,6,8,9,10,12>,
             <3,6,8,10,12,13>, <3,6,9,10,11,12>, <3,6,10,11,12,13>, <4,5,6,8,9,12>,
             <4,6,7,10,11,13>, <4,6,7,11,12,13>, <4,6,8,9,10,12>, <4,6,8,10,12,13>,
             <4,6,9,10,11,12>, <4,6,10,11,12,13>)
            sage: len(_)
            34
        """
        lex_supp = set()
        for circuit in self.circuits():
            Cplus = circuit[0]
            Cminus = circuit[2]
            s0 = min(Cplus + Cminus)
            if s0 in Cplus:
                lex_supp.add(Cplus)
            else:
                lex_supp.add(Cminus)

        lex_supp = sorted(lex_supp, key=lambda x:-len(x))
        basepts = copy(lex_supp)
        for i in range(len(lex_supp)-1):
            for j in range(i+1,len(lex_supp)):
                if set(lex_supp[j]).issubset(set(lex_supp[i])):
                    try:
                        basepts.remove(lex_supp[i])
                    except ValueError:
                        pass

        basepts = [ (len(b),)+b for b in basepts ]     # decorate
        basepts = sorted(basepts)                      # sort
        basepts = [ b[1:] for b in basepts ]           # undecorate

        def make_cotriang(basepts):
            if len(basepts) == 0:
                return [frozenset()]
            triangulation = set()
            for tail in make_cotriang(basepts[1:]):
                for head in basepts[0]:
                    triangulation.update([ frozenset([head]).union(tail) ])

            nonminimal = set()
            for rel in itertools.combinations(triangulation, 2):
                if rel[0].issubset(rel[1]):
                    nonminimal.update([rel[1]])
                if rel[1].issubset(rel[0]):
                    nonminimal.update([rel[0]])
            triangulation.difference_update(nonminimal)

            triangulation = [ [len(t)]+sorted(t) for t in triangulation ] # decorate
            triangulation = sorted(triangulation)                         # sort
            triangulation = [ frozenset(t[1:]) for t in triangulation ]   # undecorate

            return triangulation

        triangulation = make_cotriang(basepts)
        I = frozenset(range(self.n_points()))
        triangulation = [ tuple(I.difference(t)) for t in triangulation ]

        return self(triangulation)

    @cached_method
    def distance_affine(self, x, y):
        r"""
        Return the distance between two points.

        The distance function used in this method is `d_{aff}(x,y)^2`,
        the square of the usual affine distance function

        .. MATH::

            d_{aff}(x,y) = |x-y|

        INPUT:

        - ``x``, ``y`` -- two points of the point configuration

        OUTPUT:

        The metric distance-square `d_{aff}(x,y)^2`. Note that this
        distance lies in the same field as the entries of ``x``,
        ``y``. That is, the distance of rational points will be
        rational and so on.

        EXAMPLES::

            sage: pc = PointConfiguration([(0,0),(1,0),(2,1),(1,2),(0,1)])
            sage: [pc.distance_affine(pc.point(0), p) for p in pc.points()]
            [0, 1, 5, 5, 1]
        """
        self._assert_is_affine()
        d = 0
        for xi, yi in zip(x.projective(), y.projective()):
            d += (xi-yi)**2
        return d

    @cached_method
    def distance_FS(self, x, y):
        r"""
        Return the distance between two points.

        The distance function used in this method is `1-\cos
        d_{FS}(x,y)^2`, where `d_{FS}` is the Fubini-Study distance of
        projective points. Recall the Fubini-Studi distance function

        .. MATH::

            d_{FS}(x,y) = \arccos \sqrt{ \frac{(x\cdot y)^2}{|x|^2 |y|^2} }

        INPUT:

        - ``x``, ``y`` -- two points of the point configuration

        OUTPUT:

        The distance `1-\cos d_{FS}(x,y)^2`. Note that this distance
        lies in the same field as the entries of ``x``, ``y``. That
        is, the distance of rational points will be rational and so
        on.

        EXAMPLES::

            sage: pc = PointConfiguration([(0,0), (1,0), (2,1), (1,2), (0,1)])
            sage: [pc.distance_FS(pc.point(0), p) for p in pc.points()]
            [0, 1/2, 5/6, 5/6, 1/2]
        """
        x2 = y2 = xy = 0
        for xi, yi in zip(x.projective(), y.projective()):
            x2 += xi*xi
            y2 += yi*yi
            xy += xi*yi
        return 1-xy*xy/(x2*y2)

    @cached_method
    def distance(self, x, y):
        """
        Return the distance between two points.

        INPUT:

        - ``x``, ``y`` -- two points of the point configuration

        OUTPUT:

        The distance between ``x`` and ``y``, measured either with
        :meth:`distance_affine` or :meth:`distance_FS` depending on
        whether the point configuration is defined by affine or
        projective points. These are related, but not equal to the
        usual flat and Fubini-Study distance.

        EXAMPLES::

            sage: pc = PointConfiguration([(0,0), (1,0), (2,1), (1,2), (0,1)])
            sage: [pc.distance(pc.point(0), p) for p in pc.points()]
            [0, 1, 5, 5, 1]

            sage: pc = PointConfiguration([(0,0,1), (1,0,1), (2,1,1), (1,2,1), (0,1,1)],
            ....:                         projective=True)
            sage: [pc.distance(pc.point(0), p) for p in pc.points()]
            [0, 1/2, 5/6, 5/6, 1/2]
        """
        if self.is_affine():
            return self.distance_affine(x,y)
        else:
            return self.distance_FS(x,y)

    def farthest_point(self, points, among=None):
        """
        Return the point with the most distance from ``points``.

        INPUT:

        - ``points`` -- list of points

        - ``among`` -- list of points or ``None`` (default); the set
          of points from which to pick the farthest one. By default,
          all points of the configuration are considered.

        OUTPUT:

        A :class:`~sage.geometry.triangulation.base.Point` with
        largest minimal distance from all given ``points``.

        EXAMPLES::

            sage: pc = PointConfiguration([(0,0), (1,0), (1,1), (0,1)])
            sage: pc.farthest_point([pc.point(0)])
            P(1, 1)
        """
        if len(points) == 0:
            return self.point(0)
        if among is None:
            among = self.points()
        p_max = None
        for p in among:
            if p in points:
                continue
            if p_max is None:
                p_max = p
                d_max = min(self.distance(p,q) for q in points)
                continue
            d = min(self.distance(p,q) for q in points)
            if d > d_max:
                p_max = p
        return p_max

    def contained_simplex(self, large=True, initial_point=None, point_order=None):
        """
        Return a simplex contained in the point configuration.

        INPUT:

        - ``large`` -- boolean; whether to attempt to return a large
          simplex

        - ``initial_point`` -- a
          :class:`~sage.geometry.triangulation.base.Point` or ``None``
          (default). A specific point to start with when picking the
          simplex vertices.

        - ``point_order`` -- list or tuple of (some or all)
          :class:`~sage.geometry.triangulation.base.Point` s or ``None``
          (default)

        OUTPUT:

        A tuple of points that span a simplex of dimension
        :meth:`dim`. If ``large==True``, the simplex is constructed by
        successively picking the farthest point. This will ensure that
        the simplex is not unnecessarily small, but will in general
        not return a maximal simplex.
        If a ``point_order`` is specified, the simplex is greedily
        constructed by considering the points in this order.
        The ``large`` option and ``initial_point`` is ignored in this case.
        The ``point_order`` may contain only a subset of the points;
        in this case, the dimension of the simplex will be the dimension of
        this subset.

        EXAMPLES::

            sage: pc = PointConfiguration([(0,0), (1,0), (2,1), (1,1), (0,1)])
            sage: pc.contained_simplex()
            (P(0, 1), P(2, 1), P(1, 0))
            sage: pc.contained_simplex(large=False)
            (P(0, 1), P(1, 1), P(1, 0))
            sage: pc.contained_simplex(initial_point=pc.point(2))
            (P(2, 1), P(0, 0), P(1, 0))

            sage: pc = PointConfiguration([[0,0], [0,1], [1,0], [1,1], [-1,-1]])
            sage: pc.contained_simplex()
            (P(-1, -1), P(1, 1), P(0, 1))
            sage: pc.contained_simplex(point_order=[pc[1], pc[3], pc[4], pc[2], pc[0]])
            (P(0, 1), P(1, 1), P(-1, -1))

        Lower-dimensional example::

            sage: pc.contained_simplex(point_order=[pc[0], pc[3], pc[4]])
            (P(0, 0), P(1, 1))

        TESTS::

            sage: pc = PointConfiguration([[0,0], [0,1], [1,0]])
            sage: pc.contained_simplex()
            (P(1, 0), P(0, 1), P(0, 0))
            sage: pc = PointConfiguration([[0,0], [0,1]])
            sage: pc.contained_simplex()
            (P(0, 1), P(0, 0))
            sage: pc = PointConfiguration([[0,0]])
            sage: pc.contained_simplex()
            (P(0, 0),)
            sage: pc = PointConfiguration([])
            sage: pc.contained_simplex()
            ()
        """
        self._assert_is_affine()
        if point_order is None:
            points = list(self.points())
        else:
            points = list(reversed(point_order))
            # points are removed one by one from the end.
            initial_point = None
            large = False
            # If point_order is specified, the points of the
            # PointConfiguration are actually ignored.
        if not points:
            return tuple()

        if initial_point is None:
            origin = points.pop()
        else:
            origin = initial_point
            points.remove(origin)
        vertices = [origin]
        edges = []
        while points and len(vertices) <= self.dim():
            if large:
                p = self.farthest_point(vertices, points)
                points.remove(p)
            else:
                p = points.pop()
            edge = p.reduced_affine_vector() - origin.reduced_affine_vector()
            if edges and (ker * edge).is_zero():
                continue
            vertices.append(p)
            edges.append(edge)
            ker = matrix(edges).right_kernel().matrix()
        return tuple(vertices)

    def placing_triangulation(self, point_order=None):
        r"""
        Construct the placing (pushing) triangulation.

        INPUT:

        - ``point_order`` -- list of points or integers. The order in
          which the points are to be placed. If not given, the points
          will be placed in some arbitrary order that attempts to
          produce a small number of simplices.

        OUTPUT: a :class:`~sage.geometry.triangulation.triangulation.Triangulation`

        EXAMPLES::

            sage: pc = PointConfiguration([(0,0), (1,0), (2,1), (1,2), (0,1)])
            sage: pc.placing_triangulation()
            (<0,1,2>, <0,2,4>, <2,3,4>)
            sage: pc.placing_triangulation(point_order=(3,2,1,4,0))
            (<0,1,4>, <1,2,3>, <1,3,4>)
            sage: pc.placing_triangulation(point_order=[pc[1], pc[3], pc[4], pc[0]])
            (<0,1,4>, <1,3,4>)
            sage: U = matrix([
            ....:    [ 0, 0, 0, 0, 0, 2, 4,-1, 1, 1, 0, 0, 1, 0],
            ....:    [ 0, 0, 0, 1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0],
            ....:    [ 0, 2, 0, 0, 0, 0,-1, 0, 1, 0, 1, 0, 0, 1],
            ....:    [ 0, 1, 1, 0, 0, 1, 0,-2, 1, 0, 0,-1, 1, 1],
            ....:    [ 0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0]
            ....: ])
            sage: p = PointConfiguration(U.columns())
            sage: triangulation = p.placing_triangulation();  triangulation
            (<0,2,3,4,6,7>, <0,2,3,4,6,12>, <0,2,3,4,7,13>, <0,2,3,4,12,13>,
             <0,2,3,6,7,13>, <0,2,3,6,12,13>, <0,2,4,6,7,13>, <0,2,4,6,12,13>,
             <0,3,4,6,7,12>, <0,3,4,7,12,13>, <0,3,6,7,12,13>, <0,4,6,7,12,13>,
             <1,3,4,5,6,12>, <1,3,4,6,11,12>, <1,3,4,7,11,13>, <1,3,4,11,12,13>,
             <1,3,6,7,11,13>, <1,3,6,11,12,13>, <1,4,6,7,11,13>, <1,4,6,11,12,13>,
             <3,4,6,7,11,12>, <3,4,7,11,12,13>, <3,6,7,11,12,13>, <4,6,7,11,12,13>)
            sage: sum(p.volume(t) for t in triangulation)
            42
            sage: p0 = PointConfiguration([(0,0), (+1,0), (-1,0), (0,+1), (0,-1)])
            sage: p0.pushing_triangulation(point_order=[1,2,0,3,4])
            (<1,2,3>, <1,2,4>)
            sage: p0.pushing_triangulation(point_order=[0,1,2,3,4])
            (<0,1,3>, <0,1,4>, <0,2,3>, <0,2,4>)

        The same triangulation with renumbered points 0->4, 1->0, etc::

            sage: p1 = PointConfiguration([(+1,0), (-1,0), (0,+1), (0,-1), (0,0)])
            sage: p1.pushing_triangulation(point_order=[4,0,1,2,3])
            (<0,2,4>, <0,3,4>, <1,2,4>, <1,3,4>)
        """
        facet_normals = dict()

        def facets_of_simplex(simplex):
            """
            Return the facets of the simplex and store the normals in facet_normals
            """
            simplex = list(simplex)
            origin = simplex[0]
            rest = simplex[1:]
            span = matrix([origin.reduced_affine_vector()-p.reduced_affine_vector()
                           for p in rest])
            # span.inverse() linearly transforms the simplex into the unit simplex
            normals = span.inverse().columns()
            facets = []
            # The facets incident to the chosen vertex "origin"
            for opposing_vertex, normal in zip(rest, normals):
                facet = frozenset([origin] + [p for p in rest if p is not opposing_vertex])
                facets.append(facet)
                normal.set_immutable()
                facet_normals[facet] = normal
            # The remaining facet that is not incident to "origin"
            facet = frozenset(rest)
            normal = -sum(normals)
            normal.set_immutable()
            facet_normals[facet] = normal
            facets.append(facet)
            return set(facets)

        # input verification
        self._assert_is_affine()

        point_order_is_given = point_order is not None
        if point_order is None:
            point_order = list(self.points())
        elif isinstance(point_order[0], Point):
            point_order = list(point_order)
            assert all(p.point_configuration() is self for p in point_order),\
                "Some point in 'point_order' does not belong to the PointConfiguration."
        else:
            point_order = [self.point(i) for i in point_order]

        # construct the initial simplex
        if point_order_is_given:
            simplices = [frozenset(self.contained_simplex(large=False, point_order=point_order))]
        else:
            simplices = [frozenset(self.contained_simplex(large=True))]
        for s in simplices[0]:
            try:
                point_order.remove(s)
            except ValueError:
                pass
        facets = facets_of_simplex(simplices[0])

        # successively place the remaining points

        # TODO: In concordance with the heuristic to choose a LARGE starting simplex,
        # one could continue to try to pick points that are far away from the previous ones,
        # unless point_order_is_given.
        for point in point_order:
            # identify visible facets
            visible_facets = []
            for facet in facets:
                origin = next(iter(facet))
                normal = facet_normals[facet]
                v = point.reduced_affine_vector() - origin.reduced_affine_vector()
                if v * normal > 0:
                    visible_facets.append(facet)

            # construct simplices over each visible facet
            new_facets = set()
            for facet in visible_facets:
                simplex = frozenset(list(facet) + [point])
                simplices.append(simplex)
                for facet in facets_of_simplex(simplex):
                    if facet in visible_facets:
                        continue
                    if facet in new_facets:
                        new_facets.remove(facet)
                        continue
                    new_facets.add(facet)
            facets.difference_update(visible_facets)
            facets.update(new_facets)

        # construct the triangulation
        triangulation = [[p.index() for p in simplx] for simplx in simplices]
        return self(triangulation)

    pushing_triangulation = placing_triangulation

    @cached_method
    def Gale_transform(self, points=None):
        r"""
        Return the Gale transform of ``self``.

        INPUT:

        - ``points`` -- tuple of points or point indices or ``None``
          (default). A subset of points for which to compute the Gale
          transform. By default, all points are used.

        OUTPUT: a matrix over :meth:`base_ring`

        EXAMPLES::

            sage: pc = PointConfiguration([(0,0), (1,0), (2,1), (1,1), (0,1)])
            sage: pc.Gale_transform()
            [ 1 -1  0  1 -1]
            [ 0  0  1 -2  1]

            sage: pc.Gale_transform((0,1,3,4))
            [ 1 -1  1 -1]

            sage: points = (pc.point(0), pc.point(1), pc.point(3), pc.point(4))
            sage: pc.Gale_transform(points)
            [ 1 -1  1 -1]
        """
        self._assert_is_affine()
        if points is None:
            points = self.points()
        else:
            try:
                points = [ self.point(ZZ(i)) for i in points ]
            except TypeError:
                pass
        m = matrix([ (1,) + p.affine() for p in points])
        return m.left_kernel().matrix()

    def plot(self, **kwds):
        r"""
        Produce a graphical representation of the point configuration.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0], [0,1], [1,0], [1,1], [-1,-1]])
            sage: p.plot(axes=False)                                                    # needs sage.plot
            Graphics object consisting of 5 graphics primitives

        .. PLOT::
            :width: 300 px

            p = PointConfiguration([[0,0], [0,1], [1,0], [1,1], [-1,-1]])
            sphinx_plot(p.plot(axes=False))
        """
        return self.element_class([], parent=self, check=False).plot(**kwds)
