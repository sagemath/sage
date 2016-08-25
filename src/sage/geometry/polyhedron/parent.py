r"""
Parents for Polyhedra
"""
from __future__ import absolute_import

#*****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.parent import Parent
from sage.structure.element import get_coercion_model
from sage.structure.unique_representation import UniqueRepresentation
from sage.modules.free_module import is_FreeModule
from sage.misc.cachefunc import cached_method
from sage.rings.all import ZZ, QQ, RDF, CommutativeRing
from sage.categories.fields import Fields

from sage.geometry.polyhedron.base import Polyhedron_base, is_Polyhedron
from .representation import Inequality, Equation, Vertex, Ray, Line


def Polyhedra(base_ring, ambient_dim, backend=None):
    """
    Construct a suitable parent class for polyhedra

    INPUT:

    - ``base_ring`` -- A ring. Currently there are backends for `\ZZ`,
      `\QQ`, and `\RDF`.

    - ``ambient_dim`` -- integer. The ambient space dimension.

    - ``backend`` -- string. The name of the backend for computations. Currently there are two backends implemented:

         * ``backend=ppl`` uses the Parma Polyhedra Library

         * ``backend=cdd`` uses CDD

    OUTPUT:

    A parent class for polyhedra over the given base ring if the
    backend supports it. If not, the parent base ring can be larger
    (for example, `\QQ` instead of `\ZZ`). If there is no
    implementation at all, a ``ValueError`` is raised.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.parent import Polyhedra
        sage: Polyhedra(AA, 3)
        Polyhedra in AA^3
        sage: Polyhedra(ZZ, 3)
        Polyhedra in ZZ^3
        sage: type(_)
        <class 'sage.geometry.polyhedron.parent.Polyhedra_ZZ_ppl_with_category'>
        sage: Polyhedra(QQ, 3, backend='cdd')
        Polyhedra in QQ^3
        sage: type(_)
        <class 'sage.geometry.polyhedron.parent.Polyhedra_QQ_cdd_with_category'>

    CDD does not support integer polytopes directly::

        sage: Polyhedra(ZZ, 3, backend='cdd')
        Polyhedra in QQ^3
    """
    if backend is None:
        if base_ring is ZZ:
            backend = 'ppl'
        elif base_ring is QQ:
            backend = 'ppl'
        elif base_ring is RDF:
            backend = 'cdd'
        else:
            backend = 'field'
    if backend == 'ppl' and base_ring is QQ:
        return Polyhedra_QQ_ppl(base_ring, ambient_dim)
    elif backend == 'ppl' and base_ring is ZZ:
        return Polyhedra_ZZ_ppl(base_ring, ambient_dim)
    elif backend == 'cdd' and base_ring in (ZZ, QQ):
        return Polyhedra_QQ_cdd(QQ, ambient_dim)
    elif backend == 'cdd' and base_ring is RDF:
        return Polyhedra_RDF_cdd(RDF, ambient_dim)
    elif backend == 'field':
        return Polyhedra_field(base_ring.fraction_field(), ambient_dim)
    else:
        raise ValueError('No such backend (='+str(backend)+
                         ') implemented for given basering (='+str(base_ring)+').')



class Polyhedra_base(UniqueRepresentation, Parent):
    r"""
    Polyhedra in a fixed ambient space.

    INPUT:

    - ``base_ring`` -- either ``ZZ``, ``QQ``, or ``RDF``. The base
      ring of the ambient module/vector space.

    - ``ambient_dim`` -- integer. The ambient space dimension.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.parent import Polyhedra
        sage: Polyhedra(ZZ, 3)
        Polyhedra in ZZ^3
    """
    def __init__(self, base_ring, ambient_dim):
        """
        The Python constructor.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: Polyhedra(QQ, 3)
            Polyhedra in QQ^3

        TESTS::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: P = Polyhedra(QQ, 3)
            sage: TestSuite(P).run(skip='_test_pickling')
        """
        self._ambient_dim = ambient_dim
        from sage.categories.polyhedra import PolyhedralSets
        Parent.__init__(self, base=base_ring, category=PolyhedralSets(base_ring))
        self._Inequality_pool = []
        self._Equation_pool = []
        self._Vertex_pool = []
        self._Ray_pool = []
        self._Line_pool = []

    def recycle(self, polyhedron):
        """
        Recycle the H/V-representation objects of a polyhedron.

        This speeds up creation of new polyhedra by reusing
        objects. After recycling a polyhedron object, it is not in a
        consistent state any more and neither the polyhedron nor its
        H/V-representation objects may be used any more.

        INPUT:

        - ``polyhedron`` -- a polyhedron whose parent is ``self``.


        EXAMPLES::

            sage: p = Polyhedron([(0,0),(1,0),(0,1)])
            sage: p.parent().recycle(p)

        TESTS::

            sage: p = Polyhedron([(0,0),(1,0),(0,1)])
            sage: n = len(p.parent()._Vertex_pool)
            sage: p._delete()
            sage: len(p.parent()._Vertex_pool) - n
            3
        """
        if self is not polyhedron.parent():
            raise TypeError('The polyhedron has the wrong parent class.')
        self._Inequality_pool.extend(polyhedron.inequalities())
        self._Equation_pool.extend(polyhedron.equations())
        self._Vertex_pool.extend(polyhedron.vertices())
        self._Ray_pool.extend(polyhedron.rays())
        self._Line_pool.extend(polyhedron.lines())
        for Hrep in polyhedron.Hrep_generator():
            Hrep._polyhedron = None
        for Vrep in polyhedron.Vrep_generator():
            Vrep._polyhedron = None
        polyhedron._Hrepresentation = None
        polyhedron._Vrepresentation = None

    def ambient_dim(self):
        r"""
        Return the dimension of the ambient space.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: Polyhedra(QQ, 3).ambient_dim()
            3
        """
        return self._ambient_dim

    @cached_method
    def an_element(self):
        r"""
        Returns a Polyhedron.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: Polyhedra(QQ, 4).an_element()
            A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 5 vertices
        """
        zero = self.base_ring().zero()
        one  = self.base_ring().one()
        p = [zero] * self.ambient_dim()
        points = [p]
        for i in range(0,self.ambient_dim()):
            p = [zero] * self.ambient_dim()
            p[i] = one
            points.append(p)
        return self.element_class(self, [points,[],[]], None)

    @cached_method
    def some_elements(self):
        r"""
        Returns a list of some elements of the semigroup.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: Polyhedra(QQ, 4).some_elements()
            [A 3-dimensional polyhedron in QQ^4 defined as the convex hull of 4 vertices,
             A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 1 vertex and 4 rays,
             A 2-dimensional polyhedron in QQ^4 defined as the convex hull of 2 vertices and 1 ray,
             The empty polyhedron in QQ^4]
            sage: Polyhedra(ZZ,0).some_elements()
            [The empty polyhedron in ZZ^0,
             A 0-dimensional polyhedron in ZZ^0 defined as the convex hull of 1 vertex]
        """
        if self.ambient_dim() == 0:
            return [
                self.element_class(self, None, None),
                self.element_class(self, None, [[],[]]) ]
        points = []
        R = self.base_ring()
        for i in range(0,self.ambient_dim()+5):
            points.append([R(i*j^2) for j in range(0,self.ambient_dim())])
        return [
            self.element_class(self, [points[0:self.ambient_dim()+1], [], []], None),
            self.element_class(self, [points[0:1], points[1:self.ambient_dim()+1], []], None),
            self.element_class(self, [points[0:3], points[4:5], []], None),
            self.element_class(self, None, None) ]

    @cached_method
    def zero(self):
        r"""
        Return the polyhedron consisting of the origin, which is the
        neutral element for Minkowski addition.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: p = Polyhedra(QQ, 4).zero();  p
            A 0-dimensional polyhedron in QQ^4 defined as the convex hull of 1 vertex
            sage: p+p == p
            True
        """
        Vrep = [[[self.base_ring().zero()]*self.ambient_dim()], [], []]
        return self.element_class(self, Vrep, None)

    def empty(self):
        """
        Return the empty polyhedron.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: P = Polyhedra(QQ, 4)
            sage: P.empty()
            The empty polyhedron in QQ^4
            sage: P.empty().is_empty()
            True
        """
        return self(None, None)

    def universe(self):
        """
        Return the entire ambient space as polyhedron.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: P = Polyhedra(QQ, 4)
            sage: P.universe()
            A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 1 vertex and 4 lines
            sage: P.universe().is_universe()
            True
        """
        R = self.base_ring()
        return self(None, [[[R.one()]+[R.zero()]*self.ambient_dim()], []], convert=True)

    @cached_method
    def Vrepresentation_space(self):
        r"""
        Return the ambient vector space.

        This is the vector space or module containing the
        Vrepresentation vectors.

        OUTPUT:

        A free module over the base ring of dimension :meth:`ambient_dim`.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: Polyhedra(QQ, 4).Vrepresentation_space()
            Vector space of dimension 4 over Rational Field
            sage: Polyhedra(QQ, 4).ambient_space()
            Vector space of dimension 4 over Rational Field
        """
        if self.base_ring() in Fields():
            from sage.modules.free_module import VectorSpace
            return VectorSpace(self.base_ring(), self.ambient_dim())
        else:
            from sage.modules.free_module import FreeModule
            return FreeModule(self.base_ring(), self.ambient_dim())

    ambient_space = Vrepresentation_space

    @cached_method
    def Hrepresentation_space(self):
        r"""
        Return the linear space containing the H-representation vectors.

        OUTPUT:

        A free module over the base ring of dimension :meth:`ambient_dim` + 1.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: Polyhedra(ZZ, 2).Hrepresentation_space()
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
        """
        if self.base_ring() in Fields():
            from sage.modules.free_module import VectorSpace
            return VectorSpace(self.base_ring(), self.ambient_dim()+1)
        else:
            from sage.modules.free_module import FreeModule
            return FreeModule(self.base_ring(), self.ambient_dim()+1)

    def _repr_ambient_module(self):
        """
        Return an abbreviated string representation of the ambient
        space.

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: Polyhedra(QQ, 3)._repr_ambient_module()
            'QQ^3'
            sage: K.<sqrt3> = NumberField(x^2-3)
            sage: Polyhedra(K, 4)._repr_ambient_module()
            '(Number Field in sqrt3 with defining polynomial x^2 - 3)^4'
        """
        from sage.rings.qqbar import AA
        if self.base_ring() is ZZ:
            s = 'ZZ'
        elif self.base_ring() is QQ:
            s = 'QQ'
        elif self.base_ring() is RDF:
            s = 'RDF'
        elif self.base_ring() is AA:
            s = 'AA'
        else:
            s = '({0})'.format(self.base_ring())
        s += '^' + repr(self.ambient_dim())
        return s

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: Polyhedra(QQ, 3)
            Polyhedra in QQ^3
            sage: Polyhedra(QQ, 3)._repr_()
            'Polyhedra in QQ^3'
        """
        return 'Polyhedra in '+self._repr_ambient_module()

    def _element_constructor_(self, *args, **kwds):
        """
        The element (polyhedron) constructor.

        INPUT:

        - ``Vrep`` -- a list `[vertices, rays, lines]`` or ``None``.

        - ``Hrep`` -- a list `[ieqs, eqns]`` or ``None``.

        - ``convert`` -- boolean keyword argument (default:
          ``True``). Whether to convert the cooordinates into the base
          ring.

        - ``**kwds`` -- optional remaining keywords that are passed to the
          polyhedron constructor.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: P = Polyhedra(QQ, 3)
            sage: P._element_constructor_([[(0,0,0),(1,0,0),(0,1,0),(0,0,1)], [], []], None)
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 4 vertices
            sage: P([[(0,0,0),(1,0,0),(0,1,0),(0,0,1)], [], []], None)
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 4 vertices
            sage: P(0)
            A 0-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex
        """
        nargs = len(args)
        convert = kwds.pop('convert', True)
        if nargs==2:
            Vrep, Hrep = args
            def convert_base_ring(lstlst):
                return [ [self.base_ring()(x) for x in lst] for lst in lstlst]
            if convert and Hrep:
                Hrep = [convert_base_ring(_) for _ in Hrep]
            if convert and Vrep:
                Vrep = [convert_base_ring(_) for _ in Vrep]
            return self.element_class(self, Vrep, Hrep, **kwds)
        if nargs==1 and is_Polyhedron(args[0]):
            polyhedron = args[0]
            Hrep = [ polyhedron.inequality_generator(), polyhedron.equation_generator() ]
            return self.element_class(self, None, Hrep, **kwds)
        if nargs==1 and args[0]==0:
            return self.zero()
        raise ValueError('Cannot convert to polyhedron object.')

    def base_extend(self, base_ring, backend=None):
        """
        Return the base extended parent.

        INPUT:

        - ``base_ring``, ``backend`` -- see
          :func:`~sage.geometry.polyhedron.constructor.Polyhedron`.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: Polyhedra(ZZ,3).base_extend(QQ)
            Polyhedra in QQ^3
            sage: Polyhedra(ZZ,3).an_element().base_extend(QQ)
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 4 vertices
        """
        if self.base_ring().has_coerce_map_from(base_ring):
            return self
        elif base_ring.has_coerce_map_from(self.base_ring()):
            return Polyhedra(base_ring, self.ambient_dim())

    def _coerce_base_ring(self, other):
        """
        Return the common base rincg for both ``self`` and ``other``.

        This method is not part of the coercion framework, but only a
        convenience function for :class:`Polyhedra_base`.

        INPUT:

        - ``other`` -- must be either:

            * another ``Polyhedron`` object

            * `\ZZ`, `\QQ`, `RDF`, or a ring that can be coerced into them.

            * a constant that can be coerced to `\ZZ`, `\QQ`, or `RDF`.

        OUTPUT:

        Either `\ZZ`, `\QQ`, or `RDF`. Raises ``TypeError`` if
        ``other`` is not a suitable input.

        .. NOTE::

            "Real" numbers in sage are not necessarily elements of
            `RDF`. For example, the literal `1.0` is not.

        EXAMPLES::

            sage: triangle_QQ  = Polyhedron(vertices = [[1,0],[0,1],[1,1]], base_ring=QQ).parent()
            sage: triangle_RDF = Polyhedron(vertices = [[1,0],[0,1],[1,1]], base_ring=RDF).parent()
            sage: triangle_QQ._coerce_base_ring(QQ)
            Rational Field
            sage: triangle_QQ._coerce_base_ring(triangle_RDF)
            Real Double Field
            sage: triangle_RDF._coerce_base_ring(triangle_QQ)
            Real Double Field
            sage: triangle_QQ._coerce_base_ring(RDF)
            Real Double Field
            sage: triangle_QQ._coerce_base_ring(ZZ)
            Rational Field
            sage: triangle_QQ._coerce_base_ring(1/2)
            Rational Field
            sage: triangle_QQ._coerce_base_ring(0.5)
            Real Double Field
        """
        try:
            other_ring = other.base_ring()
        except AttributeError:
            try:
                # other is a constant?
                other_ring = other.parent()
            except AttributeError:
                other_ring = None
                for ring in (ZZ, QQ, RDF):
                    try:
                        ring.coerce(other)
                        other_ring = ring
                        break
                    except TypeError:
                        pass
                if other_ring is None:
                    raise TypeError('Could not coerce '+str(other)+' into ZZ, QQ, or RDF.')

        if not other_ring.is_exact():
            other_ring = RDF  # the only supported floating-point numbers for now

        cm_map, cm_ring = get_coercion_model().analyse(self.base_ring(), other_ring)
        if cm_ring is None:
            raise TypeError('Could not coerce type '+str(other)+' into ZZ, QQ, or RDF.')
        return cm_ring

    def _coerce_map_from_(self, X):
        r"""
        Return whether there is a coercion from ``X``

        INPUT:

        - ``X`` -- anything.

        OUTPUT:

        Boolean.

        EXAMPLE::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: Polyhedra(QQ,3).has_coerce_map_from( Polyhedra(ZZ,3) )   # indirect doctest
            True
            sage: Polyhedra(ZZ,3).has_coerce_map_from( Polyhedra(QQ,3) )
            False
        """
        if not isinstance(X, Polyhedra_base):
            return False
        if self.ambient_dim() != X.ambient_dim():
            return False
        return self.base_ring().has_coerce_map_from(X.base_ring())

    def _get_action_(self, other, op, self_is_left):
        """
        Register actions with the coercion model.

        The monoid actions are Minkowski sum and Cartesian product. In
        addition, we want multiplication by a scalar to be dilation
        and addition by a vector to be translation. This is
        implemented as an action in the coercion model.

        INPUT:

        - ``other`` -- a scalar or a vector.

        - ``op`` -- the operator.

        - ``self_is_left`` -- boolean. Whether ``self`` is on the left
          of the operator.

        OUTPUT:

        An action that is used by the coercion model.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: PZZ2 = Polyhedra(ZZ, 2)
            sage: PZZ2.get_action(ZZ)   # indirect doctest
            Right action by Integer Ring on Polyhedra in ZZ^2
            sage: PZZ2.get_action(QQ)
            Right action by Rational Field on Polyhedra in QQ^2
            with precomposition on left by Conversion map:
              From: Polyhedra in ZZ^2
              To:   Polyhedra in QQ^2
            with precomposition on right by Identity endomorphism of Rational Field
            sage: PQQ2 = Polyhedra(QQ, 2)
            sage: PQQ2.get_action(ZZ)
            Right action by Integer Ring on Polyhedra in QQ^2
            sage: PQQ2.get_action(QQ)
            Right action by Rational Field on Polyhedra in QQ^2

            sage: Polyhedra(ZZ,2).an_element() * 2
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices
            sage: Polyhedra(ZZ,2).an_element() * (2/3)
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices
            sage: Polyhedra(QQ,2).an_element() * 2
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices
            sage: Polyhedra(QQ,2).an_element() * (2/3)
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices

            sage: 2     * Polyhedra(ZZ,2).an_element()
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices
            sage: (2/3) * Polyhedra(ZZ,2).an_element()
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices
            sage: 2     * Polyhedra(QQ,2).an_element()
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices
            sage: (2/3) * Polyhedra(QQ,2).an_element()
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: PZZ2.get_action(ZZ^2, op=operator.add)
            Right action by Ambient free module of rank 2 over the principal ideal domain Integer Ring on Polyhedra in ZZ^2
            with precomposition on left by Identity endomorphism of Polyhedra in ZZ^2
            with precomposition on right by Generic endomorphism of Ambient free module of rank 2 over the principal ideal domain Integer Ring

        """
        import operator
        from sage.structure.coerce_actions import ActedUponAction
        from sage.categories.action import PrecomposedAction

        if op is operator.add and is_FreeModule(other):
            base_ring = self._coerce_base_ring(other)
            extended_self = self.base_extend(base_ring)
            extended_other = other.base_extend(base_ring)
            action = ActedUponAction(extended_other, extended_self, not self_is_left)
            if self_is_left:
                action = PrecomposedAction(action,
                                           extended_self._internal_coerce_map_from(self).__copy__(),
                                           extended_other._internal_coerce_map_from(other).__copy__())
            else:
                action = PrecomposedAction(action,
                                           extended_other._internal_coerce_map_from(other).__copy__(),
                                           extended_self._internal_coerce_map_from(self).__copy__())
            return action

        if op is operator.mul and isinstance(other, CommutativeRing):
            ring = self._coerce_base_ring(other)
            if ring is self.base_ring():
                return ActedUponAction(other, self, not self_is_left)
            extended = self.base_extend(ring)
            action = ActedUponAction(ring, extended, not self_is_left)
            if self_is_left:
                action = PrecomposedAction(action,
                                           extended._internal_coerce_map_from(self).__copy__(),
                                           ring._internal_coerce_map_from(other).__copy__())
            else:
                action = PrecomposedAction(action,
                                           ring._internal_coerce_map_from(other).__copy__(),
                                           extended._internal_coerce_map_from(self).__copy__())
            return action

    def _make_Inequality(self, polyhedron, data):
        """
        Create a new inequality object.

        INPUT:

        - ``polyhedron`` -- the new polyhedron.

        - ``data`` -- the H-representation data.

        OUTPUT:

        A new :class:`~sage.geometry.polyhedron.representation.Inequality` object.

        EXAMPLES::

            sage: p = Polyhedron([(1,2,3),(2/3,3/4,4/5)])   # indirect doctest
            sage: next(p.inequality_generator())
            An inequality (0, 0, -1) x + 3 >= 0
        """
        try:
            obj = self._Inequality_pool.pop()
        except IndexError:
            obj = Inequality(self)
        obj._set_data(polyhedron, data)
        return obj

    def _make_Equation(self, polyhedron, data):
        """
        Create a new equation object.

        INPUT:

        - ``polyhedron`` -- the new polyhedron.

        - ``data`` -- the H-representation data.

        OUTPUT:

        A new :class:`~sage.geometry.polyhedron.representation.Equation` object.

        EXAMPLES::

            sage: p = Polyhedron([(1,2,3),(2/3,3/4,4/5)])   # indirect doctest
            sage: next(p.equation_generator())
            An equation (0, 44, -25) x - 13 == 0
        """
        try:
            obj = self._Equation_pool.pop()
        except IndexError:
            obj = Equation(self)
        obj._set_data(polyhedron, data)
        return obj

    def _make_Vertex(self, polyhedron, data):
        """
        Create a new vertex object.

        INPUT:

        - ``polyhedron`` -- the new polyhedron.

        - ``data`` -- the V-representation data.

        OUTPUT:

        A new :class:`~sage.geometry.polyhedron.representation.Vertex` object.

        EXAMPLES::

            sage: p = Polyhedron([(1,2,3),(2/3,3/4,4/5)], rays=[(5/6,6/7,7/8)])   # indirect doctest
            sage: next(p.vertex_generator())
            A vertex at (1, 2, 3)
        """
        try:
            obj = self._Vertex_pool.pop()
        except IndexError:
            obj = Vertex(self)
        obj._set_data(polyhedron, data)
        return obj

    def _make_Ray(self, polyhedron, data):
        """
        Create a new ray object.

        INPUT:

        - ``polyhedron`` -- the new polyhedron.

        - ``data`` -- the V-representation data.

        OUTPUT:

        A new :class:`~sage.geometry.polyhedron.representation.Ray` object.

        EXAMPLES::

            sage: p = Polyhedron([(1,2,3),(2/3,3/4,4/5)], rays=[(5/6,6/7,7/8)])   # indirect doctest
            sage: next(p.ray_generator())
            A ray in the direction (140, 144, 147)
        """
        try:
            obj = self._Ray_pool.pop()
        except IndexError:
            obj = Ray(self)
        obj._set_data(polyhedron, data)
        return obj

    def _make_Line(self, polyhedron, data):
        """
        Create a new line object.

        INPUT:

        - ``polyhedron`` -- the new polyhedron.

        - ``data`` -- the V-representation data.

        OUTPUT:

        A new :class:`~sage.geometry.polyhedron.representation.Line` object.

        EXAMPLES::

            sage: p = Polyhedron([(1,2,3),(2/3,3/4,4/5)], lines=[(5/6,6/7,7/8)])   # indirect doctest
            sage: next(p.line_generator())
            A line in the direction (140, 144, 147)
        """
        try:
            obj = self._Line_pool.pop()
        except IndexError:
            obj = Line(self)
        obj._set_data(polyhedron, data)
        return obj



from sage.geometry.polyhedron.backend_cdd import Polyhedron_QQ_cdd, Polyhedron_RDF_cdd
from sage.geometry.polyhedron.backend_ppl import Polyhedron_ZZ_ppl, Polyhedron_QQ_ppl
from sage.geometry.polyhedron.backend_field import Polyhedron_field

class Polyhedra_ZZ_ppl(Polyhedra_base):
    Element = Polyhedron_ZZ_ppl

class Polyhedra_QQ_ppl(Polyhedra_base):
    Element = Polyhedron_QQ_ppl

class Polyhedra_QQ_cdd(Polyhedra_base):
    Element = Polyhedron_QQ_cdd

class Polyhedra_RDF_cdd(Polyhedra_base):
    Element = Polyhedron_RDF_cdd

class Polyhedra_field(Polyhedra_base):
    Element = Polyhedron_field

