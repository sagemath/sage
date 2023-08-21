r"""
Parents for Polyhedra
"""

# ****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.structure.parent import Parent
from sage.structure.element import get_coercion_model
from sage.structure.unique_representation import UniqueRepresentation
from sage.modules.free_module import FreeModule, is_FreeModule
from sage.misc.cachefunc import cached_method, cached_function
from sage.misc.lazy_import import lazy_import
import sage.rings.abc
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_double import RDF
from sage.categories.fields import Fields
from sage.categories.rings import Rings
from sage.categories.modules import Modules

from sage.geometry.polyhedron.base import is_Polyhedron
from .representation import Inequality, Equation, Vertex, Ray, Line


def Polyhedra(ambient_space_or_base_ring=None, ambient_dim=None, backend=None, *,
              ambient_space=None, base_ring=None):
    r"""
    Construct a suitable parent class for polyhedra

    INPUT:

    - ``base_ring`` -- A ring. Currently there are backends for `\ZZ`,
      `\QQ`, and `\RDF`.

    - ``ambient_dim`` -- integer. The ambient space dimension.

    - ``ambient_space`` -- A free module.

    - ``backend`` -- string. The name of the backend for computations. There are
       several backends implemented:

         * ``backend="ppl"`` uses the Parma Polyhedra Library

         * ``backend="cdd"`` uses CDD

         * ``backend="normaliz"`` uses normaliz

         * ``backend="polymake"`` uses polymake

         * ``backend="field"`` a generic Sage implementation

    OUTPUT:

    A parent class for polyhedra over the given base ring if the
    backend supports it. If not, the parent base ring can be larger
    (for example, `\QQ` instead of `\ZZ`). If there is no
    implementation at all, a ``ValueError`` is raised.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.parent import Polyhedra
        sage: Polyhedra(AA, 3)                                                          # needs sage.rings.number_field
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

    Using a more general form of the constructor::

        sage: V = VectorSpace(QQ, 3)
        sage: Polyhedra(V) is Polyhedra(QQ, 3)
        True
        sage: Polyhedra(V, backend='field') is Polyhedra(QQ, 3, 'field')
        True
        sage: Polyhedra(backend='field', ambient_space=V) is Polyhedra(QQ, 3, 'field')
        True

        sage: M = FreeModule(ZZ, 2)
        sage: Polyhedra(M, backend='ppl') is Polyhedra(ZZ, 2, 'ppl')
        True

    TESTS::

        sage: Polyhedra(RR, 3, backend='field')                                         # needs sage.rings.real_mpfr
        Traceback (most recent call last):
        ...
        ValueError: the 'field' backend for polyhedron cannot be used with non-exact fields
        sage: Polyhedra(RR, 3)                                                          # needs sage.rings.real_mpfr
        Traceback (most recent call last):
        ...
        ValueError: no default backend for computations with Real Field with 53 bits of precision
        sage: Polyhedra(QQ[I], 2)                                                       # needs sage.rings.number_field
        Traceback (most recent call last):
        ...
        ValueError: invalid base ring: Number Field in I
        with defining polynomial x^2 + 1 with I = 1*I
        cannot be coerced to a real field
        sage: Polyhedra(AA, 3, backend='polymake')      # optional - jupymake           # needs sage.rings.number_field
        Traceback (most recent call last):
        ...
        ValueError: the 'polymake' backend for polyhedron cannot be used with Algebraic Real Field

        sage: Polyhedra(QQ, 2, backend='normaliz')
        Polyhedra in QQ^2
        sage: Polyhedra(SR, 2, backend='normaliz')      # optional - pynormaliz         # needs sage.symbolic
        Polyhedra in (Symbolic Ring)^2
        sage: SCR = SR.subring(no_variables=True)                                       # needs sage.symbolic
        sage: Polyhedra(SCR, 2, backend='normaliz')     # optional - pynormaliz         # needs sage.symbolic
        Polyhedra in (Symbolic Constants Subring)^2

        sage: Polyhedra(SCR, 2, backend='number_field')                                 # needs sage.symbolic
        Polyhedra in (Symbolic Constants Subring)^2

    """
    if ambient_space_or_base_ring is not None:
        if ambient_space_or_base_ring in Rings():
            base_ring = ambient_space_or_base_ring
        else:
            ambient_space = ambient_space_or_base_ring
    if ambient_space is not None:
        if ambient_space not in Modules:
            # There is no category of free modules, unfortunately
            # (see https://github.com/sagemath/sage/issues/30164)...
            raise ValueError('ambient_space must be a free module')
        if base_ring is None:
            base_ring = ambient_space.base_ring()
        if ambient_dim is None:
            try:
                ambient_dim = ambient_space.rank()
            except AttributeError:
                # ... so we test whether it is free using the existence of
                # a rank method
                raise ValueError('ambient_space must be a free module')
        if ambient_space is not FreeModule(base_ring, ambient_dim):
            raise NotImplementedError('ambient_space must be a standard free module')
    if backend is None:
        if base_ring is ZZ or base_ring is QQ:
            backend = 'ppl'
        elif base_ring is RDF:
            backend = 'cdd'
        elif base_ring.is_exact():
            # TODO: find a more robust way of checking that the coefficients are indeed
            # real numbers
            if not RDF.has_coerce_map_from(base_ring):
                raise ValueError("invalid base ring: {} cannot be coerced to a real field".format(base_ring))
            backend = 'field'
        else:
            raise ValueError("no default backend for computations with {}".format(base_ring))

    if backend == 'ppl' and base_ring is QQ:
        return Polyhedra_QQ_ppl(base_ring, ambient_dim, backend)
    elif backend == 'ppl' and base_ring is ZZ:
        return Polyhedra_ZZ_ppl(base_ring, ambient_dim, backend)
    elif backend == 'normaliz' and base_ring is QQ:
        return Polyhedra_QQ_normaliz(base_ring, ambient_dim, backend)
    elif backend == 'normaliz' and base_ring is ZZ:
        return Polyhedra_ZZ_normaliz(base_ring, ambient_dim, backend)
    elif backend == 'normaliz' and (isinstance(base_ring, sage.rings.abc.SymbolicRing) or base_ring.is_exact()):
        return Polyhedra_normaliz(base_ring, ambient_dim, backend)
    elif backend == 'cdd' and base_ring in (ZZ, QQ):
        return Polyhedra_QQ_cdd(QQ, ambient_dim, backend)
    elif backend == 'cdd' and base_ring is RDF:
        return Polyhedra_RDF_cdd(RDF, ambient_dim, backend)
    elif backend == 'polymake':
        base_field = base_ring.fraction_field()
        try:
            from sage.interfaces.polymake import polymake, PolymakeElement
            polymake_base_field = polymake(base_field)
            assert isinstance(polymake_base_field, PolymakeElement)  # to muffle pyflakes
        except TypeError:
            raise ValueError(f"the 'polymake' backend for polyhedron cannot be used with {base_field}")
        return Polyhedra_polymake(base_field, ambient_dim, backend)
    elif backend == 'number_field':
        return Polyhedra_number_field(base_ring.fraction_field(), ambient_dim, backend)
    elif backend == 'field':
        if not base_ring.is_exact():
            raise ValueError("the 'field' backend for polyhedron cannot be used with non-exact fields")
        return Polyhedra_field(base_ring.fraction_field(), ambient_dim, backend)
    else:
        raise ValueError('No such backend (=' + str(backend) +
                         ') implemented for given basering (=' + str(base_ring) + ').')


class Polyhedra_base(UniqueRepresentation, Parent):
    r"""
    Polyhedra in a fixed ambient space.

    INPUT:

    - ``base_ring`` -- either ``ZZ``, ``QQ``, or ``RDF``. The base
      ring of the ambient module/vector space.

    - ``ambient_dim`` -- integer. The ambient space dimension.

    - ``backend`` -- string. The name of the backend for computations. There are
       several backends implemented:

      * ``backend="ppl"`` uses the Parma Polyhedra Library

      * ``backend="cdd"`` uses CDD

      * ``backend="normaliz"`` uses normaliz

      * ``backend="polymake"`` uses polymake

      * ``backend="field"`` a generic Sage implementation

    EXAMPLES::

        sage: from sage.geometry.polyhedron.parent import Polyhedra
        sage: Polyhedra(ZZ, 3)
        Polyhedra in ZZ^3
    """
    def __init__(self, base_ring, ambient_dim, backend):
        """
        The Python constructor.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: Polyhedra(QQ, 3)
            Polyhedra in QQ^3

        TESTS::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: P = Polyhedra(QQ, 3)
            sage: TestSuite(P).run()
            sage: P = Polyhedra(QQ, 0)
            sage: TestSuite(P).run()
        """
        self._backend = backend
        self._ambient_dim = ambient_dim
        from sage.categories.polyhedra import PolyhedralSets
        from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
        category = PolyhedralSets(base_ring)
        if ambient_dim == 0:
            category = category & FiniteEnumeratedSets()
        else:
            category = category.Infinite()

        Parent.__init__(self, base=base_ring, category=category)
        self._Inequality_pool = []
        self._Equation_pool = []
        self._Vertex_pool = []
        self._Ray_pool = []
        self._Line_pool = []

    def list(self):
        """
        Return the two polyhedra in ambient dimension 0, raise an error otherwise

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: P = Polyhedra(QQ, 3)
            sage: P.cardinality()
            +Infinity

            sage: # needs sage.rings.number_field
            sage: P = Polyhedra(AA, 0)
            sage: P.category()
            Category of finite enumerated polyhedral sets over Algebraic Real Field
            sage: P.list()
            [The empty polyhedron in AA^0,
             A 0-dimensional polyhedron in AA^0 defined as the convex hull of 1 vertex]
            sage: P.cardinality()
            2
        """
        if self.ambient_dim():
            raise NotImplementedError
        return [self.empty(), self.universe()]

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
        if polyhedron.is_mutable():
            polyhedron._dependent_objects = []

    def ambient_dim(self):
        r"""
        Return the dimension of the ambient space.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: Polyhedra(QQ, 3).ambient_dim()
            3
        """
        return self._ambient_dim

    def backend(self):
        r"""
        Return the backend.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: Polyhedra(QQ, 3).backend()
            'ppl'
        """
        return self._backend

    @cached_method
    def an_element(self):
        r"""
        Return a Polyhedron.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: Polyhedra(QQ, 4).an_element()
            A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 5 vertices
        """
        zero = self.base_ring().zero()
        one = self.base_ring().one()
        p = [zero] * self.ambient_dim()
        points = [p]
        for i in range(self.ambient_dim()):
            p = [zero] * self.ambient_dim()
            p[i] = one
            points.append(p)
        return self.element_class(self, [points, [], []], None)

    @cached_method
    def some_elements(self):
        r"""
        Return a list of some elements of the semigroup.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: Polyhedra(QQ, 4).some_elements()
            [A 3-dimensional polyhedron in QQ^4
              defined as the convex hull of 4 vertices,
             A 4-dimensional polyhedron in QQ^4
              defined as the convex hull of 1 vertex and 4 rays,
             A 2-dimensional polyhedron in QQ^4
              defined as the convex hull of 2 vertices and 1 ray,
             The empty polyhedron in QQ^4]
            sage: Polyhedra(ZZ, 0).some_elements()
            [The empty polyhedron in ZZ^0,
             A 0-dimensional polyhedron in ZZ^0 defined as the convex hull of 1 vertex]
        """
        if self.ambient_dim() == 0:
            return [
                self.element_class(self, None, None),
                self.element_class(self, None, [[], []])]
        points = []
        R = self.base_ring()
        for i in range(self.ambient_dim() + 5):
            points.append([R(i*j ^ 2) for j in range(self.ambient_dim())])
        return [
            self.element_class(self, [points[0:self.ambient_dim() + 1], [], []], None),
            self.element_class(self, [points[0:1], points[1:self.ambient_dim() + 1], []], None),
            self.element_class(self, [points[0:3], points[4:5], []], None),
            self.element_class(self, None, None)]

    @cached_method
    def zero(self):
        r"""
        Return the polyhedron consisting of the origin, which is the
        neutral element for Minkowski addition.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: p = Polyhedra(QQ, 4).zero();  p
            A 0-dimensional polyhedron in QQ^4 defined as the convex hull of 1 vertex
            sage: p + p == p
            True
        """
        Vrep = [[[self.base_ring().zero()] * self.ambient_dim()], [], []]
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
            A 4-dimensional polyhedron in QQ^4 defined as
             the convex hull of 1 vertex and 4 lines
            sage: P.universe().is_universe()
            True
        """
        R = self.base_ring()
        return self(None, [[[R.one()] + [R.zero()] * self.ambient_dim()], []], convert=True)

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
            return VectorSpace(self.base_ring(), self.ambient_dim() + 1)
        from sage.modules.free_module import FreeModule
        return FreeModule(self.base_ring(), self.ambient_dim() + 1)

    def _repr_base_ring(self):
        """
        Return an abbreviated string representation of the base ring.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: Polyhedra(QQ, 3)._repr_base_ring()
            'QQ'
            sage: x = polygen(ZZ, 'x')
            sage: K.<sqrt3> = NumberField(x^2 - 3, embedding=AA(3).sqrt())              # needs sage.rings.number_field
            sage: Polyhedra(K, 4)._repr_base_ring()                                     # needs sage.rings.number_field
            '(Number Field in sqrt3 with defining polynomial x^2 - 3 with sqrt3 = 1.732050807568878?)'
        """

        if self.base_ring() is ZZ:
            return 'ZZ'
        if self.base_ring() is QQ:
            return 'QQ'
        if self.base_ring() is RDF:
            return 'RDF'
        try:
            from sage.rings.qqbar import AA
        except ImportError:
            pass
        else:
            if self.base_ring() is AA:
                return 'AA'
        return '({0})'.format(self.base_ring())

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
            sage: x = polygen(ZZ, 'x')
            sage: K.<sqrt3> = NumberField(x^2 - 3, embedding=AA(3).sqrt())              # needs sage.rings.number_field
            sage: Polyhedra(K, 4)._repr_ambient_module()                                # needs sage.rings.number_field
            '(Number Field in sqrt3 with defining polynomial x^2 - 3 with sqrt3 = 1.732050807568878?)^4'
        """
        s = self._repr_base_ring()
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
        return 'Polyhedra in ' + self._repr_ambient_module()

    def _element_constructor_(self, *args, **kwds):
        """
        The element (polyhedron) constructor.

        INPUT:

        - ``Vrep`` -- a list ``[vertices, rays, lines]`` or ``None``.

        - ``Hrep`` -- a list ``[ieqs, eqns]`` or ``None``.

        - ``convert`` -- boolean keyword argument (default:
          ``True``). Whether to convert the coordinates into the base
          ring.

        - ``**kwds`` -- optional remaining keywords that are passed to the
          polyhedron constructor.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: P = Polyhedra(QQ, 3)
            sage: P._element_constructor_([[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0,0,1)], [], []], None)
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 4 vertices
            sage: P([[(0,0,0),(1,0,0),(0,1,0),(0,0,1)], [], []], None)
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 4 vertices
            sage: P(0)
            A 0-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex

        Check that :trac:`21270` is fixed::

            sage: # needs sage.rings.number_field
            sage: poly = polytopes.regular_polygon(7)
            sage: lp, x = poly.to_linear_program(solver='InteractiveLP',
            ....:                                return_variable=True)
            sage: lp.set_objective(x[0] + x[1])
            sage: b = lp.get_backend()
            sage: P = b.interactive_lp_problem()
            sage: p = P.plot()                                                          # needs sage.plot

            sage: Q = Polyhedron(ieqs=[[-499999, 1000000], [1499999, -1000000]])
            sage: P = Polyhedron(ieqs=[[0, 1.0], [1.0, -1.0]], base_ring=RDF)
            sage: Q.intersection(P)
            A 1-dimensional polyhedron in RDF^1 defined as the convex hull of 2 vertices
            sage: P.intersection(Q)
            A 1-dimensional polyhedron in RDF^1 defined as the convex hull of 2 vertices

        The default is not to copy an object if the parent is ``self``::

            sage: p = polytopes.cube(backend='field')
            sage: P = p.parent()
            sage: q = P._element_constructor_(p)
            sage: q is p
            True
            sage: r = P._element_constructor_(p, copy=True)
            sage: r is p
            False

        When the parent of the object is not ``self``, the default is not to copy::

            sage: # needs sage.rings.number_field
            sage: Q = P.base_extend(AA)
            sage: q = Q._element_constructor_(p)
            sage: q is p
            False
            sage: q = Q._element_constructor_(p, copy=False)
            Traceback (most recent call last):
            ...
            ValueError: you need to make a copy when changing the parent

        For mutable polyhedra either ``copy`` or ``mutable`` must be specified::

            sage: p = Polyhedron(vertices=[[0, 1], [1, 0]], mutable=True)
            sage: P = p.parent()
            sage: q = P._element_constructor_(p)
            Traceback (most recent call last):
            ...
            ValueError: must make a copy to obtain immutable object from mutable input
            sage: q = P._element_constructor_(p, mutable=True)
            sage: q is p
            True
            sage: r = P._element_constructor_(p, copy=True)
            sage: r.is_mutable()
            False
            sage: r is p
            False
        """
        nargs = len(args)
        convert = kwds.pop('convert', True)

        def convert_base_ring(lstlst):
            return [[self.base_ring()(x) for x in lst] for lst in lstlst]

        # renormalize before converting when going from QQ to RDF, see trac 21270
        def convert_base_ring_Hrep(lstlst):
            newlstlst = []
            for lst in lstlst:
                if all(c in QQ for c in lst):
                    m = max(abs(w) for w in lst)
                    if m == 0:
                        newlstlst.append(lst)
                    else:
                        newlstlst.append([q / m for q in lst])
                else:
                    newlstlst.append(lst)
            return convert_base_ring(newlstlst)
        if nargs == 2:
            Vrep, Hrep = args
            if convert and Hrep:
                if self.base_ring == RDF:
                    Hrep = [convert_base_ring_Hrep(_) for _ in Hrep]
                else:
                    Hrep = [convert_base_ring(_) for _ in Hrep]
            if convert and Vrep:
                Vrep = [convert_base_ring(_) for _ in Vrep]
            return self.element_class(self, Vrep, Hrep, **kwds)
        if nargs == 1 and is_Polyhedron(args[0]):
            copy = kwds.pop('copy', args[0].parent() is not self)
            mutable = kwds.pop('mutable', False)

            if not copy and args[0].parent() is not self:
                raise ValueError("you need to make a copy when changing the parent")
            if args[0].is_mutable() and not copy and not mutable:
                raise ValueError("must make a copy to obtain immutable object from mutable input")
            if not copy and mutable is args[0].is_mutable():
                return args[0]

            polyhedron = args[0]
            return self._element_constructor_polyhedron(polyhedron, mutable=mutable, **kwds)
        if nargs == 1 and args[0] == 0:
            return self.zero()
        raise ValueError('cannot convert to polyhedron object')

    def _element_constructor_polyhedron(self, polyhedron, **kwds):
        """
        The element (polyhedron) constructor for the case of 1 argument, a polyhedron.

        Set up the element using both representations,
        if the backend can handle it.

        Otherwise set up the element from Hrepresentation.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: P = Polyhedra(QQ, 3, backend='cdd')
            sage: p = Polyhedron(vertices=[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
            sage: p
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices
            sage: P(p)
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 4 vertices

            sage: P = Polyhedra(AA, 3, backend='field')                                 # needs sage.rings.number_field
            sage: vertices = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
            sage: p = Polyhedron(vertices=vertices)
            sage: P(p)                                                                  # needs sage.rings.number_field
            A 3-dimensional polyhedron in AA^3 defined as the convex hull of 4 vertices
        """
        Vrep = None
        if hasattr(self.Element, '_init_from_Vrepresentation_and_Hrepresentation'):
            Vrep = [polyhedron.vertex_generator(), polyhedron.ray_generator(),
                    polyhedron.line_generator()]
        Hrep = [polyhedron.inequality_generator(), polyhedron.equation_generator()]
        return self._element_constructor_(Vrep, Hrep, Vrep_minimal=True, Hrep_minimal=True, **kwds)

    def base_extend(self, base_ring, backend=None, ambient_dim=None):
        """
        Return the base extended parent.

        INPUT:

        - ``base_ring``, ``backend`` -- see
          :func:`~sage.geometry.polyhedron.constructor.Polyhedron`.
        - ``ambient_dim`` -- if not ``None`` change ambient dimension
          accordingly.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: Polyhedra(ZZ, 3).base_extend(QQ)
            Polyhedra in QQ^3
            sage: Polyhedra(ZZ, 3).an_element().base_extend(QQ)
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 4 vertices
            sage: Polyhedra(QQ, 2).base_extend(ZZ)
            Polyhedra in QQ^2

        TESTS:

        Test that :trac:`22575` is fixed::

            sage: P = Polyhedra(ZZ,3).base_extend(QQ, backend='field')
            sage: P.backend()
            'field'
        """
        if self.base_ring().has_coerce_map_from(base_ring):
            new_ring = self.base_ring()
        else:
            new_ring = self._coerce_base_ring(base_ring)

        return self.change_ring(new_ring, backend=backend, ambient_dim=ambient_dim)

    def change_ring(self, base_ring, backend=None, ambient_dim=None):
        """
        Return the parent with the new base ring.

        INPUT:

        - ``base_ring``, ``backend`` -- see
          :func:`~sage.geometry.polyhedron.constructor.Polyhedron`.
        - ``ambient_dim`` -- if not ``None`` change ambient dimension
          accordingly.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: Polyhedra(ZZ, 3).change_ring(QQ)
            Polyhedra in QQ^3
            sage: Polyhedra(ZZ, 3).an_element().change_ring(QQ)
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 4 vertices

            sage: Polyhedra(RDF, 3).change_ring(QQ).backend()
            'cdd'
            sage: Polyhedra(QQ, 3).change_ring(ZZ, ambient_dim=4)
            Polyhedra in ZZ^4
            sage: Polyhedra(QQ, 3, backend='cdd').change_ring(QQ, ambient_dim=4).backend()
            'cdd'
        """
        if ambient_dim is None:
            ambient_dim = self.ambient_dim()

        if base_ring == self.base_ring() and \
                ambient_dim == self.ambient_dim() and \
                (backend is None or backend == self.backend()):
            return self

        # if not specified try the same backend
        if backend is None and does_backend_handle_base_ring(base_ring, self.backend()):
            return Polyhedra(base_ring, ambient_dim, backend=self.backend())

        return Polyhedra(base_ring, ambient_dim, backend=backend)

    def _coerce_base_ring(self, other):
        r"""
        Return the common base ring for both ``self`` and ``other``.

        This method is not part of the coercion framework, but only a
        convenience function for :class:`Polyhedra_base`.

        INPUT:

        - ``other`` -- must be either:

          * another ``Polyhedron`` object

          * `\ZZ`, `\QQ`, `RDF`, or a ring that can be coerced into them.

          * a constant that can be coerced to `\ZZ`, `\QQ`, or `RDF`.

        OUTPUT:

        Either `\ZZ`, `\QQ`, or `RDF`. Raises :class:`TypeError` if
        ``other`` is not a suitable input.

        .. NOTE::

            "Real" numbers in sage are not necessarily elements of
            ``RDF``. For example, the literal `1.0` is not.

        EXAMPLES::

            sage: triangle_QQ  = Polyhedron(vertices=[[1,0],[0,1],[1,1]], base_ring=QQ).parent()
            sage: triangle_RDF = Polyhedron(vertices=[[1,0],[0,1],[1,1]], base_ring=RDF).parent()
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

        TESTS:

        Test that :trac:`28770` is fixed::

            sage: z = QQ['z'].0
            sage: K = NumberField(z^2 - 2, 's')                                         # needs sage.rings.number_field
            sage: triangle_QQ._coerce_base_ring(K)                                      # needs sage.rings.number_field
            Number Field in s with defining polynomial z^2 - 2
            sage: triangle_QQ._coerce_base_ring(K.gen())                                # needs sage.rings.number_field
            Number Field in s with defining polynomial z^2 - 2

            sage: z = QQ['z'].0
            sage: K = NumberField(z^2 - 2, 's')                                         # needs sage.rings.number_field
            sage: K.gen() * polytopes.simplex(backend='field')                          # needs sage.rings.number_field
            A 3-dimensional polyhedron in
             (Number Field in s with defining polynomial z^2 - 2)^4
             defined as the convex hull of 4 vertices
        """
        from sage.structure.element import Element
        if isinstance(other, Element):
            other = other.parent()
        if hasattr(other, "is_ring") and other.is_ring():
            other_ring = other
        else:
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
                        raise TypeError(f'Could not coerce {other} into ZZ, QQ, or RDF.')

        if not other_ring.is_exact():
            other_ring = RDF  # the only supported floating-point numbers for now

        cm_map, cm_ring = get_coercion_model().analyse(self.base_ring(), other_ring)
        if cm_ring is None:
            raise TypeError(f'Could not coerce type {other} into ZZ, QQ, or RDF.')
        return cm_ring

    def _coerce_map_from_(self, X):
        r"""
        Return whether there is a coercion from ``X``.

        INPUT:

        - ``X`` -- anything.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: Polyhedra(QQ, 3).has_coerce_map_from(Polyhedra(ZZ, 3))   # indirect doctest
            True
            sage: Polyhedra(ZZ, 3).has_coerce_map_from(Polyhedra(QQ, 3))
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
            with precomposition on left by Coercion map:
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
            with precomposition on left by
             Identity endomorphism of Polyhedra in ZZ^2
            with precomposition on right by
             Generic endomorphism of Ambient free module of rank 2 over the principal ideal domain Integer Ring

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

        if op is operator.mul and other in Rings().Commutative():
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

            sage: p = Polyhedron([(1,2,3), (2/3,3/4,4/5)])   # indirect doctest
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

            sage: p = Polyhedron([(1,2,3), (2/3,3/4,4/5)])   # indirect doctest
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

            sage: p = Polyhedron([(1,2,3), (2/3,3/4,4/5)], rays=[(5/6,6/7,7/8)])   # indirect doctest
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

            sage: p = Polyhedron([(1,2,3), (2/3,3/4,4/5)], rays=[(5/6,6/7,7/8)])   # indirect doctest
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

            sage: p = Polyhedron([(1,2,3), (2/3,3/4,4/5)], lines=[(5/6,6/7,7/8)])   # indirect doctest
            sage: next(p.line_generator())
            A line in the direction (140, 144, 147)
        """
        try:
            obj = self._Line_pool.pop()
        except IndexError:
            obj = Line(self)
        obj._set_data(polyhedron, data)
        return obj


from sage.geometry.polyhedron.backend_cdd import Polyhedron_QQ_cdd
lazy_import('sage.geometry.polyhedron.backend_cdd_rdf', 'Polyhedron_RDF_cdd')
from sage.geometry.polyhedron.backend_ppl import Polyhedron_ZZ_ppl, Polyhedron_QQ_ppl
from sage.geometry.polyhedron.backend_normaliz import Polyhedron_normaliz, Polyhedron_ZZ_normaliz, Polyhedron_QQ_normaliz
from sage.geometry.polyhedron.backend_polymake import Polyhedron_polymake
from sage.geometry.polyhedron.backend_field import Polyhedron_field
from sage.geometry.polyhedron.backend_number_field import Polyhedron_number_field


class Polyhedra_ZZ_ppl(Polyhedra_base):
    Element = Polyhedron_ZZ_ppl

    def _element_constructor_polyhedron(self, polyhedron, **kwds):
        """
        The element (polyhedron) constructor for the case of 1 argument, a polyhedron.

        Set up with the ``ppl_polyhedron`` of ``self``, if available.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: P = Polyhedra(ZZ, 3)
            sage: p = Polyhedron(vertices=[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)], base_ring=QQ)
            sage: p
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 4 vertices
            sage: P(p)
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices

            sage: p = Polyhedron(vertices=[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)], backend='cdd')
            sage: P(p)
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices
        """
        from copy import copy
        if polyhedron.backend() == "ppl":
            return self._element_constructor_(None, None, ppl_polyhedron=copy(polyhedron._ppl_polyhedron), **kwds)
        else:
            return Polyhedra_base._element_constructor_polyhedron(self, polyhedron, **kwds)


class Polyhedra_ZZ_normaliz(Polyhedra_base):
    Element = Polyhedron_ZZ_normaliz


class Polyhedra_QQ_ppl(Polyhedra_base):
    Element = Polyhedron_QQ_ppl

    def _element_constructor_polyhedron(self, polyhedron, **kwds):
        """
        The element (polyhedron) constructor for the case of 1 argument, a polyhedron.

        Set up with the ``ppl_polyhedron`` of ``self``, if available.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: P = Polyhedra(QQ, 3)
            sage: p = Polyhedron(vertices=[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
            sage: p
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices
            sage: P(p)
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 4 vertices

            sage: p = Polyhedron(vertices=[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)], backend='cdd')
            sage: P(p)
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 4 vertices
        """
        from copy import copy
        if polyhedron.backend() == "ppl":
            return self._element_constructor_(None, None, ppl_polyhedron=copy(polyhedron._ppl_polyhedron), **kwds)
        else:
            return Polyhedra_base._element_constructor_polyhedron(self, polyhedron, **kwds)


class Polyhedra_QQ_normaliz(Polyhedra_base):
    Element = Polyhedron_QQ_normaliz


class Polyhedra_QQ_cdd(Polyhedra_base):
    Element = Polyhedron_QQ_cdd


class Polyhedra_RDF_cdd(Polyhedra_base):
    Element = Polyhedron_RDF_cdd


class Polyhedra_normaliz(Polyhedra_base):
    Element = Polyhedron_normaliz


class Polyhedra_polymake(Polyhedra_base):
    Element = Polyhedron_polymake


class Polyhedra_field(Polyhedra_base):
    Element = Polyhedron_field


class Polyhedra_number_field(Polyhedra_base):
    Element = Polyhedron_number_field


@cached_function
def does_backend_handle_base_ring(base_ring, backend):
    r"""
    Return true, if ``backend`` can handle ``base_ring``.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.parent import does_backend_handle_base_ring
        sage: does_backend_handle_base_ring(QQ, 'ppl')
        True
        sage: does_backend_handle_base_ring(QQ[sqrt(5)], 'ppl')                         # needs sage.rings.number_field sage.symbolic
        False
        sage: does_backend_handle_base_ring(QQ[sqrt(5)], 'field')                       # needs sage.rings.number_field sage.symbolic
        True
    """
    try:
        Polyhedra(base_ring, 0, backend)
    except ValueError:
        return False
    return True
