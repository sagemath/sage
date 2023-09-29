r"""
Toric lattices

This module was designed as a part of the framework for toric varieties
(:mod:`~sage.schemes.toric.variety`,
:mod:`~sage.schemes.toric.fano_variety`).

All toric lattices are isomorphic to `\ZZ^n` for some `n`, but will prevent
you from doing "wrong" operations with objects from different lattices.

AUTHORS:

- Andrey Novoseltsev (2010-05-27): initial version.
- Andrey Novoseltsev (2010-07-30): sublattices and quotients.

EXAMPLES:

The simplest way to create a toric lattice is to specify its dimension only::

    sage: N = ToricLattice(3)
    sage: N
    3-d lattice N

While our lattice ``N`` is called exactly "N" it is a coincidence: all
lattices are called "N" by default::

    sage: another_name = ToricLattice(3)
    sage: another_name
    3-d lattice N

If fact, the above lattice is exactly the same as before as an object in
memory::

    sage: N is another_name
    True

There are actually four names associated to a toric lattice and they all must
be the same for two lattices to coincide::

    sage: N, N.dual(), latex(N), latex(N.dual())
    (3-d lattice N, 3-d lattice M, N, M)

Notice that the lattice dual to ``N`` is called "M" which is standard in toric
geometry. This happens only if you allow completely automatic handling of
names::

    sage: another_N = ToricLattice(3, "N")
    sage: another_N.dual()
    3-d lattice N*
    sage: N is another_N
    False

What can you do with toric lattices? Well, their main purpose is to allow
creation of elements of toric lattices::

    sage: n = N([1,2,3])
    sage: n
    N(1, 2, 3)
    sage: M = N.dual()
    sage: m = M(1,2,3)
    sage: m
    M(1, 2, 3)

Dual lattices can act on each other::

    sage: n * m
    14
    sage: m * n
    14

You can also add elements of the same lattice or scale them::

    sage: 2 * n
    N(2, 4, 6)
    sage: n * 2
    N(2, 4, 6)
    sage: n + n
    N(2, 4, 6)

However, you cannot "mix wrong lattices" in your expressions::

    sage: n + m
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for +:
    '3-d lattice N' and '3-d lattice M'
    sage: n * n
    Traceback (most recent call last):
    ...
    TypeError: elements of the same toric lattice cannot be multiplied!
    sage: n == m
    False

Note that ``n`` and ``m`` are not equal to each other even though they are
both "just (1,2,3)." Moreover, you cannot easily convert elements between
toric lattices::

    sage: M(n)
    Traceback (most recent call last):
    ...
    TypeError: N(1, 2, 3) cannot be converted to 3-d lattice M!

If you really need to consider elements of one lattice as elements of another,
you can either use intermediate conversion to "just a vector"::

    sage: ZZ3 = ZZ^3
    sage: n_in_M = M(ZZ3(n))
    sage: n_in_M
    M(1, 2, 3)
    sage: n == n_in_M
    False
    sage: n_in_M == m
    True

Or you can create a homomorphism from one lattice to any other::

    sage: h = N.hom(identity_matrix(3), M)
    sage: h(n)
    M(1, 2, 3)

.. WARNING::

    While integer vectors (elements of `\ZZ^n`) are printed as ``(1,2,3)``,
    in the code ``(1,2,3)`` is a :class:`tuple`, which has nothing to do
    neither with vectors, nor with toric lattices, so the following is
    probably not what you want while working with toric geometry objects::

        sage: (1,2,3) + (1,2,3)
        (1, 2, 3, 1, 2, 3)

    Instead, use syntax like ::

        sage: N(1,2,3) + N(1,2,3)
        N(2, 4, 6)
"""
# Parts of the "tutorial" above are also in toric_lattice_element.pyx.


# ****************************************************************************
#       Copyright (C) 2010 Andrey Novoseltsev <novoselt@gmail.com>
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.geometry.toric_lattice_element import (ToricLatticeElement,
                                                 is_ToricLatticeElement)
from sage.geometry.toric_plotter import ToricPlotter
from sage.misc.latex import latex
from sage.structure.all import parent
from sage.structure.richcmp import (richcmp_method, richcmp, rich_to_bool,
                                    richcmp_not_equal)
from sage.modules.fg_pid.fgp_element import FGP_Element
from sage.modules.fg_pid.fgp_module import FGP_Module_class
from sage.modules.free_module import (FreeModule_ambient_pid,
                                      FreeModule_generic_pid,
                                      FreeModule_submodule_pid,
                                      FreeModule_submodule_with_basis_pid)
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.structure.factory import UniqueFactory


def is_ToricLattice(x):
    r"""
    Check if ``x`` is a toric lattice.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    - ``True`` if ``x`` is a toric lattice and ``False`` otherwise.

    EXAMPLES::

        sage: from sage.geometry.toric_lattice import (
        ....:   is_ToricLattice)
        sage: is_ToricLattice(1)
        False
        sage: N = ToricLattice(3)
        sage: N
        3-d lattice N
        sage: is_ToricLattice(N)
        True
    """
    return isinstance(x, ToricLattice_generic)


def is_ToricLatticeQuotient(x):
    r"""
    Check if ``x`` is a toric lattice quotient.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    - ``True`` if ``x`` is a toric lattice quotient and ``False`` otherwise.

    EXAMPLES::

        sage: from sage.geometry.toric_lattice import (
        ....:   is_ToricLatticeQuotient)
        sage: is_ToricLatticeQuotient(1)
        False
        sage: N = ToricLattice(3)
        sage: N
        3-d lattice N
        sage: is_ToricLatticeQuotient(N)
        False
        sage: Q = N / N.submodule([(1,2,3), (3,2,1)])
        sage: Q
        Quotient with torsion of 3-d lattice N
        by Sublattice <N(1, 2, 3), N(0, 4, 8)>
        sage: is_ToricLatticeQuotient(Q)
        True
    """
    return isinstance(x, ToricLattice_quotient)


class ToricLatticeFactory(UniqueFactory):
    r"""
    Create a lattice for toric geometry objects.

    INPUT:

    - ``rank`` -- nonnegative integer, the only mandatory parameter;

    - ``name`` -- string;

    - ``dual_name`` -- string;

    - ``latex_name`` -- string;

    - ``latex_dual_name`` -- string.

    OUTPUT:

    - lattice.

    A toric lattice is uniquely determined by its rank and associated names.
    There are four such "associated names" whose meaning should be clear from
    the names of the corresponding parameters, but the choice of default
    values is a little bit involved. So here is the full description of the
    "naming algorithm":

    #. If no names were given at all, then this lattice will be called "N" and
       the dual one "M". These are the standard choices in toric geometry.

    #. If ``name`` was given and ``dual_name`` was not, then ``dual_name``
       will be ``name`` followed by "*".

    #. If LaTeX names were not given, they will coincide with the "usual"
       names, but if ``dual_name`` was constructed automatically, the trailing
       star will be typeset as a superscript.

    EXAMPLES:

    Let's start with no names at all and see how automatic names are given::

        sage: L1 = ToricLattice(3)
        sage: L1
        3-d lattice N
        sage: L1.dual()
        3-d lattice M

    If we give the name "N" explicitly, the dual lattice will be called "N*"::

        sage: L2 = ToricLattice(3, "N")
        sage: L2
        3-d lattice N
        sage: L2.dual()
        3-d lattice N*

    However, we can give an explicit name for it too::

        sage: L3 = ToricLattice(3, "N", "M")
        sage: L3
        3-d lattice N
        sage: L3.dual()
        3-d lattice M

    If you want, you may also give explicit LaTeX names::

        sage: L4 = ToricLattice(3, "N", "M", r"\mathbb{N}", r"\mathbb{M}")
        sage: latex(L4)
        \mathbb{N}
        sage: latex(L4.dual())
        \mathbb{M}

    While all four lattices above are called "N", only two of them are equal
    (and are actually the same)::

        sage: L1 == L2
        False
        sage: L1 == L3
        True
        sage: L1 is L3
        True
        sage: L1 == L4
        False

    The reason for this is that ``L2`` and ``L4`` have different names either
    for dual lattices or for LaTeX typesetting.
    """

    def create_key(self, rank, name=None, dual_name=None,
                   latex_name=None, latex_dual_name=None):
        """
        Create a key that uniquely identifies this toric lattice.

        See :class:`ToricLattice <ToricLatticeFactory>` for documentation.

        .. WARNING::

            You probably should not use this function directly.

        TESTS::

            sage: ToricLattice.create_key(3)
            (3, 'N', 'M', 'N', 'M')
            sage: N = ToricLattice(3)
            sage: loads(dumps(N)) is N
            True
            sage: TestSuite(N).run()
        """
        rank = int(rank)
        # Should we use standard M and N lattices?
        if name is None:
            if dual_name is not None:
                raise ValueError("you can name the dual lattice only if you "
                                 "also name the original one!")
            name = "N"
            dual_name = "M"
        if latex_name is None:
            latex_name = name
        # Now name and latex_name are set
        # The default for latex_dual_name depends on whether dual_name was
        # given or constructed, so we determine it before dual_name
        if latex_dual_name is None:
            latex_dual_name = (dual_name if dual_name is not None
                                         else latex_name + "^*")
        if dual_name is None:
            dual_name = name + "*"
        return (rank, name, dual_name, latex_name, latex_dual_name)

    def create_object(self, version, key):
        r"""
        Create the toric lattice described by ``key``.

        See :class:`ToricLattice <ToricLatticeFactory>` for documentation.

        .. WARNING::

            You probably should not use this function directly.

        TESTS::

            sage: key = ToricLattice.create_key(3)
            sage: ToricLattice.create_object(1, key)
            3-d lattice N
        """
        return ToricLattice_ambient(*key)


ToricLattice = ToricLatticeFactory("ToricLattice")


# Possible TODO's:
# - implement a better construction() method, which still will prohibit
#   operations mixing lattices by conversion to ZZ^n
# - maybe __call__ is not the right place to prohibit conversion between
#   lattices (we need it now so that morphisms behave nicely)
class ToricLattice_generic(FreeModule_generic_pid):
    r"""
    Abstract base class for toric lattices.
    """

    Element = ToricLatticeElement

    # It is not recommended to override __call__ in Parent-derived objects
    # since it may interfere with the coercion model. We do it here to allow
    # N(1,2,3) to be interpreted as N([1,2,3]). We also prohibit N(m) where
    # m is an element of another lattice. Otherwise morphisms will care only
    # about dimension of lattices.
    def __call__(self, *args, **kwds):
        r"""
        Construct a new element of ``self``.

        INPUT:

        - anything that can be interpreted as coordinates, except for elements
          of other lattices.

        OUTPUT:

        - :class:`~sage.geometry.toric_lattice_element.ToricLatticeElement`.

        TESTS::

            sage: N = ToricLattice(3)
            sage: N.__call__([1,2,3])
            N(1, 2, 3)
            sage: N([1,2,3])    # indirect test
            N(1, 2, 3)

        The point of overriding this function was to allow writing the above
        command as::

            sage: N(1,2,3)
            N(1, 2, 3)

        We also test that the special treatment of zero still works::

            sage: N(0)
            N(0, 0, 0)

        Quotients of toric lattices can be converted to a new toric
        lattice of the appropriate dimension::

            sage: N3 = ToricLattice(3, 'N3')
            sage: Q = N3 / N3.span([ N3(1,2,3) ])
            sage: Q.an_element()
            N3[1, 0, 0]
            sage: N2 = ToricLattice(2, 'N2')
            sage: N2( Q.an_element() )
            N2(1, 0)
        """
        supercall = super().__call__
        if args == (0, ):
            # Special treatment for N(0) to return (0,...,0)
            return supercall(*args, **kwds)

        if (isinstance(args[0], ToricLattice_quotient_element)
            and args[0].parent().is_torsion_free()):
            # convert a torsion free quotient lattice
            return supercall(list(args[0]), **kwds)

        try:
            coordinates = [ZZ(_) for _ in args]
        except TypeError:
            # Prohibit conversion of elements of other lattices
            if (is_ToricLatticeElement(args[0])
                and args[0].parent().ambient_module()
                is not self.ambient_module()):
                raise TypeError("%s cannot be converted to %s!"
                                % (args[0], self))
            # "Standard call"
            return supercall(*args, **kwds)
        # Coordinates were given without packing them into a list or a tuple
        return supercall(coordinates, **kwds)

    def _coerce_map_from_(self, other):
        """
        Return a coercion map from ``other`` to ``self``, or None.

        This prevents the construction of coercion maps between
        lattices with different ambient modules, so :meth:`__call__`
        is invoked instead, which prohibits conversion::

            sage: N = ToricLattice(3)
            sage: M = N.dual()
            sage: M(N(1,2,3))
            Traceback (most recent call last):
            ...
            TypeError: N(1, 2, 3) cannot be converted to 3-d lattice M!

        """
        if (is_ToricLattice(other) and
            other.ambient_module() is not self.ambient_module()):
            return None
        return super()._convert_map_from_(other)

    def __contains__(self, point):
        r"""
        Check if ``point`` is an element of ``self``.

        INPUT:

        - ``point`` -- anything.

        OUTPUT:

        - ``True`` if ``point`` is an element of ``self``, ``False``
          otherwise.

        TESTS::

            sage: N = ToricLattice(3)
            sage: M = N.dual()
            sage: L = ToricLattice(3, "L")
            sage: 1 in N
            False
            sage: (1,0) in N
            False
            sage: (1,0,0) in N
            True
            sage: N(1,0,0) in N
            True
            sage: M(1,0,0) in N
            False
            sage: L(1,0,0) in N
            False
            sage: (1/2,0,0) in N
            False
            sage: (2/2,0,0) in N
            True
        """
        try:
            self(point)
        except TypeError:
            return False
        return True

    # We need to override this function, otherwise e.g. the sum of elements of
    # different lattices of the same dimension will live in ZZ^n.
    def construction(self):
        r"""
        Return the functorial construction of ``self``.

        OUTPUT:

        - ``None``, we do not think of toric lattices as constructed from
          simpler objects since we do not want to perform arithmetic involving
          different lattices.

        TESTS::

            sage: print(ToricLattice(3).construction())
            None
        """
        return None

    def direct_sum(self, other):
        r"""
        Return the direct sum with ``other``.

        INPUT:

        - ``other`` -- a toric lattice or more general module.

        OUTPUT:

        The direct sum of ``self`` and ``other`` as `\ZZ`-modules. If
        ``other`` is a :class:`ToricLattice <ToricLatticeFactory>`,
        another toric lattice will be returned.

        EXAMPLES::

            sage: K = ToricLattice(3, 'K')
            sage: L = ToricLattice(3, 'L')
            sage: N = K.direct_sum(L); N
            6-d lattice K+L
            sage: N, N.dual(), latex(N), latex(N.dual())
            (6-d lattice K+L, 6-d lattice K*+L*, K \oplus L, K^* \oplus L^*)

        With default names::

            sage: N = ToricLattice(3).direct_sum(ToricLattice(2))
            sage: N, N.dual(), latex(N), latex(N.dual())
            (5-d lattice N+N, 5-d lattice M+M, N \oplus N, M \oplus M)

        If ``other`` is not a :class:`ToricLattice
        <ToricLatticeFactory>`, fall back to sum of modules::

            sage: ToricLattice(3).direct_sum(ZZ^2)
            Free module of degree 5 and rank 5 over Integer Ring
            Echelon basis matrix:
            [1 0 0 0 0]
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
        """
        if not isinstance(other, ToricLattice_generic):
            return super().direct_sum(other)

        def make_name(N1, N2, use_latex=False):
            if use_latex:
                return latex(N1)+ r' \oplus ' +latex(N2)
            else:
                return N1._name+ '+' +N2._name

        rank = self.rank() + other.rank()
        name = make_name(self, other, False)
        dual_name = make_name(self.dual(), other.dual(), False)
        latex_name = make_name(self, other, True)
        latex_dual_name = make_name(self.dual(), other.dual(), True)
        return ToricLattice(rank, name, dual_name, latex_name, latex_dual_name)

    def intersection(self, other):
        r"""
        Return the intersection of ``self`` and ``other``.

        INPUT:

        - ``other`` - a toric (sub)lattice.dual

        OUTPUT:

        - a toric (sub)lattice.

        EXAMPLES::

            sage: N = ToricLattice(3)
            sage: Ns1 = N.submodule([N(2,4,0), N(9,12,0)])
            sage: Ns2 = N.submodule([N(1,4,9), N(9,2,0)])
            sage: Ns1.intersection(Ns2)
            Sublattice <N(54, 12, 0)>

        Note that if one of the intersecting sublattices is a sublattice of
        another, no new lattices will be constructed::

            sage: N.intersection(N) is N
            True
            sage: Ns1.intersection(N) is Ns1
            True
            sage: N.intersection(Ns1) is Ns1
            True
        """
        # Lattice-specific input check
        if not is_ToricLattice(other):
            raise TypeError("%s is not a toric lattice!" % other)
        if self.ambient_module() != other.ambient_module():
            raise ValueError("%s and %s have different ambient lattices!" %
                             (self, other))
        # Construct a generic intersection, but make sure to return a lattice.
        I = super().intersection(other)
        if not is_ToricLattice(I):
            I = self.ambient_module().submodule(I.basis())
        return I

    def quotient(self, sub, check=True,
                 positive_point=None, positive_dual_point=None, **kwds):
        """
        Return the quotient of ``self`` by the given sublattice ``sub``.

        INPUT:

        - ``sub`` -- sublattice of self;

        - ``check`` -- (default: True) whether or not to check that ``sub`` is
          a valid sublattice.

        If the quotient is one-dimensional and torsion free, the
        following two mutually exclusive keyword arguments are also
        allowed. They decide the sign choice for the (single)
        generator of the quotient lattice:

        - ``positive_point`` -- a lattice point of ``self`` not in the
          sublattice ``sub`` (that is, not zero in the quotient
          lattice). The quotient generator will be in the same
          direction as ``positive_point``.

        - ``positive_dual_point`` -- a dual lattice point. The
          quotient generator will be chosen such that its lift has a
          positive product with ``positive_dual_point``. Note: if
          ``positive_dual_point`` is not zero on the sublattice
          ``sub``, then the notion of positivity will depend on the
          choice of lift!

        Further named arguments are passed to the constructor of a toric lattice
        quotient.

        EXAMPLES::

            sage: N = ToricLattice(3)
            sage: Ns = N.submodule([N(2,4,0), N(9,12,0)])
            sage: Q = N/Ns
            sage: Q
            Quotient with torsion of 3-d lattice N
            by Sublattice <N(1, 8, 0), N(0, 12, 0)>

        Attempting to quotient one lattice by a sublattice of another
        will result in a ``ValueError``::

            sage: N = ToricLattice(3)
            sage: M = ToricLattice(3, name='M')
            sage: Ms = M.submodule([M(2,4,0), M(9,12,0)])
            sage: N.quotient(Ms)
            Traceback (most recent call last):
            ...
            ValueError: M(1, 8, 0) cannot generate a sublattice of
            3-d lattice N

        However, if we forget the sublattice structure, then it is
        possible to quotient by vector spaces or modules constructed
        from any sublattice::

            sage: N = ToricLattice(3)
            sage: M = ToricLattice(3, name='M')
            sage: Ms = M.submodule([M(2,4,0), M(9,12,0)])
            sage: N.quotient(Ms.vector_space())
            Quotient with torsion of 3-d lattice N by Sublattice
            <N(1, 8, 0), N(0, 12, 0)>
            sage: N.quotient(Ms.sparse_module())
            Quotient with torsion of 3-d lattice N by Sublattice
            <N(1, 8, 0), N(0, 12, 0)>

        See :class:`ToricLattice_quotient` for more examples.

        TESTS:

        We check that :trac:`19603` is fixed::

            sage: K = Cone([(1,0,0),(0,1,0)])
            sage: K.lattice()
            3-d lattice N
            sage: K.orthogonal_sublattice()
            Sublattice <M(0, 0, 1)>
            sage: K.lattice().quotient(K.orthogonal_sublattice())
            Traceback (most recent call last):
            ...
            ValueError: M(0, 0, 1) cannot generate a sublattice of
            3-d lattice N

        We can quotient by the trivial sublattice::

            sage: N = ToricLattice(3)
            sage: N.quotient(N.zero_submodule())
            3-d lattice, quotient of 3-d lattice N by Sublattice <>

        We can quotient a lattice by itself::

            sage: N = ToricLattice(3)
            sage: N.quotient(N)
            0-d lattice, quotient of 3-d lattice N by Sublattice
            <N(1, 0, 0), N(0, 1, 0), N(0, 0, 1)>
        """
        return ToricLattice_quotient(self, sub, check,
                                     positive_point, positive_dual_point, **kwds)

    def saturation(self):
        r"""
        Return the saturation of ``self``.

        OUTPUT:

        - a :class:`toric lattice <ToricLatticeFactory>`.

        EXAMPLES::

            sage: N = ToricLattice(3)
            sage: Ns = N.submodule([(1,2,3), (4,5,6)])
            sage: Ns
            Sublattice <N(1, 2, 3), N(0, 3, 6)>
            sage: Ns_sat = Ns.saturation()
            sage: Ns_sat
            Sublattice <N(1, 0, -1), N(0, 1, 2)>
            sage: Ns_sat is Ns_sat.saturation()
            True
        """
        S = super().saturation()
        return S if is_ToricLattice(S) else self.ambient_module().submodule(S)

    def span(self, gens, base_ring=ZZ, *args, **kwds):
        r"""
        Return the span of the given generators.

        INPUT:

        - ``gens`` -- list of elements of the ambient vector space of
          ``self``.

        - ``base_ring`` -- (default: `\ZZ`) base ring for the generated module.

        OUTPUT:

        - submodule spanned by ``gens``.

        .. NOTE::

            The output need not be a submodule of ``self``, nor even of the
            ambient space. It must, however, be contained in the ambient
            vector space.

        See also :meth:`span_of_basis`,
        :meth:`~sage.modules.free_module.FreeModule_generic_pid.submodule`,
        and
        :meth:`~sage.modules.free_module.FreeModule_generic_pid.submodule_with_basis`,

        EXAMPLES::

            sage: N = ToricLattice(3)
            sage: Ns = N.submodule([N.gen(0)])
            sage: Ns.span([N.gen(1)])
            Sublattice <N(0, 1, 0)>
            sage: Ns.submodule([N.gen(1)])
            Traceback (most recent call last):
            ...
            ArithmeticError: argument gens (= [N(0, 1, 0)]) does not generate a submodule of self
        """
        A = self.ambient_module()
        if base_ring is ZZ and all(g in A for g in gens):
            return ToricLattice_sublattice(A, gens)
        for g in gens:
            if is_ToricLatticeElement(g) and g not in A:
                raise ValueError("%s cannot generate a sublattice of %s"
                                 % (g, A))
        return super().span(gens, base_ring, *args, **kwds)

    def span_of_basis(self, basis, base_ring=ZZ, *args, **kwds):
        r"""
        Return the submodule with the given ``basis``.

        INPUT:

        - ``basis`` -- list of elements of the ambient vector space of
          ``self``.

        - ``base_ring`` -- (default: `\ZZ`) base ring for the generated module.

        OUTPUT:

        - submodule spanned by ``basis``.

        .. NOTE::

            The output need not be a submodule of ``self``, nor even of the
            ambient space. It must, however, be contained in the ambient
            vector space.

        See also :meth:`span`,
        :meth:`~sage.modules.free_module.FreeModule_generic_pid.submodule`,
        and
        :meth:`~sage.modules.free_module.FreeModule_generic_pid.submodule_with_basis`,

        EXAMPLES::

            sage: N = ToricLattice(3)
            sage: Ns = N.span_of_basis([(1,2,3)])
            sage: Ns.span_of_basis([(2,4,0)])
            Sublattice <N(2, 4, 0)>
            sage: Ns.span_of_basis([(1/5,2/5,0), (1/7,1/7,0)])
            Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [1/5 2/5   0]
            [1/7 1/7   0]

        Of course the input basis vectors must be linearly independent::

            sage: Ns.span_of_basis([(1,2,0), (2,4,0)])
            Traceback (most recent call last):
            ...
            ValueError: The given basis vectors must be linearly independent.
        """
        A = self.ambient_module()
        if base_ring is ZZ and all(g in A for g in basis):
            return ToricLattice_sublattice_with_basis(A, basis)
        for g in basis:
            if is_ToricLatticeElement(g) and g not in A:
                raise ValueError("%s cannot generate a sublattice of %s"
                                 % (g, A))
        return super().span_of_basis(basis, base_ring, *args, **kwds)


@richcmp_method
class ToricLattice_ambient(ToricLattice_generic, FreeModule_ambient_pid):
    r"""
    Create a toric lattice.

    See :class:`ToricLattice <ToricLatticeFactory>` for documentation.

    .. WARNING::

        There should be only one toric lattice with the given rank and
        associated names. Using this class directly to create toric lattices
        may lead to unexpected results. Please, use :class:`ToricLattice
        <ToricLatticeFactory>` to create toric lattices.

    TESTS::

        sage: N = ToricLattice(3, "N", "M", "N", "M")
        sage: N
        3-d lattice N
        sage: TestSuite(N).run()
    """

    Element = ToricLatticeElement

    def __init__(self, rank, name, dual_name, latex_name, latex_dual_name):
        r"""
        See :class:`ToricLattice <ToricLatticeFactory>` for documentation.

        TESTS::

            sage: ToricLattice(3, "N", "M", "N", "M")
            3-d lattice N
        """
        super().__init__(ZZ, rank)
        self._name = name
        self._dual_name = dual_name
        self._latex_name = latex_name
        self._latex_dual_name = latex_dual_name

    def _sage_input_(self, sib, coerced):
        r"""
        Return Sage command to reconstruct ``self``.

        See :mod:`sage.misc.sage_input` for details.

        EXAMPLES::

            sage: N = ToricLattice(3, "N", "M", "N", "M")
            sage: sage_input(N, verify=True)
            # Verified
            ToricLattice(3, 'N', 'M')

            sage: N = ToricLattice(3, "N")
            sage: sage_input(N, verify=True)
            # Verified
            ToricLattice(3, 'N', 'N*', 'N', 'N^*')
        """
        args = [self.rank(), self._name, self._dual_name]
        if self._latex_name != self._name or self._latex_dual_name != self._dual_name:
            args.extend([self._latex_name, self._latex_dual_name])
        return sib.name('ToricLattice')(*args)

    def __richcmp__(self, right, op):
        r"""
        Compare ``self`` and ``right``.

        INPUT:

        - ``right`` -- anything.

        OUTPUT:

        boolean

        There is equality if ``right`` is a toric lattice of the same
        dimension as ``self`` and their associated names are the
        same.

        TESTS::

            sage: N3 = ToricLattice(3)
            sage: N4 = ToricLattice(4)
            sage: M3 = N3.dual()
            sage: N3 < N4
            True
            sage: N3 > M3
            True
            sage: N3 == ToricLattice(3)
            True
        """
        if self is right:
            return rich_to_bool(op, 0)
        if type(self) is not type(right):
            return NotImplemented

        lx = self.rank()
        rx = right.rank()
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)
        # If lattices are the same as ZZ-modules, compare associated names
        return richcmp([self._name, self._dual_name,
                        self._latex_name, self._latex_dual_name],
                       [right._name, right._dual_name,
                        right._latex_name, right._latex_dual_name], op)

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT:

        - string.

        TESTS::

            sage: L = ToricLattice(3, "L")
            sage: L.dual()._latex_()
            'L^*'
        """
        return self._latex_name

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        - string.

        TESTS::

            sage: L = ToricLattice(3, "L")
            sage: L.dual()._repr_()
            '3-d lattice L*'
        """
        return "%d-d lattice %s" % (self.dimension(), self._name)

    def ambient_module(self):
        r"""
        Return the ambient module of ``self``.

        OUTPUT:

        - :class:`toric lattice <ToricLatticeFactory>`.

        .. NOTE::

            For any ambient toric lattice its ambient module is the lattice
            itself.

        EXAMPLES::

            sage: N = ToricLattice(3)
            sage: N.ambient_module()
            3-d lattice N
            sage: N.ambient_module() is N
            True
        """
        return self

    def dual(self):
        r"""
        Return the lattice dual to ``self``.

        OUTPUT:

        - :class:`toric lattice <ToricLatticeFactory>`.

        EXAMPLES::

            sage: N = ToricLattice(3)
            sage: N
            3-d lattice N
            sage: M = N.dual()
            sage: M
            3-d lattice M
            sage: M.dual() is N
            True

        Elements of dual lattices can act on each other::

            sage: n = N(1,2,3)
            sage: m = M(4,5,6)
            sage: n * m
            32
            sage: m * n
            32
        """
        if "_dual" not in self.__dict__:
            self._dual = ToricLattice(self.rank(), self._dual_name,
                          self._name, self._latex_dual_name, self._latex_name)
        return self._dual

    def plot(self, **options):
        r"""
        Plot ``self``.

        INPUT:

        - any options for toric plots (see :func:`toric_plotter.options
          <sage.geometry.toric_plotter.options>`), none are mandatory.

        OUTPUT:

        - a plot.

        EXAMPLES::

            sage: N = ToricLattice(3)
            sage: N.plot()                                                              # needs sage.plot
            Graphics3d Object
        """
        if "show_lattice" not in options:
            # Unless user made an explicit decision, we assume that lattice
            # should be visible no matter what is the size of the bounding box.
            options["show_lattice"] = True
        tp = ToricPlotter(options, self.degree())
        tp.adjust_options()
        return tp.plot_lattice()


class ToricLattice_sublattice_with_basis(ToricLattice_generic,
                                         FreeModule_submodule_with_basis_pid):
    r"""
    Construct the sublattice of ``ambient`` toric lattice with given ``basis``.

    INPUT (same as for
    :class:`~sage.modules.free_module.FreeModule_submodule_with_basis_pid`):

    - ``ambient`` -- ambient :class:`toric lattice <ToricLatticeFactory>` for
      this sublattice;

    - ``basis`` -- list of linearly independent elements of ``ambient``, these
      elements will be used as the default basis of the constructed
      sublattice;

    - see the base class for other available options.

    OUTPUT:

    - sublattice of a toric lattice with a user-specified basis.

    See also :class:`ToricLattice_sublattice` if you do not want to specify an
    explicit basis.

    EXAMPLES:

    The intended way to get objects of this class is to use
    :meth:`submodule_with_basis` method of toric lattices::

        sage: N = ToricLattice(3)
        sage: sublattice = N.submodule_with_basis([(1,1,0), (3,2,1)])
        sage: sublattice.has_user_basis()
        True
        sage: sublattice.basis()
        [
        N(1, 1, 0),
        N(3, 2, 1)
        ]

    Even if you have provided your own basis, you still can access the
    "standard" one::

        sage: sublattice.echelonized_basis()
        [
        N(1, 0, 1),
        N(0, 1, -1)
        ]
    """

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        - string.

        TESTS::

            sage: L = ToricLattice(3, "L")
            sage: L.submodule_with_basis([(3,2,1),(1,2,3)])
            Sublattice <L(3, 2, 1), L(1, 2, 3)>
            sage: print(L.submodule([(3,2,1),(1,2,3)])._repr_())
            Sublattice <L(1, 2, 3), L(0, 4, 8)>
        """
        s = 'Sublattice '
        s += '<'
        s += ', '.join(map(str,self.basis()))
        s += '>'
        return s

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT:

        - string.

        TESTS::

            sage: L = ToricLattice(3, "L")
            sage: L.submodule_with_basis([(3,2,1),(1,2,3)])._latex_()
            '\\left\\langle\\left(3,\\,2,\\,1\\right)_{L},
             \\left(1,\\,2,\\,3\\right)_{L}\\right\\rangle'
            sage: L.submodule([(3,2,1),(1,2,3)])._latex_()
            '\\left\\langle\\left(1,\\,2,\\,3\\right)_{L},
             \\left(0,\\,4,\\,8\\right)_{L}\\right\\rangle'
        """
        s = '\\left\\langle'
        s += ', '.join([ b._latex_() for b in self.basis() ])
        s += '\\right\\rangle'
        return s

    def dual(self):
        r"""
        Return the lattice dual to ``self``.

        OUTPUT:

        - a :class:`toric lattice quotient <ToricLattice_quotient>`.

        EXAMPLES::

            sage: N = ToricLattice(3)
            sage: Ns = N.submodule([(1,1,0), (3,2,1)])
            sage: Ns.dual()
            2-d lattice, quotient of 3-d lattice M by Sublattice <M(1, -1, -1)>
        """
        if "_dual" not in self.__dict__:
            if self is not self.saturation():
                raise ValueError("only dual lattices of saturated sublattices "
                                 "can be constructed! Got %s." % self)
            self._dual = (self.ambient_module().dual() /
                          self.basis_matrix().transpose().integer_kernel())
            self._dual._dual = self
        return self._dual

    def plot(self, **options):
        r"""
        Plot ``self``.

        INPUT:

        - any options for toric plots (see :func:`toric_plotter.options
          <sage.geometry.toric_plotter.options>`), none are mandatory.

        OUTPUT:

        - a plot.

        EXAMPLES::

            sage: N = ToricLattice(3)
            sage: sublattice = N.submodule_with_basis([(1,1,0), (3,2,1)])
            sage: sublattice.plot()                                                     # needs sage.plot
            Graphics3d Object

        Now we plot both the ambient lattice and its sublattice::

            sage: N.plot() + sublattice.plot(point_color="red")                         # needs sage.plot
            Graphics3d Object
        """
        if "show_lattice" not in options:
            # Unless user made an explicit decision, we assume that lattice
            # should be visible no matter what is the size of the bounding box.
            options["show_lattice"] = True
        if "lattice_filter" in options:
            old = options["lattice_filter"]
            options["lattice_filter"] = lambda pt: pt in self and old(pt)
        else:
            options["lattice_filter"] = lambda pt: pt in self
        tp = ToricPlotter(options, self.degree())
        tp.adjust_options()
        return tp.plot_lattice()


class ToricLattice_sublattice(ToricLattice_sublattice_with_basis,
                              FreeModule_submodule_pid):
    r"""
    Construct the sublattice of ``ambient`` toric lattice generated by ``gens``.

    INPUT (same as for
    :class:`~sage.modules.free_module.FreeModule_submodule_pid`):

    - ``ambient`` -- ambient :class:`toric lattice <ToricLatticeFactory>` for
      this sublattice;

    - ``gens`` -- list of elements of ``ambient`` generating the constructed
      sublattice;

    - see the base class for other available options.

    OUTPUT:

    - sublattice of a toric lattice with an automatically chosen basis.

    See also :class:`ToricLattice_sublattice_with_basis` if you want to
    specify an explicit basis.

    EXAMPLES:

    The intended way to get objects of this class is to use
    :meth:`submodule` method of toric lattices::

        sage: N = ToricLattice(3)
        sage: sublattice = N.submodule([(1,1,0), (3,2,1)])
        sage: sublattice.has_user_basis()
        False
        sage: sublattice.basis()
        [
        N(1, 0, 1),
        N(0, 1, -1)
        ]

    For sublattices without user-specified basis, the basis obtained above is
    the same as the "standard" one::

        sage: sublattice.echelonized_basis()
        [
        N(1, 0, 1),
        N(0, 1, -1)
        ]
    """
    pass


class ToricLattice_quotient_element(FGP_Element):
    r"""
    Create an element of a toric lattice quotient.

    .. WARNING::

        You probably should not construct such elements explicitly.

    INPUT:

    - same as for :class:`~sage.modules.fg_pid.fgp_element.FGP_Element`.

    OUTPUT:

    - element of a toric lattice quotient.

    TESTS::

        sage: N = ToricLattice(3)
        sage: sublattice = N.submodule([(1,1,0), (3,2,1)])
        sage: Q = N/sublattice
        sage: e = Q(1,2,3)
        sage: e
        N[1, 2, 3]
        sage: e2 = Q(N(2,3,3))
        sage: e2
        N[2, 3, 3]
        sage: e == e2
        True
        sage: e.vector()
        (-4)
        sage: e2.vector()
        (-4)
    """

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT:

        - string.

        TESTS::

            sage: N = ToricLattice(3)
            sage: Ns = N.submodule([N(2,4,0), N(9,12,0)])
            sage: Q = N/Ns
            sage: print(Q.gen(0)._latex_())
            \left[0,\,1,\,0\right]_{N}
        """
        return latex(self.lift()).replace("(", "[", 1).replace(")", "]", 1)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        - string.

        TESTS::

            sage: N = ToricLattice(3)
            sage: Ns = N.submodule([N(2,4,0), N(9,12,0)])
            sage: Q = N/Ns
            sage: print(Q.gen(0)._repr_())
            N[0, 1, 0]
        """
        return str(self.lift()).replace("(", "[", 1).replace(")", "]", 1)

    def set_immutable(self):
        r"""
        Make ``self`` immutable.

        OUTPUT:

        - none.

        .. note:: Elements of toric lattice quotients are always immutable, so
            this method does nothing, it is introduced for compatibility
            purposes only.

        EXAMPLES::

            sage: N = ToricLattice(3)
            sage: Ns = N.submodule([N(2,4,0), N(9,12,0)])
            sage: Q = N/Ns
            sage: Q.0.set_immutable()
        """
        pass


class ToricLattice_quotient(FGP_Module_class):
    r"""
    Construct the quotient of a toric lattice ``V`` by its sublattice ``W``.

    INPUT:

    - ``V`` -- ambient toric lattice;

    - ``W`` -- sublattice of ``V``;

    - ``check`` -- (default: ``True``) whether to check correctness of input
      or not.

    If the quotient is one-dimensional and torsion free, the following
    two mutually exclusive keyword arguments are also allowed. They
    decide the sign choice for the (single) generator of the quotient
    lattice:

    - ``positive_point`` -- a lattice point of ``self`` not in the
      sublattice ``sub`` (that is, not zero in the quotient
      lattice). The quotient generator will be in the same direction
      as ``positive_point``.

    - ``positive_dual_point`` -- a dual lattice point. The quotient
      generator will be chosen such that its lift has a positive
      product with ``positive_dual_point``. Note: if
      ``positive_dual_point`` is not zero on the sublattice ``sub``,
      then the notion of positivity will depend on the choice of lift!

    Further given named arguments are passed to the constructor of an FGP
    module.

    OUTPUT:

    - quotient of ``V`` by ``W``.

    EXAMPLES:

    The intended way to get objects of this class is to use
    :meth:`quotient` method of toric lattices::

        sage: N = ToricLattice(3)
        sage: sublattice = N.submodule([(1,1,0), (3,2,1)])
        sage: Q = N/sublattice
        sage: Q
        1-d lattice, quotient of 3-d lattice N by Sublattice <N(1, 0, 1), N(0, 1, -1)>
        sage: Q.gens()
        (N[1, 0, 0],)

    Here, ``sublattice`` happens to be of codimension one in ``N``. If
    you want to prescribe the sign of the quotient generator, you can
    do either::

        sage: Q = N.quotient(sublattice, positive_point=N(0,0,-1)); Q
        1-d lattice, quotient of 3-d lattice N by Sublattice <N(1, 0, 1), N(0, 1, -1)>
        sage: Q.gens()
        (N[1, 0, 0],)

    or::

        sage: M = N.dual()
        sage: Q = N.quotient(sublattice, positive_dual_point=M(1,0,0)); Q
        1-d lattice, quotient of 3-d lattice N by Sublattice <N(1, 0, 1), N(0, 1, -1)>
        sage: Q.gens()
        (N[1, 0, 0],)

    TESTS::

        sage: loads(dumps(Q)) == Q
        True
        sage: loads(dumps(Q)).gens() == Q.gens()
        True
    """

    def __init__(self, V, W, check=True, positive_point=None, positive_dual_point=None, **kwds):
        r"""
        The constructor

        See :class:`ToricLattice_quotient` for an explanation of the arguments.

        EXAMPLES::

            sage: N = ToricLattice(3)
            sage: from sage.geometry.toric_lattice import ToricLattice_quotient
            sage: ToricLattice_quotient(N, N.span([N(1,2,3)]))
            2-d lattice, quotient of 3-d lattice N by Sublattice <N(1, 2, 3)>

        An ``ArithmeticError`` will be raised if ``W`` is not a
        sublattice of ``V``::

            sage: N = ToricLattice(3)
            sage: Ns = N.submodule([N.gen(0)])
            sage: Ns
            Sublattice <N(1, 0, 0)>
            sage: Ns.span([N.gen(1)])
            Sublattice <N(0, 1, 0)>
            sage: Ns.quotient(Ns.span([N.gen(1)]))
            Traceback (most recent call last):
            ...
            ArithmeticError: W must be a sublattice of V
        """
        if check:
            try:
                W = V.submodule(W)
            except (TypeError, ArithmeticError):
                raise ArithmeticError("W must be a sublattice of V")
        super().__init__(V, W, check, **kwds)
        if (positive_point, positive_dual_point) == (None, None):
            self._flip_sign_of_generator = False
            return

        self._flip_sign_of_generator = False
        assert self.is_torsion_free() and self.ngens()==1, \
            'You may only specify a positive direction in the codimension one case.'
        quotient_generator = self.gen(0)
        lattice = self.V().ambient_module()
        if (positive_point is not None) and (positive_dual_point is None):
            assert positive_point in lattice, 'positive_point must be a lattice point.'
            point_quotient = self(positive_point)
            scalar_product = quotient_generator.vector()[0] * point_quotient.vector()[0]
            if scalar_product==0:
                raise ValueError(str(positive_point)+' is zero in the quotient.')
        elif (positive_point is None) and (positive_dual_point is not None):
            assert positive_dual_point in lattice.dual(), 'positive_dual_point must be a dual lattice point.'
            scalar_product = quotient_generator.lift() * positive_dual_point
            if scalar_product==0:
                raise ValueError(str(positive_dual_point)+' is zero on the lift of the quotient generator.')
        else:
            raise ValueError('You may not specify both positive_point and positive_dual_point.')
        self._flip_sign_of_generator = (scalar_product<0)

    def gens(self):
        """
        Return the generators of the quotient.

        OUTPUT:

        A tuple of :class:`ToricLattice_quotient_element` generating
        the quotient.

        EXAMPLES::

            sage: N = ToricLattice(3)
            sage: Q = N.quotient(N.span([N(1,2,3), N(0,2,1)]), positive_point=N(0,-1,0))
            sage: Q.gens()
            (N[0, -1, 0],)
        """
        gens = self.smith_form_gens()
        if self._flip_sign_of_generator:
            assert len(gens)==1
            return (-gens[0],)
        else:
            return gens

    # Should be overridden in derived classes.
    Element = ToricLattice_quotient_element

    def _element_constructor_(self, *x, **kwds):
        r"""
        Construct an element of ``self``.

        INPUT:

        - element of a compatible toric object (lattice, sublattice, quotient)
          or something that defines such an element (list, generic vector,
          etc.).

        OUTPUT:

        - :class:`toric lattice quotient element
          <ToricLattice_quotient_element>`.

        EXAMPLES::

            sage: N = ToricLattice(3)
            sage: Ns = N.submodule([N(2,4,0), N(9,12,0)])
            sage: Q = N/Ns
            sage: x = Q(1,2,3)  # indirect doctest
            sage: x
            N[1, 2, 3]
            sage: type(x)
            <class 'sage.geometry.toric_lattice.ToricLattice_quotient_with_category.element_class'>
            sage: x is Q(x)
            True
            sage: x.parent() is Q
            True
            sage: x == Q(N(1,2,3))
            True
            sage: y = Q(3,6,3)
            sage: y
            N[3, 6, 3]
            sage: x == y
            True
        """
        if len(x) == 1 and (x[0] not in ZZ or x[0] == 0):
            x = x[0]
        if parent(x) is self:
            return x
        try:
            x = x.lift()
        except AttributeError:
            pass
        try:
            return self.element_class(self, self._V(x), **kwds)
        except TypeError:
            return self.linear_combination_of_smith_form_gens(x)

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT:

        - string.

        TESTS::

            sage: N = ToricLattice(3)
            sage: Ns = N.submodule([N(2,4,0), N(9,12,0)])
            sage: Q = N/Ns
            sage: print(Q._latex_())
            N / \left\langle\left(1,\,8,\,0\right)_{N}, \left(0,\,12,\,0\right)_{N}\right\rangle
            sage: Ns = N.submodule([N(1,4,0)])
            sage: Q = N/Ns
            sage: print(Q._latex_())
            N / \left\langle\left(1,\,4,\,0\right)_{N}\right\rangle
        """
        return "%s / %s" % (latex(self.V()), latex(self.W()))

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        - string.

        TESTS::

            sage: N = ToricLattice(3)
            sage: Ns = N.submodule([N(2,4,0), N(9,12,0)])
            sage: Q = N/Ns
            sage: print(Q._repr_())
            Quotient with torsion of 3-d lattice N
            by Sublattice <N(1, 8, 0), N(0, 12, 0)>
            sage: Ns = N.submodule([N(1,4,0)])
            sage: Q = N/Ns
            sage: print(Q._repr_())
            2-d lattice, quotient of 3-d lattice N
            by Sublattice <N(1, 4, 0)>
        """
        if self.is_torsion_free():
            return "%d-d lattice, quotient of %s by %s" % (self.rank(),
                                                           self.V(), self.W())
        else:
            return "Quotient with torsion of %s by %s" % (self.V(), self.W())

    def _module_constructor(self, V, W, check=True):
        r"""
        Construct new quotient modules.

        INPUT:

        - ``V`` -- ambient toric lattice;

        - ``W`` -- sublattice of ``V``;

        - ``check`` -- (default: ``True``) whether to check
          correctness of input or not.

        TESTS::

            sage: N = ToricLattice(3)
            sage: Ns = N.submodule([N(2,4,0), N(9,12,0)])
            sage: Q = N/Ns; Q
            Quotient with torsion of 3-d lattice N by Sublattice <N(1, 8, 0), N(0, 12, 0)>
            sage: Q._module_constructor(N,Ns)
            Quotient with torsion of 3-d lattice N by Sublattice <N(1, 8, 0), N(0, 12, 0)>
        """
        return ToricLattice_quotient(V,W,check)

    def base_extend(self, R):
        r"""
        Return the base change of ``self`` to the ring ``R``.

        INPUT:

        - ``R`` -- either `\ZZ` or `\QQ`.

        OUTPUT:

        - ``self`` if `R=\ZZ`, quotient of the base extension of the ambient
          lattice by the base extension of the sublattice if `R=\QQ`.

        EXAMPLES::

            sage: N = ToricLattice(3)
            sage: Ns = N.submodule([N(2,4,0), N(9,12,0)])
            sage: Q = N/Ns
            sage: Q.base_extend(ZZ) is Q
            True
            sage: Q.base_extend(QQ)
            Vector space quotient V/W of dimension 1 over Rational Field where
            V: Vector space of dimension 3 over Rational Field
            W: Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [1 0 0]
            [0 1 0]
        """
        if R is ZZ:
            return self
        if R is QQ:
            return self.V().base_extend(R) / self.W().base_extend(R)
        raise NotImplementedError("quotients of toric lattices can only be "
                                  "extended to ZZ or QQ, not %s!" % R)

    def is_torsion_free(self):
        r"""
        Check if ``self`` is torsion-free.

        OUTPUT:

        - ``True`` is ``self`` has no torsion and ``False`` otherwise.

        EXAMPLES::

            sage: N = ToricLattice(3)
            sage: Ns = N.submodule([N(2,4,0), N(9,12,0)])
            sage: Q = N/Ns
            sage: Q.is_torsion_free()
            False
            sage: Ns = N.submodule([N(1,4,0)])
            sage: Q = N/Ns
            sage: Q.is_torsion_free()
            True
        """
        return sum(self.invariants()) == 0

    def dual(self):
        r"""
        Return the lattice dual to ``self``.

        OUTPUT:

        - a :class:`toric lattice quotient <ToricLattice_quotient>`.

        EXAMPLES::

            sage: N = ToricLattice(3)
            sage: Ns = N.submodule([(1, -1, -1)])
            sage: Q = N / Ns
            sage: Q.dual()
            Sublattice <M(1, 0, 1), M(0, 1, -1)>
        """
        if "_dual" not in self.__dict__:
            self._dual = self.V().dual().submodule(
                    self.W().basis_matrix().transpose().integer_kernel().gens())
            self._dual._dual = self
        return self._dual

    def rank(self):
        r"""
        Return the rank of ``self``.

        OUTPUT:

        Integer. The dimension of the free part of the quotient.

        EXAMPLES::

            sage: N = ToricLattice(3)
            sage: Ns = N.submodule([N(2,4,0), N(9,12,0)])
            sage: Q = N/Ns
            sage: Q.ngens()
            2
            sage: Q.rank()
            1
            sage: Ns = N.submodule([N(1,4,0)])
            sage: Q = N/Ns
            sage: Q.ngens()
            2
            sage: Q.rank()
            2
        """
        return self.V().rank() - self.W().rank()

    dimension = rank

    def coordinate_vector(self, x, reduce=False):
        """
        Return coordinates of ``x`` with respect to the optimized
        representation of ``self``.

        INPUT:

        - ``x`` -- element of ``self`` or convertible to ``self``

        - ``reduce`` -- (default: ``False``); if ``True``, reduce coefficients
          modulo invariants

        OUTPUT:

        The coordinates as a vector.

        EXAMPLES::

            sage: N = ToricLattice(3)
            sage: Q = N.quotient(N.span([N(1,2,3), N(0,2,1)]), positive_point=N(0,-1,0))
            sage: q = Q.gen(0); q
            N[0, -1, 0]
            sage: q.vector()  # indirect test
            (1)
            sage: Q.coordinate_vector(q)
            (1)
        """
        coordinates = super().coordinate_vector(x, reduce)
        if self._flip_sign_of_generator:
            assert len(coordinates) == 1, "Sign flipped for a multi-dimensional quotient!"
            return -coordinates
        else:
            return coordinates
