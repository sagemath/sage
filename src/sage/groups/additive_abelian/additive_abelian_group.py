r"""
Additive Abelian Groups

Additive abelian groups are just modules over `\ZZ`. Hence the classes in this
module derive from those in the module :mod:`sage.modules.fg_pid`. The only
major differences are in the way elements are printed.
"""

from sage.groups.old import AbelianGroup
from sage.modules.fg_pid.fgp_module import FGP_Module_class
from sage.modules.fg_pid.fgp_element import FGP_Element
from sage.rings.integer_ring import ZZ

def AdditiveAbelianGroup(invs, remember_generators = True):
    r"""
    Construct a finitely-generated additive abelian group.

    INPUT:

    - ``invs`` (list of integers): the invariants.
      These should all be greater than or equal to zero.

    - ``remember_generators`` (boolean): whether or not to fix a set of
      generators (corresponding to the given invariants, which need not be in
      Smith form).

    OUTPUT:

    The abelian group `\bigoplus_i \ZZ / n_i \ZZ`, where `n_i` are the invariants.

    EXAMPLES::

        sage: AdditiveAbelianGroup([0, 2, 4])
        Additive abelian group isomorphic to Z + Z/2 + Z/4

    An example of the ``remember_generators`` switch::

        sage: G = AdditiveAbelianGroup([0, 2, 3]); G
        Additive abelian group isomorphic to Z + Z/2 + Z/3
        sage: G.gens()
        ((1, 0, 0), (0, 1, 0), (0, 0, 1))

        sage: H = AdditiveAbelianGroup([0, 2, 3], remember_generators = False); H
        Additive abelian group isomorphic to Z/6 + Z
        sage: H.gens()
        ((0, 1, 1), (1, 0, 0))

    There are several ways to create elements of an additive abelian group.
    Realize that there are two sets of generators:  the "obvious" ones composed
    of zeros and ones, one for each invariant given to construct the group, the
    other being a set of minimal generators.  Which set is the default varies
    with the use of the ``remember_generators`` switch.

    First with "obvious" generators.  Note that a raw list will use the
    minimal generators and a vector (a module element) will use the generators
    that pair up naturally with the invariants.  We create the same element
    repeatedly. ::

        sage: H=AdditiveAbelianGroup([3,2,0], remember_generators=True)
        sage: H.gens()
        ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        sage: [H.0, H.1, H.2]
        [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        sage: p=H.0+H.1+6*H.2; p
        (1, 1, 6)

        sage: H.smith_form_gens()
        ((2, 1, 0), (0, 0, 1))
        sage: q=H.linear_combination_of_smith_form_gens([5,6]); q
        (1, 1, 6)
        sage: p==q
        True

        sage: r=H(vector([1,1,6])); r
        (1, 1, 6)
        sage: p==r
        True

        sage: s=H(p)
        sage: p==s
        True

    Again, but now where the generators are the minimal set.  Coercing a
    list or a vector works as before, but the default generators are different. ::

        sage: G=AdditiveAbelianGroup([3,2,0], remember_generators=False)
        sage: G.gens()
        ((2, 1, 0), (0, 0, 1))
        sage: [G.0, G.1]
        [(2, 1, 0), (0, 0, 1)]
        sage: p=5*G.0+6*G.1; p
        (1, 1, 6)

        sage: H.smith_form_gens()
        ((2, 1, 0), (0, 0, 1))
        sage: q=G.linear_combination_of_smith_form_gens([5,6]); q
        (1, 1, 6)
        sage: p==q
        True

        sage: r=G(vector([1,1,6])); r
        (1, 1, 6)
        sage: p==r
        True

        sage: s=H(p)
        sage: p==s
        True
    """
    invs = [ZZ(x) for x in invs]
    if not all(x >= 0 for x in invs):
        raise ValueError("Invariants must be nonnegative")
    A, B = cover_and_relations_from_invariants(invs)
    if remember_generators:
        G = AdditiveAbelianGroup_fixed_gens(A, B, A.gens())
    else:
        G = AdditiveAbelianGroup_class(A, B)
    return G


def cover_and_relations_from_invariants(invs):
    r"""
    A utility function to construct modules required to initialize the super class.

    Given a list of integers, this routine constructs the obvious pair of
    free modules such that the quotient of the two free modules over `\ZZ`
    is naturally isomorphic to the corresponding product of cyclic modules
    (and hence isomorphic to a direct sum of cyclic groups).

    EXAMPLES::

        sage: from sage.groups.additive_abelian.additive_abelian_group import cover_and_relations_from_invariants as cr
        sage: cr([0,2,3])
        (Ambient free module of rank 3 over the principal ideal domain Integer Ring, Free module of degree 3 and rank 2 over Integer Ring
        Echelon basis matrix:
        [0 2 0]
        [0 0 3])
    """
    n = len(invs)
    A = ZZ**n
    B = A.span([A.gen(i) * invs[i] for i in range(n)])
    return (A, B)



class AdditiveAbelianGroupElement(FGP_Element):
    """
    An element of an :class:`AdditiveAbelianGroup_class`.
    """

    def _hermite_lift(self):
        r"""
        This gives a certain canonical lifting of elements of this group
        (represented as a quotient `G/H` of free abelian groups) to `G`, using
        the Hermite normal form of the matrix of relations.

        Mainly used by the ``_repr_`` method.

        EXAMPLES::

            sage: A = AdditiveAbelianGroup([2, 3])
            sage: v = 3000001 * A.0
            sage: v.lift()
            (3000001, 0)
            sage: v._hermite_lift()
            (1, 0)
        """
        y = self.lift()
        H = self.parent().W().basis_matrix()
        pivot_rows = H.pivot_rows()
        pivots = H.pivots()

        for i in range(H.nrows()):
            if i in pivot_rows:
                j = pivots[i]
                N = H[i,j]
                a = (y[j] - (y[j] % N)) // N
                y = y - a*H.row(i)
        return y

    def _repr_(self):
        r"""
        String representation. This uses a canonical lifting of elements of
        this group (represented as a quotient `G/H` of free abelian groups) to
        `G`, using the Hermite normal form of the matrix of relations.

        EXAMPLES::

            sage: G = AdditiveAbelianGroup([2,3])
            sage: repr(G.gen(0)) # indirect doctest
            '(1, 0)'
            sage: a = 13*G.gen(0); repr(a) # indirect doctest
            '(1, 0)'
            sage: a._x
            (13, 0)
        """
        return repr(self._hermite_lift())


# Note: It's important that the class inherits from FGP_Module_class first,
# since we want to inherit things like __hash__ from there rather than the
# hyper-generic implementation for abstract abelian groups.

class AdditiveAbelianGroup_class(FGP_Module_class, AbelianGroup):
    r"""
    An additive abelian group, implemented using the `\ZZ`-module machinery.

    INPUT:

    - ``cover`` -- the covering group as `\ZZ`-module.

    - ``relations`` -- the relations as submodule of ``cover``.
    """

    # The element class must be overridden in derived classes
    Element = AdditiveAbelianGroupElement

    def __init__(self, cover, relations):
        r"""
        EXAMPLES::

            sage: G = AdditiveAbelianGroup([0]); G # indirect doctest
            Additive abelian group isomorphic to Z
            sage: G == loads(dumps(G))
            True
        """
        FGP_Module_class.__init__(self, cover, relations)

    def _repr_(self):
        r"""
        String representation of this group.

        EXAMPLES::

            sage: AdditiveAbelianGroup([0, 2, 3])._repr_()
            'Additive abelian group isomorphic to Z + Z/2 + Z/3'
        """
        if self.V().rank() == 0:
            return "Trivial group"
        else:
            return "Additive abelian group isomorphic to %s" % self.short_name()

    def _latex_(self):
        r"""
        Returns a Latex representation of the group, using the invariants.

        EXAMPLES::

            sage: G=AdditiveAbelianGroup([66, 77, 0, 0])
            sage: G._latex_()
            '\\frac{\\ZZ}{11\\ZZ} \\oplus \\frac{\\ZZ}{462\\ZZ} \\oplus \\ZZ \\oplus \\ZZ'

        A trivial group is represented as zero, rather than Z/1Z. ::

            sage: G=AdditiveAbelianGroup([1])
            sage: G._latex_()
            '0'
        """
        inv = self.invariants()
        if not inv:
            inv = (1,)
        terms=[]
        for i in range(len(inv)):
            if inv[i] == 0:
                terms.append('\\ZZ')
            elif inv[i] == 1:
                terms.append('0')
            else:
                terms.append('\\frac{\\ZZ}{' + str(inv[i]) + '\\ZZ}')
        return ' \\oplus '.join(terms)

    def short_name(self):
        r"""
        Return a name for the isomorphism class of this group.

        EXAMPLES::

            sage: AdditiveAbelianGroup([0, 2,4]).short_name()
            'Z + Z/2 + Z/4'
            sage: AdditiveAbelianGroup([0, 2, 3]).short_name()
            'Z + Z/2 + Z/3'
        """
        from sage.rings.infinity import Infinity as oo
        invs = [j.additive_order() for j in self.gens()]
        if not invs:
            return "Trivial group"
        return " + ".join("Z" if j == +oo else "Z/%s"%j for j in invs)

    def _module_constructor(self, cover, relations, check=True):
        r"""
        Construct quotients of groups.

        INPUT:

        - ``cover`` -- the covering group as `\ZZ`-module.

        - ``relations`` -- the relations as submodule of ``cover``.

        - ``check`` -- ignored, present for compatibility with ``fg_pid`` code.

        EXAMPLES::

            sage: G = AdditiveAbelianGroup([0, 4, 2]); G
            Additive abelian group isomorphic to Z + Z/4 + Z/2
            sage: H = G.submodule([G.1]); H
            Additive abelian group isomorphic to Z/4
            sage: G/H    # indirect test
            Additive abelian group isomorphic to Z/2 + Z
            sage: G._module_constructor(G.cover(),H.cover()+G.relations())
            Additive abelian group isomorphic to Z/2 + Z

        TESTS:

        Check that :trac:`21027` is fixed::

            sage: G = AdditiveAbelianGroup([2,2,2])
            sage: phi = G.hom([G.0, G.0, G.0])
            sage: phi.image()
            Additive abelian group isomorphic to Z/2
        """
        return AdditiveAbelianGroup_class(cover, relations)

    def order(self):
        r"""
        Return the order of this group (an integer or infinity)

        EXAMPLES::

            sage: AdditiveAbelianGroup([2,4]).order()
            8
            sage: AdditiveAbelianGroup([0, 2,4]).order()
            +Infinity
            sage: AdditiveAbelianGroup([]).order()
            1
        """
        return self.cardinality()

    def exponent(self):
        r"""
        Return the exponent of this group (the smallest positive integer `N`
        such that `Nx = 0` for all `x` in the group). If there is no such
        integer, return 0.

        EXAMPLES::

            sage: AdditiveAbelianGroup([2,4]).exponent()
            4
            sage: AdditiveAbelianGroup([0, 2,4]).exponent()
            0
            sage: AdditiveAbelianGroup([]).exponent()
            1
        """
        if not self.invariants():
            return 1
        else:
            ann =  self.annihilator().gen()
            if ann:
                return ann
            return ZZ(0)

    def is_multiplicative(self):
        r"""
        Return False since this is an additive group.

        EXAMPLES::

            sage: AdditiveAbelianGroup([0]).is_multiplicative()
            False
        """
        return False

    def is_cyclic(self):
        r"""
        Returns ``True`` if the group is cyclic.

        EXAMPLES:

        With no common factors between the orders of the generators,
        the group will be cyclic. ::

            sage: G=AdditiveAbelianGroup([6, 7, 55])
            sage: G.is_cyclic()
            True

        Repeating primes in the orders will create a non-cyclic group. ::

            sage: G=AdditiveAbelianGroup([6, 15, 21, 33])
            sage: G.is_cyclic()
            False

        A trivial group is trivially cyclic. ::

            sage: T=AdditiveAbelianGroup([1])
            sage: T.is_cyclic()
            True
        """
        # One invariant is characteristic of a cyclic group
        # while zero invariants is characteristic of the trivial group
        return len(self.invariants()) < 2


class AdditiveAbelianGroup_fixed_gens(AdditiveAbelianGroup_class):
    r"""
    A variant which fixes a set of generators, which need not be in Smith form
    (or indeed independent).
    """
    def __init__(self, cover, rels, gens):
        r"""
        Standard initialisation function

        EXAMPLES::

            sage: AdditiveAbelianGroup([3]) # indirect doctest
            Additive abelian group isomorphic to Z/3
        """
        AdditiveAbelianGroup_class.__init__(self, cover, rels)
        self._orig_gens = tuple(self(x) for x in gens)

    def gens(self):
        r"""
        Return the specified generators for self (as a tuple). Compare
        ``self.smithform_gens()``.

        EXAMPLES::

            sage: G = AdditiveAbelianGroup([2,3])
            sage: G.gens()
            ((1, 0), (0, 1))
            sage: G.smith_form_gens()
            ((1, 2),)
        """
        return self._orig_gens

    def identity(self):
        r"""
        Return the identity (zero) element of this group.

        EXAMPLES::

            sage: G = AdditiveAbelianGroup([2, 3])
            sage: G.identity()
            (0, 0)
        """
        return self(0)

    def permutation_group(self):
        r"""
        Return the permutation group attached to this group.

        EXAMPLES::

            sage: G = AdditiveAbelianGroup([2, 3])
            sage: G.permutation_group()
            Permutation Group with generators [(3,4,5), (1,2)]

        TESTS:

        Check that :trac:`25692` is fixed::

            sage: G = AdditiveAbelianGroup([0])
            sage: G.permutation_group()
            Traceback (most recent call last):
            ...
            TypeError: Additive Abelian group must be finite
        """
        # GAP does not support infinite permutation groups
        if not self.is_finite():
            raise TypeError('Additive Abelian group must be finite')
        from sage.groups.perm_gps.permgroup import PermutationGroup
        s = 'Image(IsomorphismPermGroup(AbelianGroup(%s)))'%(list(self.invariants()),)
        return PermutationGroup(gap_group=s)

