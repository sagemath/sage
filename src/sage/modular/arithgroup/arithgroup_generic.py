r"""
Arithmetic subgroups (finite index subgroups of `{\rm SL}_2(\ZZ)`)
"""
################################################################################
#
#       Copyright (C) 2009, The Sage Group -- http://www.sagemath.org/
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
#
################################################################################

from sage.groups.old import Group
from sage.categories.groups import Groups
from sage.rings.integer_ring import ZZ
from sage.arith.all import lcm
from sage.misc.cachefunc import cached_method
from copy import copy # for making copies of lists of cusps
from sage.modular.modsym.p1list import lift_to_sl2z
from sage.modular.cusps import Cusp

from sage.misc.lazy_import import lazy_import
lazy_import('sage.modular.arithgroup.congroup_sl2z', 'SL2Z')
from sage.structure.element import parent

from .arithgroup_element import ArithmeticSubgroupElement, M2Z as Mat2Z


def is_ArithmeticSubgroup(x):
    r"""
    Return ``True`` if ``x`` is of type :class:`ArithmeticSubgroup`.

    EXAMPLES::

        sage: from sage.modular.arithgroup.all import is_ArithmeticSubgroup
        sage: is_ArithmeticSubgroup(GL(2, GF(7)))
        False
        sage: is_ArithmeticSubgroup(Gamma0(4))
        True
    """
    return isinstance(x, ArithmeticSubgroup)


class ArithmeticSubgroup(Group):
    r"""
    Base class for arithmetic subgroups of `{\rm SL}_2(\ZZ)`. Not
    intended to be used directly, but still includes quite a few
    general-purpose routines which compute data about an arithmetic subgroup
    assuming that it has a working element testing routine.
    """

    Element = ArithmeticSubgroupElement

    def __init__(self):
        r"""
        Standard init routine.

        EXAMPLES::

            sage: G = Gamma1(7)
            sage: G.category() # indirect doctest
            Category of infinite groups
        """
        Group.__init__(self, category=Groups().Infinite())

    def _repr_(self):
        r"""
        Return the string representation of self.

        NOTE: This function should be overridden by all subclasses.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup()._repr_()
            'Generic arithmetic subgroup of SL2Z'
        """
        return "Generic arithmetic subgroup of SL2Z"

    def _repr_option(self, key):
        """
        Metadata about the :meth:`_repr_` output.

        See :meth:`sage.structure.parent._repr_option` for details.

        EXAMPLES::

            sage: Gamma1(7)._repr_option('element_ascii_art')
            True
        """
        if key == 'element_ascii_art':
            return True
        return super(ArithmeticSubgroup, self)._repr_option(key)

    def __reduce__(self):
        r"""
        Used for pickling self.

        NOTE: This function should be overridden by all subclasses.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup().__reduce__()
            Traceback (most recent call last):
            ...
            NotImplementedError: all subclasses must define a __reduce__ method
        """
        raise NotImplementedError("all subclasses must define a __reduce__ method")

    def _element_constructor_(self, x, check=True):
        r"""
        Create an element of this congruence subgroup from x.

        If the optional flag check is True (default), check whether
        x actually gives an element of self.

        EXAMPLES::

            sage: G = Gamma(5)
            sage: G([1, 0, -10, 1]) # indirect doctest
            [ 1   0]
            [-10  1]
            sage: G(matrix(ZZ, 2, [26, 5, 5, 1]))
            [26  5]
            [ 5  1]
            sage: G([1, 1, 6, 7])
            Traceback (most recent call last):
            ...
            TypeError: matrix [1 1]
            [6 7] is not an element of Congruence Subgroup Gamma(5)
        """
        # Do not override this function! Derived classes should override
        # _contains_sl2.
        x = SL2Z(x, check)
        if not check or x in self:
            return x
        raise TypeError("matrix %s is not an element of %s" % (x, self))

    def __contains__(self, x):
        r"""
        Test if x is an element of this group. This checks that x defines (is?) a 2x2 integer matrix of determinant 1, and
        then hands over to the routine _contains_sl2, which derived classes should implement.

        EXAMPLES::

            sage: [1,2] in SL2Z # indirect doctest
            False
            sage: [1,2,0,1] in SL2Z # indirect doctest
            True
            sage: SL2Z([1,2,0,1]) in Gamma(3) # indirect doctest
            False
            sage: -1 in SL2Z
            True
            sage: 2 in SL2Z
            False
        """
        # Do not override this function! Derived classes should override
        # _contains_sl2.
        if isinstance(x, list) and len(x) == 4:
            if not (x[0] in ZZ and x[1] in ZZ and x[2] in ZZ and x[3] in ZZ):
                return False
            a,b,c,d = map(ZZ, x)
            if a*d - b*c != 1:
                return False
            return self._contains_sl2(a,b,c,d)
        else:
            if parent(x) is not SL2Z:
                try:
                    y = SL2Z(x)
                except TypeError:
                    return False
                x = y
            return self._contains_sl2(x.a(),x.b(),x.c(),x.d())

    def _contains_sl2(self, a,b,c,d):
        r"""
        Test whether the matrix [a,b;c,d], which may be assumed to have
        determinant 1, is an element of self. This must be overridden by all
        subclasses.

        EXAMPLES::

            sage: G = sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup()
            sage: 1 in G
            Traceback (most recent call last):
            ...
            NotImplementedError: Please implement _contains_sl2 for <class 'sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup_with_category'>
        """
        raise NotImplementedError("Please implement _contains_sl2 for %s" % self.__class__)

    def __hash__(self):
        r"""
        Return a hash of self.

        EXAMPLES::

            sage: h = hash(Gamma0(11))
            sage: h = hash(Gamma1(11))

        TESTS:

        We test that :trac:`18743` is fixed::

            sage: G1 = GammaH(37,[4]); G1
            Congruence Subgroup Gamma_H(37) with H generated by [4]
            sage: G2 = GammaH(37,[4,16]); G2
            Congruence Subgroup Gamma_H(37) with H generated by [4, 7]
            sage: G1 == G2
            True
            sage: hash(G1) == hash(G2)
            True
            sage: set([G1,G2])
            {Congruence Subgroup Gamma_H(37) with H generated by [4]}
        """
        return hash((self.level(), self.index()))

    def matrix_space(self):
        """
        Return the parent space of the matrices, which is always
        ``MatrixSpace(ZZ, 2)``.

        EXAMPLES::

            sage: Gamma(3).matrix_space()
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
        """
        return Mat2Z

    def is_parent_of(self, x):
        r"""
        Check whether this group is a valid parent for the element x. Required
        by Sage's testing framework.

        EXAMPLES::

            sage: Gamma(3).is_parent_of(ZZ(1))
            False
            sage: Gamma(3).is_parent_of([1,0,0,1])
            False
            sage: Gamma(3).is_parent_of(SL2Z([1,1,0,1]))
            False
            sage: Gamma(3).is_parent_of(SL2Z(1))
            True
        """
        return (parent(x) == SL2Z and x in self)

    def coset_reps(self, G=None):
        r"""
        Return right coset representatives for self \\ G, where G is another
        arithmetic subgroup that contains self.  If G = None, default to G =
        SL2Z.

        For generic arithmetic subgroups G this is carried out by Todd-Coxeter
        enumeration; here G is treated as a black box, implementing nothing but
        membership testing.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup().coset_reps()
            Traceback (most recent call last):
            ...
            NotImplementedError: Please implement _contains_sl2 for <class 'sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup_with_category'>
            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup.coset_reps(Gamma0(3))
            [
            [1 0]  [ 0 -1]  [ 0 -1]  [ 0 -1]
            [0 1], [ 1  0], [ 1  1], [ 1  2]
            ]
        """
        return self.todd_coxeter(G)[0]

    @cached_method
    def todd_coxeter(self, G=None, on_right=True):
        r"""
        Compute coset representatives for self \\ G and action of standard
        generators on them via Todd-Coxeter enumeration.

        If ``G`` is ``None``, default to ``SL2Z``. The method also computes
        generators of the subgroup at same time.

        INPUT:

        - ``G`` - intermediate subgroup (currently not implemented if different
          from SL(2,Z))

        - ``on_right`` - boolean (default: True) - if True return right coset
          enumeration, if False return left one.

        This is *extremely* slow in general.

        OUTPUT:

        - a list of coset representatives

        - a list of generators for the group

        - ``l`` - list of integers that correspond to the action of the
          standard parabolic element [[1,1],[0,1]] of `SL(2,\ZZ)` on the cosets
          of self.

        - ``s`` - list of integers that correspond to the action of the standard
          element of order `2` [[0,-1],[1,0]] on the cosets of self.

        EXAMPLES::

            sage: L = SL2Z([1,1,0,1])
            sage: S = SL2Z([0,-1,1,0])

            sage: G = Gamma(2)
            sage: reps, gens, l, s = G.todd_coxeter()
            sage: len(reps) == G.index()
            True
            sage: all(reps[i] * L * ~reps[l[i]] in G for i in range(6))
            True
            sage: all(reps[i] * S * ~reps[s[i]] in G for i in range(6))
            True

            sage: G = Gamma0(7)
            sage: reps, gens, l, s = G.todd_coxeter()
            sage: len(reps) == G.index()
            True
            sage: all(reps[i] * L * ~reps[l[i]] in G for i in range(8))
            True
            sage: all(reps[i] * S * ~reps[s[i]] in G for i in range(8))
            True

            sage: G = Gamma1(3)
            sage: reps, gens, l, s = G.todd_coxeter(on_right=False)
            sage: len(reps) == G.index()
            True
            sage: all(~reps[l[i]] * L * reps[i] in G for i in range(8))
            True
            sage: all(~reps[s[i]] * S * reps[i] in G for i in range(8))
            True

            sage: G = Gamma0(5)
            sage: reps, gens, l, s = G.todd_coxeter(on_right=False)
            sage: len(reps) == G.index()
            True
            sage: all(~reps[l[i]] * L * reps[i] in G for i in range(6))
            True
            sage: all(~reps[s[i]] * S * reps[i] in G for i in range(6))
            True
        """
        if G is None:
            G = SL2Z
        if G != SL2Z:
            raise NotImplementedError("Don't know how to compute coset reps for subgroups yet")

        id = SL2Z([1,0,0,1])
        l = SL2Z([1,1,0,1])
        s = SL2Z([0,-1,1,0])

        reps = [id]       # coset representatives
        reps_inv = {id:0} # coset representatives index

        l_wait_back = [id] # rep with no incoming s_edge
        s_wait_back = [id] # rep with no incoming l_edge
        l_wait = [id]      # rep with no outgoing l_edge
        s_wait = [id]      # rep with no outgoing s_edge

        l_edges = [None]    # edges for l
        s_edges = [None]    # edges for s

        gens = []

        while l_wait or s_wait:
            if l_wait:
                x = l_wait.pop(0)
                y = x
                not_end = True
                while not_end:
                    if on_right:
                        y = y*l
                    else:
                        y = l*y
                    for i in range(len(l_wait_back)):
                        v = l_wait_back[i]
                        if on_right:
                            yy = y*~v
                        else:
                            yy = ~v*y
                        if yy in self:
                            l_edges[reps_inv[x]] = reps_inv[v]
                            del l_wait_back[i]
                            if yy != id:
                                gens.append(self(yy))
                            not_end = False
                            break
                    else:
                        reps_inv[y] = len(reps)
                        l_edges[reps_inv[x]] = len(reps)
                        reps.append(y)
                        l_edges.append(None)
                        s_edges.append(None)
                        s_wait_back.append(y)
                        s_wait.append(y)
                    x = y

            if s_wait:
                x = s_wait.pop(0)
                y = x
                not_end = True
                while not_end:
                    if on_right:
                        y = y*s
                    else:
                        y = s*y
                    for i in range(len(s_wait_back)):
                        v = s_wait_back[i]
                        if on_right:
                            yy = y*~v
                        else:
                            yy = ~v*y
                        if yy in self:
                            s_edges[reps_inv[x]] = reps_inv[v]
                            del s_wait_back[i]
                            if yy != id:
                                gens.append(self(yy))
                            not_end = False
                            break
                    else:
                        reps_inv[y] = len(reps)
                        s_edges[reps_inv[x]] = len(reps)
                        reps.append(y)
                        l_edges.append(None)
                        s_edges.append(None)
                        l_wait_back.append(y)
                        l_wait.append(y)
                    x = y

        return reps, gens, l_edges, s_edges

    def nu2(self):
        r"""
        Return the number of orbits of elliptic points of order 2 for this
        arithmetic subgroup.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup().nu2()
            Traceback (most recent call last):
            ...
            NotImplementedError: Please implement _contains_sl2 for <class 'sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup_with_category'>
            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup.nu2(Gamma0(1105)) == 8
            True
        """

        # Subgroups not containing -1 have no elliptic points of order 2.

        if not self.is_even():
            return 0

        # Cheap trick: if self is a subgroup of something with no elliptic points,
        # then self has no elliptic points either.

        from .all import Gamma0, is_CongruenceSubgroup
        if is_CongruenceSubgroup(self):
            if self.is_subgroup(Gamma0(self.level())) and Gamma0(self.level()).nu2() == 0:
                return 0

        # Otherwise, the number of elliptic points is the number of g in self \
        # SL2Z such that the stabiliser of g * i in self is not trivial. (Note
        # that the points g*i for g in the coset reps are not distinct, but it
        # still works, since the failure of these points to be distinct happens
        # precisely when the preimages are not elliptic.)

        count = 0
        for g in self.coset_reps():
            if g * SL2Z([0,1,-1,0]) * (~g) in self:
                count += 1
        return count

    def nu3(self):
        r"""
        Return the number of orbits of elliptic points of order 3 for this
        arithmetic subgroup.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup().nu3()
            Traceback (most recent call last):
            ...
            NotImplementedError: Please implement _contains_sl2 for <class 'sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup_with_category'>
            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup.nu3(Gamma0(1729)) == 8
            True

        We test that a bug in handling of subgroups not containing -1 is fixed::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup.nu3(GammaH(7, [2]))
            2
        """

        # Cheap trick: if self is a subgroup of something with no elliptic points,
        # then self has no elliptic points either.

        from .all import Gamma0, is_CongruenceSubgroup
        if is_CongruenceSubgroup(self):
            if self.is_subgroup(Gamma0(self.level())) and Gamma0(self.level()).nu3() == 0:
                return 0

        count = 0
        for g in self.coset_reps():
            if g * SL2Z([0,1,-1,-1]) * (~g) in self:
                count += 1

        if self.is_even():
            return count
        else:
            return count // 2

    def is_abelian(self):
        r"""
        Return True if this arithmetic subgroup is abelian.

        Since arithmetic subgroups are always nonabelian, this always
        returns False.

        EXAMPLES::

            sage: SL2Z.is_abelian()
            False
            sage: Gamma0(3).is_abelian()
            False
            sage: Gamma1(12).is_abelian()
            False
            sage: GammaH(4, [3]).is_abelian()
            False
        """
        return False

    def is_finite(self):
        r"""
        Return True if this arithmetic subgroup is finite.

        Since arithmetic subgroups are always infinite, this always
        returns False.

        EXAMPLES::

            sage: SL2Z.is_finite()
            False
            sage: Gamma0(3).is_finite()
            False
            sage: Gamma1(12).is_finite()
            False
            sage: GammaH(4, [3]).is_finite()
            False
        """
        return False

    def is_subgroup(self, right):
        r"""
        Return True if self is a subgroup of right, and False otherwise. For
        generic arithmetic subgroups this is done by the absurdly slow
        algorithm of checking all of the generators of self to see if they are
        in right.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup().is_subgroup(SL2Z)
            Traceback (most recent call last):
            ...
            NotImplementedError: Please implement _contains_sl2 for <class 'sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup_with_category'>
            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup.is_subgroup(Gamma1(18), Gamma0(6))
            True
        """
        # ridiculously slow generic algorithm

        w = self.gens()
        for g in w:
            if not (g in right):
                return False
        return True

    def is_normal(self):
        r"""
        Return True precisely if this subgroup is a normal subgroup of SL2Z.

        EXAMPLES::

            sage: Gamma(3).is_normal()
            True
            sage: Gamma1(3).is_normal()
            False
        """
        for x in self.gens():
            for y in SL2Z.gens():
                if y*SL2Z(x)*(~y) not in self:
                    return False
        return True

    def is_odd(self):
        r"""
        Return True precisely if this subgroup does not contain the
        matrix -1.

        EXAMPLES::

            sage: SL2Z.is_odd()
            False
            sage: Gamma0(20).is_odd()
            False
            sage: Gamma1(5).is_odd()
            True
            sage: GammaH(11, [3]).is_odd()
            True
        """
        return not self.is_even()

    def is_even(self):
        r"""
        Return True precisely if this subgroup contains the matrix -1.

        EXAMPLES::

            sage: SL2Z.is_even()
            True
            sage: Gamma0(20).is_even()
            True
            sage: Gamma1(5).is_even()
            False
            sage: GammaH(11, [3]).is_even()
            False
        """
        return [-1, 0, 0, -1] in self

    def to_even_subgroup(self):
        r"""
        Return the smallest even subgroup of `SL(2, \ZZ)` containing self.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup().to_even_subgroup()
            Traceback (most recent call last):
            ...
            NotImplementedError: Please implement _contains_sl2 for <class 'sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup_with_category'>
        """
        if self.is_even():
            return self
        else:
            raise NotImplementedError

    def order(self):
        r"""
        Return the number of elements in this arithmetic subgroup.

        Since arithmetic subgroups are always infinite, this always returns
        infinity.

        EXAMPLES::

            sage: SL2Z.order()
            +Infinity
            sage: Gamma0(5).order()
            +Infinity
            sage: Gamma1(2).order()
            +Infinity
            sage: GammaH(12, [5]).order()
            +Infinity
        """
        from sage.rings.infinity import infinity
        return infinity

    def reduce_cusp(self, c):
        r"""
        Given a cusp `c \in \mathbb{P}^1(\QQ)`, return the unique reduced cusp
        equivalent to c under the action of self, where a reduced cusp is an
        element `\tfrac{r}{s}` with r,s coprime non-negative integers, s as
        small as possible, and r as small as possible for that s.

        NOTE: This function should be overridden by all subclasses.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup().reduce_cusp(1/4)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def cusps(self, algorithm='default'):
        r"""
        Return a sorted list of inequivalent cusps for self, i.e. a set of
        representatives for the orbits of self on `\mathbb{P}^1(\QQ)`.

        These should be returned in a reduced form where this makes sense.

        INPUT:

        - ``algorithm`` -- which algorithm to use to compute the cusps of self.
          ``'default'`` finds representatives for a known complete set of
          cusps. ``'modsym'`` computes the boundary map on the space of weight
          two modular symbols associated to self, which finds the cusps for
          self in the process.

        EXAMPLES::

            sage: Gamma0(36).cusps()
            [0, 1/18, 1/12, 1/9, 1/6, 1/4, 1/3, 5/12, 1/2, 2/3, 5/6, Infinity]
            sage: Gamma0(36).cusps(algorithm='modsym') == Gamma0(36).cusps()
            True
            sage: GammaH(36, [19,29]).cusps() == Gamma0(36).cusps()
            True
            sage: Gamma0(1).cusps()
            [Infinity]
        """
        try:
            return copy(self._cusp_list[algorithm])
        except (AttributeError, KeyError):
            self._cusp_list = {}

        from .congroup_sl2z import is_SL2Z
        if algorithm == 'default':
            if is_SL2Z(self):
                s = [Cusp(1, 0)]
            else:
                s = self._find_cusps()
        elif algorithm == 'modsym':
            s = sorted(self.reduce_cusp(c)
                       for c in self.modular_symbols().cusps())
        else:
            raise ValueError("unknown algorithm: %s" % algorithm)

        self._cusp_list[algorithm] = s
        return copy(s)

    def _find_cusps(self):
        r"""
        Calculate a list of inequivalent cusps.

        EXAMPLES::

            sage: sage.modular.arithgroup.congroup_generic.CongruenceSubgroup(5)._find_cusps()
            Traceback (most recent call last):
            ...
            NotImplementedError

        NOTE: There is a generic algorithm implemented at the top level that
        uses the coset representatives of self. This is *very slow* and for all
        the standard congruence subgroups there is a quicker way of doing it,
        so this should usually be overridden in subclasses; but it doesn't have
        to be.
        """
        i = Cusp([1,0])
        L = [i]
        for a in self.coset_reps():
            ai = i.apply([a.a(), a.b(), a.c(), a.d()])
            new = 1
            for v in L:
                if self.are_equivalent(ai, v):
                    new = 0
                    break
            if new == 1:
                L.append(ai)
        return L

    def are_equivalent(self, x, y, trans = False):
        r"""
        Test whether or not cusps x and y are equivalent modulo self.  If self
        has a reduce_cusp() method, use that; otherwise do a slow explicit
        test.

        If trans = False, returns True or False. If trans = True, then return
        either False or an element of self mapping x onto y.

        EXAMPLES::

            sage: Gamma0(7).are_equivalent(Cusp(1/3), Cusp(0), trans=True)
            [  3  -1]
            [-14   5]
            sage: Gamma0(7).are_equivalent(Cusp(1/3), Cusp(1/7))
            False
        """
        x = Cusp(x)
        y = Cusp(y)
        if not trans:
            try:
                xr = self.reduce_cusp(x)
                yr = self.reduce_cusp(y)
                if xr != yr:
                    return False
                if xr == yr:
                    return True
            except NotImplementedError:
                pass

        vx = lift_to_sl2z(x.numerator(),x.denominator(), 0)
        dx = SL2Z([vx[2], -vx[0], vx[3], -vx[1]])
        vy = lift_to_sl2z(y.numerator(),y.denominator(), 0)
        dy = SL2Z([vy[2], -vy[0], vy[3], -vy[1]])

        for i in range(self.index()):
            # Note that the width of any cusp is bounded above by the index of self.
            # If self is congruence, then the level of self is a much better bound, but
            # this method is written to work with non-congruence subgroups as well,
            if dy * SL2Z([1,i,0,1])*(~dx) in self:
                if trans:
                    return dy * SL2Z([1,i,0,1]) * ~dx
                else:
                    return True
            elif (self.is_odd() and dy * SL2Z([-1,-i,0,-1]) * ~dx in self):
                if trans:
                    return dy * SL2Z([-1,-i,0,-1]) * ~dx
                else:
                    return True
        return False

    def cusp_data(self, c):
        r"""
        Return a triple (g, w, t) where g is an element of self generating the
        stabiliser of the given cusp, w is the width of the cusp, and t is 1 if
        the cusp is regular and -1 if not.

        EXAMPLES::

            sage: Gamma1(4).cusp_data(Cusps(1/2))
            (
            [ 1 -1]
            [ 4 -3], 1, -1
            )
        """
        c = Cusp(c)

        # first find an element of SL2Z sending infinity to the given cusp
        w = lift_to_sl2z(c.denominator(), c.numerator(), 0)
        g = SL2Z([w[3], w[1], w[2],w[0]])

        for d in range(1,1+self.index()):
            if g * SL2Z([1,d,0,1]) * (~g) in self:
                return (g * SL2Z([1,d,0,1]) * (~g), d, 1)
            elif g * SL2Z([-1,-d,0,-1]) * (~g) in self:
                return (g * SL2Z([-1,-d,0,-1]) * (~g), d, -1)
        raise ArithmeticError("Can't get here!")

    def is_regular_cusp(self, c):
        r"""
        Return True if the orbit of the given cusp is a regular cusp for self,
        otherwise False. This is automatically true if -1 is in self.

        EXAMPLES::

            sage: Gamma1(4).is_regular_cusp(Cusps(1/2))
            False
            sage: Gamma1(4).is_regular_cusp(Cusps(oo))
            True
        """
        if self.is_even():
            return True
        return (self.cusp_data(c)[2] == 1)

    def cusp_width(self, c):
        r"""
        Return the width of the orbit of cusps represented by c.

        EXAMPLES::

            sage: Gamma0(11).cusp_width(Cusps(oo))
            1
            sage: Gamma0(11).cusp_width(0)
            11
            sage: [Gamma0(100).cusp_width(c) for c in Gamma0(100).cusps()]
            [100, 1, 4, 1, 1, 1, 4, 25, 1, 1, 4, 1, 25, 4, 1, 4, 1, 1]
        """
        return self.cusp_data(c)[1]

    def index(self):
        r"""
        Return the index of self in the full modular group.

        EXAMPLES::

            sage: Gamma0(17).index()
            18
            sage: sage.modular.arithgroup.congroup_generic.CongruenceSubgroup(5).index()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """

        return len(list(self.coset_reps()))

    def generalised_level(self):
        r"""
        Return the generalised level of self, i.e. the least common multiple of
        the widths of all cusps.

        If self is *even*, Wohlfart's theorem tells us that this is equal to
        the (conventional) level of self when self is a congruence subgroup.
        This can fail if self is odd, but the actual level is at most twice the
        generalised level. See the paper by Kiming, Schuett and Verrill for
        more examples.

        EXAMPLES::

            sage: Gamma0(18).generalised_level()
            18
            sage: sage.modular.arithgroup.arithgroup_perm.HsuExample18().generalised_level()
            24

        In the following example, the actual level is twice the generalised
        level. This is the group `G_2` from Example 17 of K-S-V.

        ::

            sage: G = CongruenceSubgroup(8, [ [1,1,0,1], [3,-1,4,-1] ])
            sage: G.level()
            8
            sage: G.generalised_level()
            4
        """
        return lcm([self.cusp_width(c) for c in self.cusps()])

    def projective_index(self):
        r"""
        Return the index of the image of self in `{\rm PSL}_2(\ZZ)`. This is equal
        to the index of self if self contains -1, and half of this otherwise.

        This is equal to the degree of the natural map from the modular curve
        of self to the `j`-line.

        EXAMPLES::

            sage: Gamma0(5).projective_index()
            6
            sage: Gamma1(5).projective_index()
            12
        """

        if self.is_even():
            return self.index()
        else:
            return self.index() // 2

    def is_congruence(self):
        r"""
        Return True if self is a congruence subgroup.

        EXAMPLES::

            sage: Gamma0(5).is_congruence()
            True
            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup().is_congruence()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """

        raise NotImplementedError

    def genus(self):
        r"""
        Return the genus of the modular curve of ``self``.

        EXAMPLES::

            sage: Gamma1(5).genus()
            0
            sage: Gamma1(31).genus()
            26
            sage: from sage.modular.dims import dimension_cusp_forms
            sage: Gamma1(157).genus() == dimension_cusp_forms(Gamma1(157), 2)
            True
            sage: GammaH(7, [2]).genus()
            0
            sage: [Gamma0(n).genus() for n in [1..23]]
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 2, 2]
            sage: [n for n in [1..200] if Gamma0(n).genus() == 1]
            [11, 14, 15, 17, 19, 20, 21, 24, 27, 32, 36, 49]
        """
        return ZZ(1 + (self.projective_index()) / ZZ(12)  - (self.nu2())/ZZ(4) - (self.nu3())/ZZ(3) - self.ncusps()/ZZ(2))

    def farey_symbol(self):
        r"""
        Return the Farey symbol associated to this subgroup. See the
        :mod:`~sage.modular.arithgroup.farey_symbol` module for more
        information.

        EXAMPLES::

            sage: Gamma1(4).farey_symbol()
            FareySymbol(Congruence Subgroup Gamma1(4))
        """
        from .farey_symbol import Farey
        return Farey(self)

    @cached_method
    def generators(self, algorithm="farey"):
        r"""
        Return a list of generators for this congruence subgroup. The result is cached.

        INPUT:

        - ``algorithm`` (string): either ``farey`` or ``todd-coxeter``.

        If ``algorithm`` is set to ``"farey"``, then the generators will be
        calculated using Farey symbols, which will always return a *minimal*
        generating set. See :mod:`~sage.modular.arithgroup.farey_symbol` for
        more information.

        If ``algorithm`` is set to ``"todd-coxeter"``, a simpler algorithm
        based on Todd-Coxeter enumeration will be used. This is *exceedingly*
        slow for general subgroups, and the list of generators will be far from
        minimal (indeed it may contain repetitions).

        EXAMPLES::

            sage: Gamma(2).generators()
            [
            [1 2]  [ 3 -2]  [-1  0]
            [0 1], [ 2 -1], [ 0 -1]
            ]
            sage: Gamma(2).generators(algorithm="todd-coxeter")
            [
            [1 2]  [-1  0]  [ 1  0]  [-1  0]  [-1  2]  [-1  0]  [1 0]
            [0 1], [ 0 -1], [-2  1], [ 0 -1], [-2  3], [ 2 -1], [2 1]
            ]
        """
        if algorithm=="farey":
            return self.farey_symbol().generators()
        elif algorithm == "todd-coxeter":
            return self.todd_coxeter()[1]
        else:
            raise ValueError("Unknown algorithm '%s' (should be either 'farey' or 'todd-coxeter')" % algorithm)

    def gens(self, *args, **kwds):
        r"""
        Return a tuple of generators for this congruence subgroup.

        The generators need not be minimal. For arguments, see :meth:`~generators`.

        EXAMPLES::

            sage: SL2Z.gens()
            (
            [ 0 -1]  [1 1]
            [ 1  0], [0 1]
            )
        """
        return tuple(self.generators(*args, **kwds))

    def gen(self, i):
        r"""
        Return the i-th generator of self, i.e. the i-th element of the
        tuple self.gens().

        EXAMPLES::

            sage: SL2Z.gen(1)
            [1 1]
            [0 1]
        """
        return self.generators()[i]

    def ngens(self):
        r"""
        Return the size of the minimal generating set of self returned by
        :meth:`generators`.

        EXAMPLES::

            sage: Gamma0(22).ngens()
            8
            sage: Gamma1(14).ngens()
            13
            sage: GammaH(11, [3]).ngens()
            3
            sage: SL2Z.ngens()
            2
        """
        return len(self.generators())

    def ncusps(self):
        r"""
        Return the number of cusps of this arithmetic subgroup. This is
        provided as a separate function since for dimension formulae in even
        weight all we need to know is the number of cusps, and this can be
        calculated very quickly, while enumerating all cusps is much slower.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup.ncusps(Gamma0(7))
            2
        """

        return ZZ(len(self.cusps()))

    def nregcusps(self):
        r"""
        Return the number of cusps of self that are "regular", i.e. their
        stabiliser has a generator with both eigenvalues +1 rather than -1. If
        the group contains -1, every cusp is clearly regular.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup.nregcusps(Gamma1(4))
            2
        """
        return self.ncusps() - self.nirregcusps()

    def nirregcusps(self):
        r"""
        Return the number of cusps of self that are "irregular", i.e. their
        stabiliser can only be generated by elements with both eigenvalues -1
        rather than +1. If the group contains -1, every cusp is clearly
        regular.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup.nirregcusps(Gamma1(4))
            1
        """
        if self.is_even():
            return 0
        else:
            return ZZ(len([c for c in self.cusps() if not self.is_regular_cusp(c)]))

    def dimension_modular_forms(self, k=2):
        r"""
        Return the dimension of the space of weight k modular forms for this
        group. This is given by a standard formula in terms of k and various
        invariants of the group; see Diamond + Shurman, "A First Course in
        Modular Forms", section 3.5 and 3.6. If k is not given, defaults to k =
        2.

        For dimensions of spaces of modular forms with character for Gamma1, use
        the dimension_modular_forms method of the Gamma1 class, or the standalone
        function dimension_modular_forms().

        For weight 1 modular forms this generic implementation only works in
        cases where one can prove solely via Riemann-Roch theory that there
        aren't any cusp forms (i.e. when the number of regular cusps is
        strictly greater than the degree of the canonical divisor). Otherwise a
        NotImplementedError is raised.

        EXAMPLES::

            sage: Gamma1(31).dimension_modular_forms(2)
            55
            sage: Gamma1(3).dimension_modular_forms(1)
            1
            sage: Gamma1(4).dimension_modular_forms(1) # irregular cusp
            1
            sage: Gamma(13).dimension_modular_forms(1)
            Traceback (most recent call last):
            ...
            NotImplementedError: Computation of dimensions of weight 1 cusp forms spaces not implemented in general
        """
        return self.dimension_cusp_forms(k) + self.dimension_eis(k)

    def dimension_cusp_forms(self, k=2):
        r"""
        Return the dimension of the space of weight k cusp forms for this
        group. For `k \ge 2`, this is given by a standard formula in terms of k
        and various invariants of the group; see Diamond + Shurman, "A First
        Course in Modular Forms", section 3.5 and 3.6. If k is not given,
        default to k = 2.

        For dimensions of spaces of cusp forms with character for Gamma1, use
        the dimension_cusp_forms method of the Gamma1 class, or the standalone
        function dimension_cusp_forms().

        For weight 1 cusp forms this generic implementation only works in cases
        where one can prove solely via Riemann-Roch theory that there aren't
        any cusp forms (i.e. when the number of regular cusps is strictly
        greater than the degree of the canonical divisor). Otherwise a
        NotImplementedError is raised.

        EXAMPLES::

            sage: Gamma1(31).dimension_cusp_forms(2)
            26
            sage: Gamma(5).dimension_cusp_forms(1)
            0
            sage: Gamma1(4).dimension_cusp_forms(1) # irregular cusp
            0
            sage: Gamma(13).dimension_cusp_forms(1)
            Traceback (most recent call last):
            ...
            NotImplementedError: Computation of dimensions of weight 1 cusp forms spaces not implemented in general
        """
        k = ZZ(k)
        if k <= 0:
            return ZZ(0)

        if not (k % 2):
            # k even

            if k == 2:
                return self.genus()

            else:
                return (k-1) * (self.genus() - 1) + (k // ZZ(4))*self.nu2() + (k // ZZ(3))*self.nu3() + (k // ZZ(2) - 1)*self.ncusps()

        else:
            # k odd

            if self.is_even():
                return ZZ(0)

            else:
                e_reg = self.nregcusps()
                e_irr = self.nirregcusps()

                if k > 1:
                    return (k-1)*(self.genus()-1) + (k // ZZ(3)) * self.nu3() + (k-2)/ZZ(2) * e_reg + (k-1)/ZZ(2) * e_irr
                else:
                    if e_reg > 2*self.genus() - 2:
                        return ZZ(0)
                    else:
                        raise NotImplementedError("Computation of dimensions of weight 1 cusp forms spaces not implemented in general")

    def dimension_eis(self, k=2):
        r"""
        Return the dimension of the space of weight k Eisenstein series for
        this group, which is a subspace of the space of modular forms
        complementary to the space of cusp forms.

        INPUT:

        - ``k`` - an integer (default 2).

        EXAMPLES::

            sage: GammaH(33,[2]).dimension_eis()
            7
            sage: GammaH(33,[2]).dimension_eis(3)
            0
            sage: GammaH(33, [2,5]).dimension_eis(2)
            3
            sage: GammaH(33, [4]).dimension_eis(1)
            4
        """
        if k < 0:
            return ZZ(0)
        if k == 0:
            return ZZ(1)

        if not (k % 2): # k even
            if k > 2:
                return self.ncusps()
            else: # k = 2
                return self.ncusps() - 1

        else: # k odd
            if self.is_even():
                return 0
            if k > 1:
                return self.nregcusps()
            else: # k = 1
                return ZZ(self.nregcusps()/ ZZ(2))

    def as_permutation_group(self):
        r"""
        Return self as an arithmetic subgroup defined in terms of the
        permutation action of `SL(2,\ZZ)` on its right cosets.

        This method uses Todd-Coxeter enumeration (via the method
        :meth:`~todd_coxeter`) which can be extremely slow for arithmetic
        subgroups with relatively large index in `SL(2,\ZZ)`.

        EXAMPLES::

            sage: G = Gamma(3)
            sage: P = G.as_permutation_group(); P
            Arithmetic subgroup of index 24
            sage: G.ncusps() == P.ncusps()
            True
            sage: G.nu2() == P.nu2()
            True
            sage: G.nu3() == P.nu3()
            True
            sage: G.an_element() in P
            True
            sage: P.an_element() in G
            True
        """
        _,_,l_edges,s2_edges=self.todd_coxeter()
        n = len(l_edges)
        s3_edges = [None] * n
        r_edges = [None] * n
        for i in range(n):
            ii = s2_edges[l_edges[i]]
            s3_edges[ii] = i
            r_edges[ii] = s2_edges[i]
        if self.is_even():
            from sage.modular.arithgroup.arithgroup_perm import EvenArithmeticSubgroup_Permutation
            g=EvenArithmeticSubgroup_Permutation(S2=s2_edges,S3=s3_edges,L=l_edges,R=r_edges)
        else:
            from sage.modular.arithgroup.arithgroup_perm import OddArithmeticSubgroup_Permutation
            g=OddArithmeticSubgroup_Permutation(S2=s2_edges,S3=s3_edges,L=l_edges,R=r_edges)
        g.relabel()
        return g

    def sturm_bound(self, weight=2):
        r"""
        Returns the Sturm bound for modular forms of the given weight and level
        this subgroup.

        INPUT:

        -  ``weight`` - an integer `\geq 2` (default: 2)

        EXAMPLES::

            sage: Gamma0(11).sturm_bound(2)
            2
            sage: Gamma0(389).sturm_bound(2)
            65
            sage: Gamma0(1).sturm_bound(12)
            1
            sage: Gamma0(100).sturm_bound(2)
            30
            sage: Gamma0(1).sturm_bound(36)
            3
            sage: Gamma0(11).sturm_bound()
            2
            sage: Gamma0(13).sturm_bound()
            3
            sage: Gamma0(16).sturm_bound()
            4
            sage: GammaH(16,[13]).sturm_bound()
            8
            sage: GammaH(16,[15]).sturm_bound()
            16
            sage: Gamma1(16).sturm_bound()
            32
            sage: Gamma1(13).sturm_bound()
            28
            sage: Gamma1(13).sturm_bound(5)
            70

        FURTHER DETAILS: This function returns a positive integer
        `n` such that the Hecke operators
        `T_1,\ldots, T_n` acting on *cusp forms* generate the
        Hecke algebra as a `\ZZ`-module when the character
        is trivial or quadratic. Otherwise, `T_1,\ldots,T_n`
        generate the Hecke algebra at least as a
        `\ZZ[\varepsilon]`-module, where
        `\ZZ[\varepsilon]` is the ring generated by the
        values of the Dirichlet character `\varepsilon`.
        Alternatively, this is a bound such that if two cusp forms
        associated to this space of modular symbols are congruent modulo
        `(\lambda, q^n)`, then they are congruent modulo
        `\lambda`.

        REFERENCES:

        - See the Agashe-Stein appendix to Lario and Schoof,
          *Some computations with Hecke rings and deformation rings*,
          Experimental Math., 11 (2002), no. 2, 303-311.

        - This result originated in the paper Sturm,
          *On the congruence of modular forms*,
          Springer LNM 1240, 275-280, 1987.

        REMARK: Kevin Buzzard pointed out to me (William Stein) in Fall
        2002 that the above bound is fine for `\Gamma_1(N)` with
        character, as one sees by taking a power of `f`. More
        precisely, if `f \cong 0 \pmod{p}` for first
        `s` coefficients, then `f^r \cong 0 \pmod{p}` for
        first `sr` coefficients. Since the weight of `f^r`
        is `r\cdot k(f)`, it follows that if
        `s \geq b`, where `b` is the Sturm bound for
        `\Gamma_0(N)` at weight `k(f)`, then `f^r`
        has valuation large enough to be forced to be `0` at
        `r*k(f)` by Sturm bound (which is valid if we choose
        `r` correctly). Thus `f \cong 0 \pmod{p}`.
        Conclusion: For `\Gamma_1(N)` with fixed character, the
        Sturm bound is *exactly* the same as for `\Gamma_0(N)`.

        A key point is that we are finding
        `\ZZ[\varepsilon]` generators for the Hecke algebra
        here, not `\ZZ`-generators. So if one wants
        generators for the Hecke algebra over `\ZZ`, this
        bound must be suitably modified (and I'm not sure what the
        modification is).

        AUTHORS:

        - William Stein
        """
        return ZZ((self.index() * weight / ZZ(12)).ceil())
