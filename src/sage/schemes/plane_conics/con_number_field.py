# sage.doctest: needs sage.rings.number_field
r"""
Projective plane conics over a number field

AUTHORS:

- Marco Streng (2010-07-20)
"""
# ****************************************************************************
#       Copyright (C) 2009/2010 Marco Streng <marco.streng@gmail.com>
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.rational_field import RationalField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from .con_field import ProjectiveConic_field


class ProjectiveConic_number_field(ProjectiveConic_field):
    r"""
    Create a projective plane conic curve over a number field.
    See ``Conic`` for full documentation.

    EXAMPLES::

        sage: K.<a> = NumberField(x^3 - 2, 'a')
        sage: P.<X, Y, Z> = K[]
        sage: Conic(X^2 + Y^2 - a*Z^2)
        Projective Conic Curve over Number Field in a with defining polynomial x^3 - 2
         defined by X^2 + Y^2 + (-a)*Z^2

    TESTS::

        sage: K.<a> = NumberField(x^3 - 3, 'a')
        sage: Conic([a, 1, -1])._test_pickling()
    """
    def __init__(self, A, f):
        r"""
        See ``Conic`` for full documentation.

        EXAMPLES::

            sage: Conic([1, 1, 1])
            Projective Conic Curve over Rational Field defined by x^2 + y^2 + z^2
        """
        ProjectiveConic_field.__init__(self, A, f)

        # a single prime such that self has no point over the completion
        self._local_obstruction = None
        # all finite primes such that self has no point over the completion
        self._finite_obstructions = None
        # all infinite primes such that self has no point over the completion
        self._infinite_obstructions = None

    def has_rational_point(self, point=False, obstruction=False,
                           algorithm='default', read_cache=True):
        r"""
        Return ``True`` if and only if ``self`` has a point
        defined over its base field `B`.

        If ``point`` and ``obstruction`` are both False (default),
        then the output is a boolean ``out`` saying whether ``self``
        has a rational point.

        If ``point`` or ``obstruction`` is ``True``, then the output is
        a pair ``(out, S)``, where ``out`` is as above and:

        - if ``point`` is ``True`` and ``self`` has a rational point,
          then ``S`` is a rational point,

        - if ``obstruction`` is ``True``, ``self`` has no rational point,
          then ``S`` is a prime or infinite place of `B` such that no
          rational point exists over the completion at ``S``.

        Points and obstructions are cached whenever they are found.
        Cached information is used for the output if available, but only
        if ``read_cache`` is ``True``.

        ALGORITHM:

        The parameter ``algorithm``
        specifies the algorithm to be used:

        - ``'rnfisnorm'`` -- use PARI's ``rnfisnorm``
          (cannot be combined with ``obstruction = True``)

        - ``'local'`` -- check if a local solution exists for all primes
          and infinite places of `B` and apply the Hasse principle
          (cannot be combined with ``point = True``)

        - ``'default'`` -- use algorithm ``'rnfisnorm'`` first.
          Then, if no point exists and obstructions are requested, use
          algorithm ``'local'`` to find an obstruction.

        - ``'magma'`` (requires Magma to be installed) --
          delegates the task to the Magma computer algebra
          system

        EXAMPLES:

        An example over `\QQ` ::

            sage: C = Conic(QQ, [1, 113922743, -310146482690273725409])
            sage: C.has_rational_point(point=True)
            (True, (-76842858034579/5424 : -5316144401/5424 : 1))
            sage: C.has_rational_point(algorithm='local', read_cache=False)
            True

        Examples over number fields::

            sage: K.<i> = QuadraticField(-1)
            sage: C = Conic(K, [1, 3, -5])
            sage: C.has_rational_point(point=True, obstruction=True)
            (False, Fractional ideal (i + 2))
            sage: C.has_rational_point(algorithm='rnfisnorm')
            False
            sage: C.has_rational_point(algorithm='rnfisnorm', obstruction=True,
            ....:                      read_cache=False)
            Traceback (most recent call last):
            ...
            ValueError: Algorithm rnfisnorm cannot be combined with
            obstruction = True in has_rational_point

            sage: P.<x> = QQ[]
            sage: L.<b> = NumberField(x^3 - 5)
            sage: C = Conic(L, [1, 2, -3])
            sage: C.has_rational_point(point=True, algorithm='rnfisnorm')
            (True, (-5/3 : 1/3 : 1))

            sage: K.<a> = NumberField(x^4+2)
            sage: Conic(QQ, [4,5,6]).has_rational_point()
            False
            sage: Conic(K, [4,5,6]).has_rational_point()
            True
            sage: Conic(K, [4,5,6]).has_rational_point(algorithm='magma',   # optional - magma
            ....:                                      read_cache=False)
            True

            sage: P.<a> = QuadraticField(2)
            sage: C = Conic(P, [1,1,1])
            sage: C.has_rational_point()
            False
            sage: C.has_rational_point(point=True)
            (False, None)
            sage: C.has_rational_point(obstruction=True)
            (False,
             Ring morphism:
               From: Number Field in a with defining polynomial x^2 - 2
                     with a = 1.414213562373095?
               To:   Algebraic Real Field
               Defn: a |--> -1.414213562373095?)
            sage: C.has_rational_point(point=True, obstruction=True)
            (False,
             Ring morphism:
               From: Number Field in a with defining polynomial x^2 - 2
                     with a = 1.414213562373095?
               To:   Algebraic Real Field
               Defn: a |--> -1.414213562373095?)

        TESTS:

        Create a bunch of conics over number fields and check whether
        ``has_rational_point`` runs without errors for algorithms
        ``'rnfisnorm'`` and ``'local'``. Check if all points returned are
        valid. If Magma is available, then also check if the output agrees with
        Magma. ::

            sage: P.<X> = QQ[]
            sage: Q = P.fraction_field()
            sage: c = [1, X/2, 1/X]
            sage: l = Sequence(cartesian_product_iterator([c for i in range(3)]))
            sage: l = l + [[X, 1, 1, 1, 1, 1]] + [[X, 1/5, 1, 1, 2, 1]]
            sage: K.<a> = QuadraticField(-23)
            sage: L.<b> = QuadraticField(19)
            sage: M.<c> = NumberField(X^3+3*X+1)
            sage: m = [[Q(b)(F.gen()) for b in a] for a in l for F in [K, L, M]]
            sage: d = []
            sage: c = []
            sage: c = [Conic(a) for a in m if a != [0,0,0]]
            sage: d = [C.has_rational_point(algorithm='rnfisnorm', point=True) for C in c] # long time: 3.3 seconds
            sage: all(c[k].defining_polynomial()(Sequence(d[k][1])) == 0 for k in range(len(d)) if d[k][0])
            True
            sage: [C.has_rational_point(algorithm='local', read_cache=False) for C in c] == [o[0] for o in d] # long time: 5 seconds
            True
            sage: [C.has_rational_point(algorithm='magma', read_cache=False) for C in c] == [o[0] for o in d] # long time: 3 seconds, optional - magma
            True

        Create a bunch of conics that are known to have rational points
        already over `\QQ` and check if points are found by
        ``has_rational_point``. ::

            sage: l = Sequence(cartesian_product_iterator([[-1, 0, 1] for i in range(3)]))
            sage: K.<a> = QuadraticField(-23)
            sage: L.<b> = QuadraticField(19)
            sage: M.<c> = NumberField(x^5+3*x+1)
            sage: m = [[F(b) for b in a] for a in l for F in [K, L, M]]
            sage: c = [Conic(a) for a in m if a != [0,0,0] and a != [1,1,1] and a != [-1,-1,-1]]
            sage: assert all(C.has_rational_point(algorithm='rnfisnorm') for C in c)
            sage: assert all(C.defining_polynomial()(Sequence(C.has_rational_point(point=True)[1])) == 0 for C in c)
            sage: assert all(C.has_rational_point(algorithm='local', read_cache=False) for C in c) # long time: 1 second
        """
        if read_cache:
            if self._rational_point is not None:
                # a rational point is already known, return True
                if point or obstruction:
                    return True, self._rational_point
                else:
                    return True
            if self._local_obstruction is not None:
                # a local obstruction is already known, return False
                if point or obstruction:
                    return False, self._local_obstruction
                else:
                    return False
            # `_(in)finite_obstructions` is ``None`` if the cache is empty,
            # so we explicitly check against a list:
            if (not point) and self._finite_obstructions == [] and \
               self._infinite_obstructions == []:
                # list of local obstructions is computed and empty, return True
                if obstruction:
                    return True, None
                return True

        # cache is empty or user specifically asks not to check the cache

        if self.has_singular_point():
            if point:
                return self.has_singular_point(point=True)
            if obstruction:
                return True, None
            return True
        B = self.base_ring()

        if algorithm == 'default':
            ret = self.has_rational_point(point=True, obstruction=False,
                                          algorithm='rnfisnorm',
                                          read_cache=False)
            if ret[0]:
                if point or obstruction:
                    return ret
                return True
            if obstruction:
                ret = self.has_rational_point(point=False, obstruction=True,
                                              algorithm='local',
                                              read_cache=False)
                if ret[0]:
                    raise RuntimeError("Outputs of algorithms in "
                                       "has_rational_point disagree "
                                       "for conic %s" % self)
                return ret
            if point:
                return False, None
            return False

        if algorithm == 'local':
            if point:
                raise ValueError("Algorithm 'local' cannot be combined "
                                 "with point = True in has_rational_point")
            obs = self.local_obstructions(infinite=True, finite=False,
                                          read_cache=read_cache)
            if obs:
                if obstruction:
                    return False, obs[0]
                return False
            obs = self.local_obstructions(read_cache=read_cache)
            if not obs:
                if obstruction:
                    return True, None
                return True
            if obstruction:
                return False, obs[0]
            return False
        if algorithm == 'rnfisnorm':
            from sage.modules.free_module_element import vector
            if obstruction:
                raise ValueError("Algorithm rnfisnorm cannot be combined "
                                 "with obstruction = True in "
                                 "has_rational_point")
            D, T = self.diagonal_matrix()
            abc = [D[0, 0], D[1, 1], D[2, 2]]
            for j in range(3):
                if abc[j] == 0:
                    pt = self.point(T * vector({2: 0, j: 1}))
                    if point or obstruction:
                        return True, pt
                    return True
            if (-abc[1] / abc[0]).is_square():
                pt = self.point(T * vector([(-abc[1] / abc[0]).sqrt(), 1, 0]))
                if point or obstruction:
                    return True, pt
                return True
            if (-abc[2] / abc[0]).is_square():
                pt = self.point(T * vector([(-abc[2] / abc[0]).sqrt(), 0, 1]))
                if point or obstruction:
                    return True, pt
                return True
            if isinstance(B, RationalField):
                K = B
                [KtoB, BtoK] = [K.hom(K) for i in range(2)]
            else:
                K = B.absolute_field('Y')
                [KtoB, BtoK] = K.structure()
            X = PolynomialRing(K, 'X').gen()
            d = BtoK(-abc[1] / abc[0])
            den = d.denominator()
            L = K.extension(X**2 - d * den**2, names='y')
            isnorm = BtoK(-abc[2] / abc[0]).is_norm(L, element=True)
            if isnorm[0]:

                pt = self.point(T * vector([KtoB(isnorm[1][0]),
                                          KtoB(isnorm[1][1] * den), 1]))
                if point:
                    return True, pt
                return True
            if point:
                return False, None
            return False
        if algorithm == 'qfsolve':
            raise TypeError("Algorithm qfsolve in has_rational_point only "
                            "for conics over QQ, not over %s" % B)
        if obstruction:
            raise ValueError("Invalid combination: obstruction=True and "
                             "algorithm=%s" % algorithm)

        return ProjectiveConic_field.has_rational_point(self, point=point,
                           algorithm=algorithm, read_cache=False)

    def is_locally_solvable(self, p):
        r"""
        Return ``True`` if and only if ``self`` has a solution over the
        completion of the base field `B` of ``self`` at ``p``. Here ``p``
        is a finite prime or infinite place of `B`.

        EXAMPLES::

            sage: P.<x> = QQ[]
            sage: K.<a> = NumberField(x^3 + 5)
            sage: C = Conic(K, [1, 2, 3 - a])
            sage: [p1, p2] = K.places()
            sage: C.is_locally_solvable(p1)
            False

            sage: C.is_locally_solvable(p2)
            True

            sage: f = (2*K).factor()
            sage: C.is_locally_solvable(f[0][0])
            True

            sage: C.is_locally_solvable(f[1][0])
            False
        """
        D, T = self.diagonal_matrix()
        abc = [D[j, j] for j in range(3)]
        for a in abc:
            if a == 0:
                return True
        a = -abc[0] / abc[2]
        b = -abc[1] / abc[2]

        ret = self.base_ring().hilbert_symbol(a, b, p)

        if ret == -1:
            if self._local_obstruction is None:
                from sage.categories.map import Map
                from sage.categories.rings import Rings
                from sage.rings.qqbar import AA
                from sage.rings.real_lazy import RLF

                if not (isinstance(p, Map) and p.category_for().is_subcategory(Rings())) or p.codomain() is AA or p.codomain() is RLF:
                    self._local_obstruction = p
            return False

        return True

    def local_obstructions(self, finite=True, infinite=True, read_cache=True):
        r"""
        Return the sequence of finite primes and/or infinite places
        such that ``self`` is locally solvable at those primes and places.

        If the base field is `\QQ`, then the infinite place is denoted `-1`.

        The parameters ``finite`` and ``infinite`` (both ``True`` by default) are
        used to specify whether to look at finite and/or infinite places.
        Note that ``finite = True`` involves factorization of the determinant
        of ``self``, hence may be slow.

        Local obstructions are cached. The parameter ``read_cache``
        specifies whether to look at the cache before computing anything.

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: Conic(K, [1, 2, 3]).local_obstructions()
            []

            sage: L.<a> = QuadraticField(5)
            sage: Conic(L, [1, 2, 3]).local_obstructions()
            [Ring morphism:
               From: Number Field in a with defining polynomial x^2 - 5
                     with a = 2.236067977499790?
               To:   Algebraic Real Field
               Defn: a |--> -2.236067977499790?,
             Ring morphism:
               From: Number Field in a with defining polynomial x^2 - 5
                     with a = 2.236067977499790?
               To:   Algebraic Real Field
               Defn: a |--> 2.236067977499790?]
        """
        obs0 = []
        obs1 = []
        B = self.base_ring()
        if infinite:
            if read_cache and self._infinite_obstructions is not None:
                obs0 = self._infinite_obstructions
            else:
                from sage.rings.qqbar import AA
                for b in B.embeddings(AA):
                    if not self.is_locally_solvable(b):
                        obs0.append(b)
                self._infinite_obstructions = obs0
        if finite:
            if read_cache and self._finite_obstructions is not None:
                obs1 = self._finite_obstructions
            else:
                candidates = []
                if self.determinant() != 0:
                    O = B.maximal_order()
                    for a in self.symmetric_matrix().list():
                        if a != 0:
                            for f in O.fractional_ideal(a).factor():
                                if f[1] < 0 and f[0] not in candidates:
                                    candidates.append(f[0])
                    for f in O.fractional_ideal(2 * self.determinant()).factor():
                        if f[1] > 0 and f[0] not in candidates:
                            candidates.append(f[0])
                for b in candidates:
                    if not self.is_locally_solvable(b):
                        obs1.append(b)
                self._infinite_obstructions = obs1
        obs = obs1 + obs0
        if finite and infinite:
            assert len(obs) % 2 == 0
        return obs
