# sage.doctest: needs sage.rings.padics
r"""
Manin map

Represents maps from a set of right coset representatives to a
coefficient module.

This is a basic building block for implementing modular symbols, and
provides basic arithmetic and right action of matrices.

EXAMPLES::

    sage: E = EllipticCurve('11a')
    sage: phi = E.pollack_stevens_modular_symbol()
    sage: phi
    Modular symbol of level 11 with values in Sym^0 Q^2
    sage: phi.values()
    [-1/5, 1, 0]

    sage: from sage.modular.pollack_stevens.manin_map import ManinMap, M2Z
    sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
    sage: D = OverconvergentDistributions(0, 11, 10)
    sage: MR = ManinRelations(11)
    sage: data  = {M2Z([1,0,0,1]):D([1,2]), M2Z([0,-1,1,3]):D([3,5]), M2Z([-1,-1,3,2]):D([1,1])}
    sage: f = ManinMap(D, MR, data)
    sage: f(M2Z([1,0,0,1]))
    (1 + O(11^2), 2 + O(11))

    sage: S = Symk(0,QQ)
    sage: MR = ManinRelations(37)
    sage: data  = {M2Z([-2,-3,5,7]): S(0), M2Z([1,0,0,1]): S(0), M2Z([-1,-2,3,5]): S(0), M2Z([-1,-4,2,7]): S(1), M2Z([0,-1,1,4]): S(1), M2Z([-3,-1,7,2]): S(-1), M2Z([-2,-3,3,4]): S(0), M2Z([-4,-3,7,5]): S(0), M2Z([-1,-1,4,3]): S(0)}
    sage: f = ManinMap(S,MR,data)
    sage: f(M2Z([2,3,4,5]))
    1
"""

# ****************************************************************************
#       Copyright (C) 2012 Robert Pollack <rpollack@math.bu.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.continued_fraction import convergents
from .sigma0 import Sigma0
from .fund_domain import t00, t10, t01, t11, M2Z
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.integer_ring import ZZ
from sage.structure.element import coercion_model


def unimod_matrices_to_infty(r, s):
    r"""
    Return a list of matrices whose associated unimodular paths connect `0` to ``r/s``.

    INPUT:

    - ``r``, ``s`` -- rational numbers

    OUTPUT:

    - a list of matrices in `SL_2(\ZZ)`

    EXAMPLES::

        sage: v = sage.modular.pollack_stevens.manin_map.unimod_matrices_to_infty(19,23); v
        [
        [1 0]  [ 0  1]  [1 4]  [-4  5]  [ 5 19]
        [0 1], [-1  1], [1 5], [-5  6], [ 6 23]
        ]
        sage: [a.det() for a in v]
        [1, 1, 1, 1, 1]

        sage: sage.modular.pollack_stevens.manin_map.unimod_matrices_to_infty(11,25)
        [
        [1 0]  [ 0  1]  [1 3]  [-3  4]  [ 4 11]
        [0 1], [-1  2], [2 7], [-7  9], [ 9 25]
        ]


    ALGORITHM:

    This is Manin's continued fraction trick, which gives an expression
    `\{0,r/s\} = \{0,\infty\} + ... + \{a,b\} + ... + \{*,r/s\}`, where each `\{a,b\}` is
    the image of `\{0,\infty\}` under a matrix in `SL_2(\ZZ)`.
    """
    if s == 0:
        return []
    # the function contfrac_q in
    # https://github.com/williamstein/psage/blob/master/psage/modform/rational/modular_symbol_map.pyx
    # is very, very relevant to massively optimizing this.
    L = convergents(r / s)
    # Computes the continued fraction convergents of r/s
    v = [M2Z([1, L[0].numerator(), 0, L[0].denominator()])]
    # Initializes the list of matrices
    for j in range(len(L) - 1):
        a = L[j].numerator()
        c = L[j].denominator()
        b = L[j + 1].numerator()
        d = L[j + 1].denominator()
        v.append(M2Z([(-1) ** (j + 1) * a, b, (-1) ** (j + 1) * c, d]))
        # The matrix connecting two consecutive convergents is added on
    return v


def unimod_matrices_from_infty(r, s):
    r"""
    Return a list of matrices whose associated unimodular paths connect `\infty` to ``r/s``.

    INPUT:

    - ``r``, ``s`` -- rational numbers

    OUTPUT:

    - a list of `SL_2(\ZZ)` matrices

    EXAMPLES::

        sage: v = sage.modular.pollack_stevens.manin_map.unimod_matrices_from_infty(19,23); v
        [
        [ 0  1]  [-1  0]  [-4  1]  [-5 -4]  [-19   5]
        [-1  0], [-1 -1], [-5  1], [-6 -5], [-23   6]
        ]
        sage: [a.det() for a in v]
        [1, 1, 1, 1, 1]

        sage: sage.modular.pollack_stevens.manin_map.unimod_matrices_from_infty(11,25)
        [
        [ 0  1]  [-1  0]  [-3  1]  [-4 -3]  [-11   4]
        [-1  0], [-2 -1], [-7  2], [-9 -7], [-25   9]
        ]


    ALGORITHM:

    This is Manin's continued fraction trick, which gives an expression
    `\{\infty,r/s\} = \{\infty,0\} + ... + \{a,b\} + ... + \{*,r/s\}`, where each
    `\{a,b\}` is the image of `\{0,\infty\}` under a matrix in `SL_2(\ZZ)`.
    """
    if s != 0:
        L = convergents(r / s)
        # Computes the continued fraction convergents of r/s
        v = [M2Z([-L[0].numerator(), 1, -L[0].denominator(), 0])]
        # Initializes the list of matrices
        # the function contfrac_q in https://github.com/williamstein/psage/blob/master/psage/modform/rational/modular_symbol_map.pyx
        # is very, very relevant to massively optimizing this.
        for j in range(len(L) - 1):
            a = L[j].numerator()
            c = L[j].denominator()
            b = L[j + 1].numerator()
            d = L[j + 1].denominator()
            v.append(M2Z([-b, (-1) ** (j + 1) * a, -d, (-1) ** (j + 1) * c]))
            # The matrix connecting two consecutive convergents is added on
        return v
    else:
        return []


class ManinMap:
    r"""
    Map from a set of right coset representatives of `\Gamma_0(N)` in
    `SL_2(\ZZ)` to a coefficient module that satisfies the Manin
    relations.

    INPUT:

    - ``codomain`` -- coefficient module
    - ``manin_relations`` -- a :class:`sage.modular.pollack_stevens.fund_domain.ManinRelations`
      object
    - ``defining_data`` -- dictionary whose keys are a superset of
      ``manin_relations.gens()`` and a subset of ``manin_relations.reps()``,
      and whose values are in the codomain
    - ``check`` -- do numerous (slow) checks and transformations to
      ensure that the input data is perfect

    EXAMPLES::

        sage: from sage.modular.pollack_stevens.manin_map import M2Z, ManinMap
        sage: D = OverconvergentDistributions(0, 11, 10)
        sage: manin = sage.modular.pollack_stevens.fund_domain.ManinRelations(11)
        sage: data  = {M2Z([1,0,0,1]):D([1,2]), M2Z([0,-1,1,3]):D([3,5]), M2Z([-1,-1,3,2]):D([1,1])}
        sage: f = ManinMap(D, manin, data); f # indirect doctest
        Map from the set of right cosets of Gamma0(11) in SL_2(Z) to Space of 11-adic distributions with k=0 action and precision cap 10
        sage: f(M2Z([1,0,0,1]))
        (1 + O(11^2), 2 + O(11))
    """
    def __init__(self, codomain, manin_relations, defining_data, check=True):
        """
        INPUT:

        - ``codomain`` -- coefficient module
        - ``manin_relations`` -- a :class:`ManinRelations` object
        - ``defining_data`` -- dictionary whose keys are a superset of
          :meth:`manin_relations.gens()` and a subset of manin_relations.reps(),
          and whose values are in the codomain
        - ``check`` -- do numerous (slow) checks and transformations to
          ensure that the input data is perfect

        TESTS:

        Test that it fails gracefully on some bogus inputs::

            sage: from sage.modular.pollack_stevens.manin_map import  ManinMap
            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: rels = ManinRelations(37)
            sage: ManinMap(ZZ, rels, {})
            Traceback (most recent call last):
            ...
            ValueError: Codomain must have an action of Sigma0(N)
            sage: ManinMap(Symk(0), rels, [])
            Traceback (most recent call last):
            ...
            ValueError: length of defining data must be the same as number of Manin generators
        """
        self._codomain = codomain
        self._manin = manin_relations
        if check:
            if coercion_model.get_action(codomain, Sigma0(manin_relations._N)) is None:
                raise ValueError("Codomain must have an action of Sigma0(N)")
            self._dict = {}
            if isinstance(defining_data, (list, tuple)):
                if len(defining_data) != manin_relations.ngens():
                    raise ValueError("length of defining data must be the same as number of Manin generators")
                for i in range(len(defining_data)):
                    self._dict[manin_relations.gen(i)] = codomain(defining_data[i])
            elif isinstance(defining_data, dict):
                for g in manin_relations.gens():
                    self._dict[g] = codomain(defining_data[g])
            else:
                # constant function
                try:
                    c = codomain(defining_data)
                except TypeError:
                    raise TypeError("unrecognized type %s for defining_data" % type(defining_data))
                g = manin_relations.gens()
                self._dict = dict(zip(g, [c] * len(g)))
        else:
            self._dict = defining_data

    def extend_codomain(self, new_codomain, check=True):
        r"""
        Extend the codomain of ``self`` to ``new_codomain``. There must be a
        valid conversion operation from the old to the new codomain. This is
        most often used for extension of scalars from `\QQ` to `\QQ_p`.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.manin_map import ManinMap, M2Z
            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: S = Symk(0,QQ)
            sage: MR = ManinRelations(37)
            sage: data  = {M2Z([-2,-3,5,7]): S(0), M2Z([1,0,0,1]): S(0), M2Z([-1,-2,3,5]): S(0), M2Z([-1,-4,2,7]): S(1), M2Z([0,-1,1,4]): S(1), M2Z([-3,-1,7,2]): S(-1), M2Z([-2,-3,3,4]): S(0), M2Z([-4,-3,7,5]): S(0), M2Z([-1,-1,4,3]): S(0)}
            sage: m = ManinMap(S, MR, data); m
            Map from the set of right cosets of Gamma0(37) in SL_2(Z) to Sym^0 Q^2
            sage: m.extend_codomain(Symk(0, Qp(11)))
            Map from the set of right cosets of Gamma0(37) in SL_2(Z) to Sym^0 Q_11^2
        """
        new_dict = {}
        for g in self._manin.gens():
            new_dict[g] = new_codomain(self._dict[g])
        return ManinMap(new_codomain, self._manin, new_dict, check)

    def _compute_image_from_gens(self, B):
        r"""
        Compute the image of ``B`` under ``self``.

        INPUT:

        - ``B`` -- generator of Manin relations

        OUTPUT:

        An element in the codomain of ``self`` (e.g. a distribution), the image
        of ``B`` under ``self``.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.manin_map import M2Z, ManinMap
            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: D = OverconvergentDistributions(0, 11, 10)
            sage: MR = ManinRelations(11)
            sage: data  = {M2Z([1,0,0,1]):D([1,2]), M2Z([0,-1,1,3]):D([3,5]), M2Z([-1,-1,3,2]):D([1,1])}
            sage: f = ManinMap(D, MR, data)
            sage: f._compute_image_from_gens(MR.reps()[1])
            (10 + 10*11 + O(11^2), 8 + O(11))
        """
        L = self._manin.relations(B)
        # could raise KeyError if B is not a generator
        t = self._codomain(0)
        for c, A, g in L:
            g1 = self._dict[self._manin.reps(g)] * A
            t += g1 * c
        return t.normalize()

    def __getitem__(self, B):
        r"""

        Compute the image of ``B`` under ``self``.

        INPUT:

        - ``B`` -- coset representative of Manin relations

        OUTPUT:

        An element in the codomain of ``self`` (e.g. a distribution), the image
        of ``B`` under ``self``.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.manin_map import M2Z, ManinMap
            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: S = Symk(0,QQ)
            sage: MR = ManinRelations(37); MR.gens()
            [
            [1 0]  [ 0 -1]  [-1 -1]  [-1 -2]  [-2 -3]  [-3 -1]  [-1 -4]  [-4 -3]
            [0 1], [ 1  4], [ 4  3], [ 3  5], [ 5  7], [ 7  2], [ 2  7], [ 7  5],
            <BLANKLINE>
            [-2 -3]
            [ 3  4]
            ]

            sage: data  = {M2Z([-2,-3,5,7]): S(0), M2Z([1,0,0,1]): S(0), M2Z([-1,-2,3,5]): S(0), M2Z([-1,-4,2,7]): S(1), M2Z([0,-1,1,4]): S(1), M2Z([-3,-1,7,2]): S(-1), M2Z([-2,-3,3,4]): S(0), M2Z([-4,-3,7,5]): S(0), M2Z([-1,-1,4,3]): S(0)}
            sage: D = OverconvergentDistributions(2, 37, 40)
            sage: f = ManinMap(D, MR, data)
            sage: f.__getitem__(MR.gens()[1])
            1 + O(37)
            sage: f.__getitem__(MR.gens()[3])
            O(37^40)
            sage: f.__getitem__(MR.gens()[5])
            36 + O(37)
            sage: f[MR.gens()[5]]
            36 + O(37)
        """
        try:
            return self._dict[B]
        except KeyError:
            # To prevent memory overflow
            return self._compute_image_from_gens(B)
            # self._dict[B] = self._compute_image_from_gens(B)
            # return self._dict[B]

    def compute_full_data(self):
        r"""
        Compute the values of ``self`` on all coset reps from its values on our
        generating set.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.manin_map import M2Z, ManinMap
            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: S = Symk(0,QQ)
            sage: MR = ManinRelations(37); MR.gens()
            [
            [1 0]  [ 0 -1]  [-1 -1]  [-1 -2]  [-2 -3]  [-3 -1]  [-1 -4]  [-4 -3]
            [0 1], [ 1  4], [ 4  3], [ 3  5], [ 5  7], [ 7  2], [ 2  7], [ 7  5],
            <BLANKLINE>
            [-2 -3]
            [ 3  4]
            ]

            sage: data  = {M2Z([-2,-3,5,7]): S(0), M2Z([1,0,0,1]): S(0), M2Z([-1,-2,3,5]): S(0), M2Z([-1,-4,2,7]): S(1), M2Z([0,-1,1,4]): S(1), M2Z([-3,-1,7,2]): S(-1), M2Z([-2,-3,3,4]): S(0), M2Z([-4,-3,7,5]): S(0), M2Z([-1,-1,4,3]): S(0)}
            sage: f = ManinMap(S,MR,data)
            sage: len(f._dict)
            9
            sage: f.compute_full_data()
            sage: len(f._dict)
            38
        """
        for B in self._manin.reps():
            if B not in self._dict:
                self._dict[B] = self._compute_image_from_gens(B)

    def __add__(self, right):
        r"""
        Return sum ``self + right``, where ``self`` and ``right`` are
        assumed to have identical codomains and Manin relations.

        INPUT:

        - ``self``, ``right`` -- two Manin maps with the same codomain and
          Manin relations

        OUTPUT: the sum of ``self`` and ``right`` -- a Manin map

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.manin_map import M2Z, ManinMap
            sage: D = OverconvergentDistributions(0, 11, 10); D
            Space of 11-adic distributions with k=0 action and precision cap 10
            sage: manin = sage.modular.pollack_stevens.fund_domain.ManinRelations(11)
            sage: data  = {M2Z([1,0,0,1]):D([1,2]), M2Z([0,-1,1,3]):D([3,5]), M2Z([-1,-1,3,2]):D([1,1])}
            sage: f = ManinMap(D, manin, data); f
            Map from the set of right cosets of Gamma0(11) in SL_2(Z) to Space of 11-adic distributions with k=0 action and precision cap 10
            sage: f(M2Z([1,0,0,1]))
            (1 + O(11^2), 2 + O(11))
            sage: f+f # indirect doctest
            Map from the set of right cosets of Gamma0(11) in SL_2(Z) to Space of 11-adic distributions with k=0 action and precision cap 10
            sage: (f+f)(M2Z([1,0,0,1]))
            (2 + O(11^2), 4 + O(11))
        """
        D = {}
        sd = self._dict
        rd = right._dict
        for ky, val in sd.items():
            if ky in rd:
                D[ky] = val + rd[ky]
        return self.__class__(self._codomain, self._manin, D, check=False)

    def __sub__(self, right):
        """
        Return difference ``self`` - right, where ``self`` and ``right`` are
        assumed to have identical codomains and Manin relations.

        INPUT:

        - ``self``, ``right`` -- two Manin maps with the same codomain and
          Manin relations

        OUTPUT: the difference of ``self`` and ``right`` -- a Manin map

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.manin_map import M2Z, ManinMap
            sage: D = OverconvergentDistributions(0, 11, 10); D
            Space of 11-adic distributions with k=0 action and precision cap 10
            sage: manin = sage.modular.pollack_stevens.fund_domain.ManinRelations(11)
            sage: data  = {M2Z([1,0,0,1]):D([1,2]), M2Z([0,-1,1,3]):D([3,5]), M2Z([-1,-1,3,2]):D([1,1])}
            sage: f = ManinMap(D, manin, data); f
            Map from the set of right cosets of Gamma0(11) in SL_2(Z) to Space of 11-adic distributions with k=0 action and precision cap 10
            sage: f(M2Z([1,0,0,1]))
            (1 + O(11^2), 2 + O(11))
            sage: f-f
            Map from the set of right cosets of Gamma0(11) in SL_2(Z) to Space of 11-adic distributions with k=0 action and precision cap 10
            sage: (f-f)(M2Z([1,0,0,1]))
            (O(11^2), O(11))
        """
        D = {}
        sd = self._dict
        rd = right._dict
        for ky, val in sd.items():
            if ky in rd:
                D[ky] = val - rd[ky]
        return self.__class__(self._codomain, self._manin, D, check=False)

    def __mul__(self, right):
        """
        Return scalar multiplication ``self * right``, where ``right`` is in
        the base ring of the codomain.

        INPUT:

        - ``self`` -- a Manin map
        - ``right`` -- an element of the base ring of the codomain of self

        OUTPUT: the sum ``self`` and ``right`` -- a Manin map

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.manin_map import M2Z, ManinMap
            sage: D = OverconvergentDistributions(0, 11, 10)
            sage: manin = sage.modular.pollack_stevens.fund_domain.ManinRelations(11)
            sage: data  = {M2Z([1,0,0,1]):D([1,2]), M2Z([0,-1,1,3]):D([3,5]), M2Z([-1,-1,3,2]):D([1,1])}
            sage: f = ManinMap(D, manin, data)
            sage: f(M2Z([1,0,0,1]))
            (1 + O(11^2), 2 + O(11))
            sage: f*2
            Map from the set of right cosets of Gamma0(11) in SL_2(Z) to Space of 11-adic distributions with k=0 action and precision cap 10
            sage: (f*2)(M2Z([1,0,0,1]))
            (2 + O(11^2), 4 + O(11))
        """
        tp = Sigma0(self._manin.level())(MatrixSpace(ZZ, 2, 2)([1, 0, 0, 1]))
        if isinstance(right, type(tp)):
            return self._right_action(right)

        D = {}
        for ky, val in self._dict.items():
            D[ky] = val * right
        return self.__class__(self._codomain, self._manin, D, check=False)

    def __repr__(self):
        """
        Return string representation of ``self``.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.manin_map import M2Z, ManinMap
            sage: D = OverconvergentDistributions(0, 11, 10)
            sage: manin = sage.modular.pollack_stevens.fund_domain.ManinRelations(11)
            sage: data  = {M2Z([1,0,0,1]):D([1,2]), M2Z([0,-1,1,3]):D([3,5]), M2Z([-1,-1,3,2]):D([1,1])}
            sage: f = ManinMap(D, manin, data)
            sage: f.__repr__()
            'Map from the set of right cosets of Gamma0(11) in SL_2(Z) to Space of 11-adic distributions with k=0 action and precision cap 10'
        """
        return "Map from the set of right cosets of Gamma0(%s) in SL_2(Z) to %s" % (self._manin.level(), self._codomain)

    def _eval_sl2(self, A):
        r"""
        Return the value of ``self`` on the unimodular divisor corresponding to `A`.

        Note that `A` must be in `SL_2(Z)` for this to work.

        INPUT:

        - ``A`` -- an element of `SL_2(Z)`

        OUTPUT:

        The value of ``self`` on the divisor corresponding to `A` -- i.e. on
        the divisor `\{A(0)\} - \{A(\infty)\}`.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.manin_map import M2Z, ManinMap
            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: D = OverconvergentDistributions(0, 11, 10)
            sage: MR = ManinRelations(11)
            sage: data  = {M2Z([1,0,0,1]):D([1,2]), M2Z([0,-1,1,3]):D([3,5]), M2Z([-1,-1,3,2]):D([1,1])}
            sage: f = ManinMap(D, MR, data)
            sage: A = MR.reps()[1]
            sage: f._eval_sl2(A)
            (10 + 10*11 + O(11^2), 8 + O(11))
        """
        SN = Sigma0(self._manin._N)
        A = M2Z(A)
        B = self._manin.equivalent_rep(A)
        gaminv = SN(B * M2Z(A).adjugate())
        return (self[B] * gaminv).normalize()

    def __call__(self, A):
        """
        Evaluate ``self`` at A.

        INPUT:

        - ``A`` -- a `2 \times 2` matrix

        OUTPUT: the value of ``self`` on the divisor corresponding to ``A`` --
        an element of the codomain of self

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.manin_map import M2Z, ManinMap
            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: D = OverconvergentDistributions(0, 11, 10); D
            Space of 11-adic distributions with k=0 action and precision cap 10
            sage: manin = ManinRelations(11)
            sage: data  = {M2Z([1,0,0,1]):D([1,2]), M2Z([0,-1,1,3]):D([3,5]), M2Z([-1,-1,3,2]):D([1,1])}
            sage: f = ManinMap(D, manin, data); f
            Map from the set of right cosets of Gamma0(11) in SL_2(Z) to Space of 11-adic distributions with k=0 action and precision cap 10
            sage: f(M2Z([1,0,0,1]))
            (1 + O(11^2), 2 + O(11))

            sage: S = Symk(0,QQ)
            sage: MR = ManinRelations(37)
            sage: data  = {M2Z([-2,-3,5,7]): S(0), M2Z([1,0,0,1]): S(0), M2Z([-1,-2,3,5]): S(0), M2Z([-1,-4,2,7]): S(1), M2Z([0,-1,1,4]): S(1), M2Z([-3,-1,7,2]): S(-1), M2Z([-2,-3,3,4]): S(0), M2Z([-4,-3,7,5]): S(0), M2Z([-1,-1,4,3]): S(0)}
            sage: f = ManinMap(S,MR,data)
            sage: f(M2Z([2,3,4,5]))
            1
        """
        a = A[t00]
        b = A[t01]
        c = A[t10]
        d = A[t11]
        # v1: a list of unimodular matrices whose divisors add up to {b/d} - {infty}
        v1 = unimod_matrices_to_infty(b, d)
        # v2: a list of unimodular matrices whose divisors add up to {a/c} - {infty}
        v2 = unimod_matrices_to_infty(a, c)
        # ans: the value of self on A
        ans = self._codomain(0)
        # This loop computes self({b/d}-{infty}) by adding up the values of self on elements of v1
        for B in v1:
            ans = ans + self._eval_sl2(B)

        # This loops subtracts away the value self({a/c}-{infty}) from ans by subtracting away the values of self on elements of v2
        # and so in the end ans becomes self({b/d}-{a/c}) = self({A(0)} - {A(infty)}
        for B in v2:
            ans = ans - self._eval_sl2(B)
        return ans.normalize()

    def apply(self, f, codomain=None, to_moments=False):
        r"""
        Return Manin map given by `x \mapsto f(self(x))`, where `f` is
        anything that can be called with elements of the coefficient
        module.

        This might be used to normalize, reduce modulo a prime, change
        base ring, etc.

        INPUT:

        - ``f`` -- anything that can be called with elements of the coefficient
          module
        - ``codomain`` -- (default: ``None``) the codomain of the return map
        - ``to_moments`` -- boolean (default: ``False``); if ``True``, will
          apply ``f`` to each of the moments instead

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.manin_map import M2Z, ManinMap
            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: S = Symk(0,QQ)
            sage: MR = ManinRelations(37)
            sage: data  = {M2Z([-2,-3,5,7]): S(0), M2Z([1,0,0,1]): S(0), M2Z([-1,-2,3,5]): S(0), M2Z([-1,-4,2,7]): S(1), M2Z([0,-1,1,4]): S(1), M2Z([-3,-1,7,2]): S(-1), M2Z([-2,-3,3,4]): S(0), M2Z([-4,-3,7,5]): S(0), M2Z([-1,-1,4,3]): S(0)}
            sage: f = ManinMap(S,MR,data)
            sage: list(f.apply(lambda t:2*t))
            [0, 2, 0, 0, 0, -2, 2, 0, 0]
        """
        D = {}
        sd = self._dict
        if codomain is None:
            codomain = self._codomain
        for ky, val in sd.items():
            if to_moments:
                D[ky] = codomain([f(val.moment(a))
                                  for a in range(val.precision_absolute())])
            else:
                D[ky] = f(val)
        return self.__class__(codomain, self._manin, D, check=False)

    def __iter__(self):
        r"""
        Return iterator over the values of this map on the reduced
        representatives.

        This might be used to compute the valuation.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.manin_map import M2Z, ManinMap
            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: S = Symk(0,QQ)
            sage: MR = ManinRelations(37)
            sage: data  = {M2Z([-2,-3,5,7]): S(0), M2Z([1,0,0,1]): S(0), M2Z([-1,-2,3,5]): S(0), M2Z([-1,-4,2,7]): S(1), M2Z([0,-1,1,4]): S(1), M2Z([-3,-1,7,2]): S(-1), M2Z([-2,-3,3,4]): S(0), M2Z([-4,-3,7,5]): S(0), M2Z([-1,-1,4,3]): S(0)}
            sage: f = ManinMap(S,MR,data)
            sage: [a for a in f]
            [0, 1, 0, 0, 0, -1, 1, 0, 0]
        """
        for A in self._manin.gens():
            yield self._dict[A]

    def _right_action(self, gamma):
        r"""
        Return `self | \gamma`, where `\gamma` is a `2 \times 2` integer matrix.

        The action is defined by `(self | \gamma)(D) = self(\gamma D)|\gamma`

        For the action by a single element `\gamma` to be a modular symbol, `\gamma`
        must normalize `\Gamma_0(N)`.  However, this right action
        can also be used to define Hecke operators, in which case each
        individual `self | \gamma` is not a modular symbol on `\Gamma_0(N)`, but
        the sum over acting by the appropriate double coset representatives is.

        INPUT:

        - ``gamma`` -- `2 \times 2` integer matrix of nonzero determinant, with a
          well-defined action on the coefficient module

        OUTPUT: the image of ``self`` under the action of `\gamma` -- a Manin map

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.manin_map import ManinMap, M2Z, Sigma0
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_simple_modsym_space
            sage: S01 = Sigma0(1)
            sage: f = Newforms(7, 4)[0]
            sage: f.modular_symbols(1)
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(7) of weight 4 with sign 1 over Rational Field
            sage: phi = ps_modsym_from_simple_modsym_space(f.modular_symbols(1))._map
            sage: psi = phi._right_action(S01([2,3,4,5])); psi
            Map from the set of right cosets of Gamma0(7) in SL_2(Z) to Sym^2 Q^2

            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_simple_modsym_space
            sage: M = ModularSymbols(17,4,1).cuspidal_subspace()
            sage: A = M.decomposition()
            sage: f = ps_modsym_from_simple_modsym_space(A[0])._map
            sage: g = f._right_action(S01([1,2,0,1]))
            sage: g
            Map from the set of right cosets of Gamma0(17) in SL_2(Z) to Sym^2 Q^2
            sage: x = sage.modular.pollack_stevens.fund_domain.M2Z([2,3,1,0])
            sage: g(x)
            (17, -34, 69)
        """
        D = {}
        # we should eventually replace the for loop with a call to apply_many
        for ky in self._dict:
            D[ky] = self(gamma * ky) * gamma
        return self.__class__(self._codomain, self._manin, D, check=False)

    def normalize(self):
        r"""
        Normalize every value of ``self`` -- e.g., reduce each value's
        `j`-th moment modulo `p^{N-j}`.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.manin_map import M2Z, ManinMap
            sage: D = OverconvergentDistributions(0, 11, 10)
            sage: manin = sage.modular.pollack_stevens.fund_domain.ManinRelations(11)
            sage: data  = {M2Z([1,0,0,1]):D([1,2]), M2Z([0,-1,1,3]):D([3,5]), M2Z([-1,-1,3,2]):D([1,1])}
            sage: f = ManinMap(D, manin, data)
            sage: f._dict[M2Z([1,0,0,1])]
            (1 + O(11^2), 2 + O(11))
            sage: g = f.normalize()
            sage: g._dict[M2Z([1,0,0,1])]
            (1 + O(11^2), 2 + O(11))
        """
        sd = self._dict
        for val in sd.values():
            val.normalize()
        return self

    def reduce_precision(self, M):
        r"""
        Reduce the precision of all the values of the Manin map.

        INPUT:

        - ``M`` -- integer; the new precision

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.manin_map import M2Z, ManinMap
            sage: D = OverconvergentDistributions(0, 11, 10)
            sage: manin = sage.modular.pollack_stevens.fund_domain.ManinRelations(11)
            sage: data  = {M2Z([1,0,0,1]):D([1,2]), M2Z([0,-1,1,3]):D([3,5]), M2Z([-1,-1,3,2]):D([1,1])}
            sage: f = ManinMap(D, manin, data)
            sage: f._dict[M2Z([1,0,0,1])]
            (1 + O(11^2), 2 + O(11))
            sage: g = f.reduce_precision(1)
            sage: g._dict[M2Z([1,0,0,1])]
            1 + O(11^2)
        """
        D = {}
        for ky, val in self._dict.items():
            D[ky] = val.reduce_precision(M)
        return self.__class__(self._codomain, self._manin, D, check=False)

    def specialize(self, *args):
        r"""
        Specialize all the values of the Manin map to a new coefficient
        module. Assumes that the codomain has a ``specialize`` method, and
        passes all its arguments to that method.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.manin_map import M2Z, ManinMap
            sage: D = OverconvergentDistributions(0, 11, 10)
            sage: manin = sage.modular.pollack_stevens.fund_domain.ManinRelations(11)
            sage: data  = {M2Z([1,0,0,1]):D([1,2]), M2Z([0,-1,1,3]):D([3,5]), M2Z([-1,-1,3,2]):D([1,1])}
            sage: f = ManinMap(D, manin, data)
            sage: g = f.specialize()
            sage: g._codomain
            Sym^0 Z_11^2
        """
        D = {}
        for ky, val in self._dict.items():
            D[ky] = val.specialize(*args)
        return self.__class__(self._codomain.specialize(*args), self._manin,
                              D, check=False)

    def hecke(self, ell, algorithm='prep'):
        r"""
        Return the image of this Manin map under the Hecke operator `T_{\ell}`.

        INPUT:

        - ``ell`` -- a prime

        - ``algorithm`` -- string; either ``'prep'`` (default) or ``'naive'``

        OUTPUT: the image of this ManinMap under the Hecke operator `T_{\ell}`

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: phi.values()
            [-1/5, 1, 0]
            sage: phi.is_Tq_eigensymbol(7,7,10)
            True
            sage: phi.hecke(7).values()
            [2/5, -2, 0]
            sage: phi.Tq_eigenvalue(7,7,10)
            -2
        """
        self.compute_full_data()
        self.normalize()
        M = self._manin

        if algorithm == 'prep':
            # psi will denote self | T_ell
            psi = {}
            for g in M.gens():
                psi_g = sum((self[h] * A for h, A in M.prep_hecke_on_gen_list(ell, g)), self._codomain(0))
                psi_g.normalize()
                psi[g] = psi_g
            return self.__class__(self._codomain, self._manin,
                                  psi, check=False).normalize()
        elif algorithm == 'naive':
            S0N = Sigma0(self._manin.level())
            psi = self._right_action(S0N([1, 0, 0, ell]))
            for a in range(1, ell):
                psi += self._right_action(S0N([1, a, 0, ell]))
            if self._manin.level() % ell != 0:
                psi += self._right_action(S0N([ell, 0, 0, 1]))
            return psi.normalize()
        else:
            raise ValueError('Algorithm must be either "naive" or "prep"')

    def p_stabilize(self, p, alpha, V):
        r"""
        Return the `p`-stabilization of ``self`` to level `N*p` on which
        `U_p` acts by `\alpha`.

        INPUT:

        - ``p`` -- a prime

        - ``alpha`` -- a `U_p`-eigenvalue

        - ``V`` -- a space of modular symbols

        OUTPUT: the image of this ManinMap under the Hecke operator `T_{\ell}`

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: f = phi._map
            sage: V = phi.parent()
            sage: f.p_stabilize(5,1,V)
            Map from the set of right cosets of Gamma0(11) in SL_2(Z) to Sym^0 Q^2
        """
        manin = V.source()
        S0 = Sigma0(self._codomain._act._Np)
        pmat = S0([p, 0, 0, 1])
        D = {}
        scalar = 1 / alpha
        W = self._codomain.change_ring(scalar.parent())
        for g in map(M2Z, manin.gens()):
            # we use scale here so that we do not need to define a
            # construction functor in order to scale by something
            # outside the base ring.
            D[g] = W(self._eval_sl2(g) - (self(pmat * g) * pmat).scale(scalar))
        ans = self.__class__(W, manin, D, check=False)
        return ans
