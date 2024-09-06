"""
Artin Groups

Artin groups are implemented as a particular case of finitely presented
groups. For finite-type Artin groups, there is a specific left normal
form using the Garside structure associated to the lift the long element
of the corresponding Coxeter group.

AUTHORS:

- Travis Scrimshaw (2018-02-05): Initial version
"""

# ****************************************************************************
#       Copyright (C) 2018 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.groups.free_group import FreeGroup
from sage.groups.finitely_presented import FinitelyPresentedGroup, FinitelyPresentedGroupElement
from sage.combinat.root_system.coxeter_matrix import CoxeterMatrix
from sage.combinat.root_system.coxeter_group import CoxeterGroup
from sage.rings.infinity import Infinity
from sage.structure.richcmp import richcmp, rich_to_bool


class ArtinGroupElement(FinitelyPresentedGroupElement):
    """
    An element of an Artin group.

    It is a particular case of element of a finitely presented group.

    EXAMPLES::

        sage: # needs sage.rings.number_field
        sage: A.<s1,s2,s3> = ArtinGroup(['B',3])
        sage: A
        Artin group of type ['B', 3]
        sage: s1 * s2 / s3 / s2
        s1*s2*s3^-1*s2^-1
        sage: A((1, 2, -3, -2))
        s1*s2*s3^-1*s2^-1
    """
    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT: string; a valid LaTeX math command sequence

        TESTS::

            sage: A = ArtinGroup(['B', 3])                                              # needs sage.rings.number_field
            sage: b = A([1, 2, 3, -1, 2, -3])                                           # needs sage.rings.number_field
            sage: b._latex_()                                                           # needs sage.rings.number_field
            '\\sigma_{1}\\sigma_{2}\\sigma_{3}\\sigma_{1}^{-1}\\sigma_{2}\\sigma_{3}^{-1}'

            sage: B = BraidGroup(4)
            sage: b = B([1, 2, 3, -1, 2, -3])
            sage: b._latex_()
            '\\sigma_{1}\\sigma_{2}\\sigma_{3}\\sigma_{1}^{-1}\\sigma_{2}\\sigma_{3}^{-1}'
            sage: B.one()._latex_()
            '1'
        """
        word = self.Tietze()
        if not word:
            return '1'
        return ''.join(r"\sigma_{%s}^{-1}" % (-i) if i < 0 else r"\sigma_{%s}" % i
                       for i in word)

    def exponent_sum(self):
        """
        Return the exponent sum of ``self``.

        OUTPUT: integer

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: A = ArtinGroup(['E',6])
            sage: b = A([1, 4, -3, 2])
            sage: b.exponent_sum()
            2
            sage: b = A([])
            sage: b.exponent_sum()
            0

            sage: B = BraidGroup(5)
            sage: b = B([1, 4, -3, 2])
            sage: b.exponent_sum()
            2
            sage: b = B([])
            sage: b.exponent_sum()
            0
        """
        return sum(s.sign() for s in self.Tietze())

    def coxeter_group_element(self, W=None):
        """
        Return the corresponding Coxeter group element under the natural
        projection.

        INPUT:

        - ``W`` -- (default: ``self.parent().coxeter_group()``) the image
          Coxeter group

        OUTPUT: an element of the Coxeter group  ``W``

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: B.<s1,s2,s3> = ArtinGroup(['B',3])
            sage: b = s1 * s2 / s3 / s2
            sage: b1 = b.coxeter_group_element(); b1
            [ 1 -1  0]
            [ 2 -1  0]
            [ a -a  1]
            sage: b.coxeter_group_element().reduced_word()
            [1, 2, 3, 2]
            sage: A.<s1,s2,s3> = ArtinGroup(['A',3])
            sage: c = s1 * s2 *s3
            sage: c1 = c.coxeter_group_element(); c1
            [4, 1, 2, 3]
            sage: c1.reduced_word()
            [3, 2, 1]
            sage: c.coxeter_group_element(W=SymmetricGroup(4))
            (1,4,3,2)
            sage: A.<s1,s2,s3> = BraidGroup(4)
            sage: c = s1 * s2 * s3^-1
            sage: c0 = c.coxeter_group_element(); c0
            [4, 1, 2, 3]
            sage: c1 = c.coxeter_group_element(W=SymmetricGroup(4)); c1
            (1,4,3,2)

        From an element of the Coxeter group it is possible to recover
        the image by the standard section to the Artin group::

            sage: # needs sage.rings.number_field
            sage: B(b1)
            s1*s2*s3*s2
            sage: A(c0)
            s1*s2*s3
            sage: A(c0) == A(c1)
            True
        """
        if W is None:
            W = self.parent().coxeter_group()
        s = W.simple_reflections()
        In = W.index_set()
        return W.prod(s[In[abs(i)-1]] for i in self.Tietze())


class FiniteTypeArtinGroupElement(ArtinGroupElement):
    """
    An element of a finite-type Artin group.
    """
    def _richcmp_(self, other, op):
        """
        Compare ``self`` and ``other``.

        TESTS::

            sage: # needs sage.rings.number_field
            sage: A = ArtinGroup(['B',3])
            sage: x = A([1, 2, 1])
            sage: y = A([2, 1, 2])
            sage: x == y
            True
            sage: x < y^(-1)
            True
            sage: A([]) == A.one()
            True
            sage: x = A([2, 3, 2, 3])
            sage: y = A([3, 2, 3, 2])
            sage: x == y
            True
            sage: x < y^(-1)
            True
        """
        if self.Tietze() == other.Tietze():
            return rich_to_bool(op, 0)
        nfself = [i.Tietze() for i in self.left_normal_form()]
        nfother = [i.Tietze() for i in other.left_normal_form()]
        return richcmp(nfself, nfother, op)

    def __hash__(self):
        r"""
        Return a hash value for ``self``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: B.<s1,s2,s3> = ArtinGroup(['B',3])
            sage: hash(s1*s3) == hash(s3*s1)
            True
            sage: hash(s1*s2) == hash(s2*s1)
            False
            sage: hash(s1*s2*s1) == hash(s2*s1*s2)
            True
            sage: hash(s2*s3*s2) == hash(s3*s2*s3)
            False
            sage: hash(s2*s3*s2*s3) == hash(s3*s2*s3*s2)
            True
        """
        return hash(tuple(i.Tietze() for i in self.left_normal_form()))

    @cached_method
    def left_normal_form(self):
        r"""
        Return the left normal form of ``self``.

        OUTPUT:

        A tuple of simple generators in the left normal form. The first
        element is a power of `\Delta`, and the rest are elements of the
        natural section lift from the corresponding Coxeter group.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: A = ArtinGroup(['B',3])
            sage: A([1]).left_normal_form()
            (1, s1)
            sage: A([-1]).left_normal_form()
            (s1^-1*(s2^-1*s1^-1*s3^-1)^2*s2^-1*s3^-1, s3*(s2*s3*s1)^2*s2)
            sage: A([1, 2, 2, 1, 2]).left_normal_form()
            (1, s1*s2*s1, s2*s1)
            sage: A([3, 3, -2]).left_normal_form()
            (s1^-1*(s2^-1*s1^-1*s3^-1)^2*s2^-1*s3^-1,
             s3*s1*s2*s3*s2*s1, s3, s3*s2*s3)
            sage: A([1, 2, 3, -1, 2, -3]).left_normal_form()
            (s1^-1*(s2^-1*s1^-1*s3^-1)^2*s2^-1*s3^-1,
             (s3*s1*s2)^2*s1, s1*s2*s3*s2)
            sage: A([1,2,1,3,2,1,3,2,3,3,2,3,1,2,3,1,2,3,1,2]).left_normal_form()
            ((s3*(s2*s3*s1)^2*s2*s1)^2, s3*s2)

            sage: B = BraidGroup(4)
            sage: b = B([1, 2, 3, -1, 2, -3])
            sage: b.left_normal_form()
            (s0^-1*s1^-1*s0^-1*s2^-1*s1^-1*s0^-1, s0*s1*s2*s1*s0, s0*s2*s1)
            sage: c = B([1])
            sage: c.left_normal_form()
            (1, s0)
        """
        lnfp = self._left_normal_form_coxeter()
        P = self.parent()
        return tuple([P.delta() ** lnfp[0]] +
                     [P._standard_lift(w) for w in lnfp[1:]])

    def _left_normal_form_coxeter(self):
        r"""
        Return the left normal form of the element, in the `\Delta`
        exponent and Coxeter group element form.

        OUTPUT: tuple whose first element is the power of `\Delta`, and the
        rest are the Coxeter elements corresponding to the simple factors

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: A = ArtinGroup(['E',6])
            sage: A([2, -4, 2, 3, 1, 3, 2, 1, -2])._left_normal_form_coxeter()
            (
                [ 0  0  0  0  0 -1]  [ 0  0 -1  1  0  0]  [-1  0  1  0  0  0]
                [ 0  1  0 -1  0  0]  [ 0 -1  0  1  0  0]  [ 0 -1  0  1  0  0]
                [ 0  0  0  0 -1  0]  [-1  0  0  1  0  0]  [ 0  0  1  0  0  0]
                [ 0  1 -1  0 -1  0]  [-1 -1  0  1  1  0]  [ 0  0  0  1  0  0]
                [ 0  0 -1  0  0  0]  [ 0  0  0  0  1  0]  [ 0  0  0  0  1  0]
            -1, [-1  0  0  0  0  0], [ 0  0  0  0  0  1], [ 0  0  0  0  0  1]
            )
            sage: A = ArtinGroup(['F',4])
            sage: A([2, 3, -4, 2, 3, -2, 1, -2, 3, 4, 1, -2])._left_normal_form_coxeter()
            (
                [-1  0  0  0]  [ -1   1   0   0]  [-1  0  0  0]
                [-1 -1  a -a]  [ -2   3  -a   0]  [-1  1 -a  0]
                [-a  0  1 -2]  [ -a 2*a  -1  -1]  [ 0  0 -1  0]
            -3, [-a  0  1 -1], [ -a   a   0  -1], [ 0  0  0 -1],
            <BLANKLINE>
            [  -1    1    0   -a]  [  -1    0    a   -a]
            [   0    1    0 -2*a]  [  -1    1    a -2*a]
            [   0    a   -1   -2]  [   0    0    2   -3]
            [   0    a   -1   -1], [   0    0    1   -1]
            )
        """
        delta = 0
        Delta = self.parent().coxeter_group().long_element()
        sr = self.parent().coxeter_group().simple_reflections()
        tz = self.Tietze()
        if tz == ():
            return (0,)
        form = []
        for i in tz:
            if i > 0:
                form.append(sr[i])
            else:
                delta += 1
                form = [Delta * a * Delta for a in form]
                form.append(Delta * sr[-i])
        i = j = 0
        while j < len(form):
            while i < len(form) - j - 1:
                e = form[i].descents(side='right')
                s = form[i + 1].descents(side='left')
                S = set(s).difference(set(e))
                while S:
                    a = list(S)[0]
                    form[i] = form[i] * sr[a]
                    form[i + 1] = sr[a] * form[i+1]
                    e = form[i].descents(side='right')
                    s = form[i + 1].descents(side='left')
                    S = set(s).difference(set(e))
                if form[i+1].length() == 0:
                    form.pop(i+1)
                    i = 0
                else:
                    i += 1
            j += 1
            i = 0
        form = [elt for elt in form if elt.length()]
        while form and form[0] == Delta:
            form.pop(0)
            delta -= 1
        return tuple([-delta] + form)


class ArtinGroup(FinitelyPresentedGroup):
    r"""
    An Artin group.

    Fix an index set `I`. Let `M = (m_{ij})_{i,j \in I}` be a
    :class:`Coxeter matrix
    <sage.combinat.root_system.coxeter_matrix.CoxeterMatrix>`.
    An *Artin group* is a group `A_M` that has a presentation
    given by generators `\{ s_i \mid i \in I \}` and relations

    .. MATH::

        \underbrace{s_i s_j s_i \cdots}_{m_{ij}}
        = \underbrace{s_j s_i s_j \cdots}_{\text{$m_{ji}$ factors}}

    for all `i,j \in I` with the usual convention that `m_{ij} = \infty`
    implies no relation between `s_i` and `s_j`. There is a natural
    corresponding Coxeter group `W_M` by imposing the additional
    relations `s_i^2 = 1` for all `i \in I`. Furthermore, there is
    a natural section of `W_M` by sending a reduced word
    `s_{i_1} \cdots s_{i_{\ell}} \mapsto s_{i_1} \cdots s_{i_{\ell}}`.

    Artin groups `A_M` are classified based on the Coxeter data:

    - `A_M` is of *finite type* or *spherical* if `W_M` is finite;
    - `A_M` is of *affine type* if `W_M` is of affine type;
    - `A_M` is of *large type* if `m_{ij} \geq 4` for all `i,j \in I`;
    - `A_M` is of *extra-large type* if `m_{ij} \geq 5` for all `i,j \in I`;
    - `A_M` is *right-angled* if `m_{ij} \in \{2,\infty\}` for all `i,j \in I`.

    Artin groups are conjectured to have many nice properties:

    - Artin groups are torsion free.
    - Finite type Artin groups have `Z(A_M) = \ZZ` and infinite type
      Artin groups have trivial center.
    - Artin groups have solvable word problems.
    - `H_{W_M} / W_M` is a `K(A_M, 1)`-space, where `H_W` is the
      hyperplane complement of the Coxeter group `W` acting on `\CC^n`.

    These conjectures are known when the Artin group is finite type and a
    number of other cases. See, e.g., [GP2012]_ and references therein.

    INPUT:

    - ``coxeter_data`` -- data defining a Coxeter matrix

    - ``names`` -- string or list/tuple/iterable of strings
      (default: ``'s'``); the generator names or name prefix

    EXAMPLES::

        sage: A.<a,b,c> = ArtinGroup(['B',3]);  A                                       # needs sage.rings.number_field
        Artin group of type ['B', 3]
        sage: ArtinGroup(['B',3])                                                       # needs sage.rings.number_field
        Artin group of type ['B', 3]

    The input must always include the Coxeter data, but the ``names``
    can either be a string representing the prefix of the names or
    the explicit names of the generators. Otherwise the default prefix
    of ``'s'`` is used::

        sage: ArtinGroup(['B',2]).generators()                                          # needs sage.rings.number_field
        (s1, s2)
        sage: ArtinGroup(['B',2], 'g').generators()                                     # needs sage.rings.number_field
        (g1, g2)
        sage: ArtinGroup(['B',2], 'x,y').generators()                                   # needs sage.rings.number_field
        (x, y)

    REFERENCES:

    - :wikipedia:`Artin_group`

    .. SEEALSO::

        :class:`~sage.groups.raag.RightAngledArtinGroup`
    """
    @staticmethod
    def __classcall_private__(cls, coxeter_data, names=None):
        """
        Normalize input to ensure a unique representation.

        TESTS::

            sage: # needs sage.rings.number_field
            sage: A1 = ArtinGroup(['B',3])
            sage: A2 = ArtinGroup(['B',3], 's')
            sage: A3 = ArtinGroup(['B',3],  ['s1','s2','s3'])
            sage: A1 is A2 and A2 is A3
            True

            sage: # needs sage.rings.number_field
            sage: A1 = ArtinGroup(['B',2], 'a,b')
            sage: A2 = ArtinGroup([[1,4],[4,1]], 'a,b')
            sage: A3.<a,b> = ArtinGroup('B2')
            sage: A1 is A2 and A2 is A3
            True

            sage: ArtinGroup(['A',3]) is BraidGroup(4, 's1,s2,s3')                      # needs sage.rings.number_field
            True

            sage: G = graphs.PathGraph(3)
            sage: CM = CoxeterMatrix([[1,-1,2],[-1,1,-1],[2,-1,1]], index_set=G.vertices(sort=True))
            sage: A = groups.misc.Artin(CM)                                             # needs sage.rings.number_field
            sage: Ap = groups.misc.RightAngledArtin(G, 's')                             # needs sage.rings.number_field
            sage: A is Ap                                                               # needs sage.rings.number_field
            True
        """
        coxeter_data = CoxeterMatrix(coxeter_data)
        if names is None:
            names = 's'
        if isinstance(names, str):
            if ',' in names:
                names = [x.strip() for x in names.split(',')]
            else:
                names = [names + str(i) for i in coxeter_data.index_set()]
        names = tuple(names)
        if len(names) != coxeter_data.rank():
            raise ValueError("the number of generators must match"
                             " the rank of the Coxeter type")
        if all(m == Infinity for m in coxeter_data.coxeter_graph().edge_labels()):
            from sage.groups.raag import RightAngledArtinGroup
            return RightAngledArtinGroup(coxeter_data.coxeter_graph(), names)
        if not coxeter_data.is_finite():
            raise NotImplementedError
        if coxeter_data.coxeter_type().cartan_type().type() == 'A':
            from sage.groups.braid import BraidGroup
            return BraidGroup(coxeter_data.rank()+1, names)
        return FiniteTypeArtinGroup(coxeter_data, names)

    def __init__(self, coxeter_matrix, names):
        """
        Initialize ``self``.

        TESTS::

            sage: # needs sage.rings.number_field
            sage: A = ArtinGroup(['D',4])
            sage: TestSuite(A).run()
            sage: A = ArtinGroup(['B',3], ['x','y','z'])
            sage: TestSuite(A).run()
        """
        self._coxeter_group = CoxeterGroup(coxeter_matrix)
        free_group = FreeGroup(names)
        rels = []
        # Generate the relations based on the Coxeter graph
        I = coxeter_matrix.index_set()
        for ii, i in enumerate(I):
            for j in I[ii + 1:]:
                m = coxeter_matrix[i, j]
                if m == Infinity:  # no relation
                    continue
                elt = [i, j] * m
                for ind in range(m, 2*m):
                    elt[ind] = -elt[ind]
                rels.append(free_group(elt))
        FinitelyPresentedGroup.__init__(self, free_group, tuple(rels))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: ArtinGroup(['B',3])                                                   # needs sage.rings.number_field
            Artin group of type ['B', 3]
            sage: ArtinGroup(['D',4], 'g')                                              # needs sage.rings.number_field
            Artin group of type ['D', 4]
        """
        try:
            data = self.coxeter_type().cartan_type()
            return "Artin group of type {}".format(data)
        except AttributeError:
            pass
        return "Artin group with Coxeter matrix:\n{}".format(self.coxeter_matrix())

    def cardinality(self):
        """
        Return the number of elements of ``self``.

        OUTPUT: infinity

        EXAMPLES::

            sage: Gamma = graphs.CycleGraph(5)
            sage: G = RightAngledArtinGroup(Gamma)
            sage: G.cardinality()
            +Infinity

            sage: A = ArtinGroup(['A',1])                                               # needs sage.rings.number_field
            sage: A.cardinality()                                                       # needs sage.rings.number_field
            +Infinity
        """
        from sage.rings.infinity import Infinity
        return Infinity

    order = cardinality

    def as_permutation_group(self):
        """
        Return an isomorphic permutation group.

        This raises a :exc:`ValueError` error since Artin groups are
        infinite and have no corresponding permutation group.

        EXAMPLES::

            sage: Gamma = graphs.CycleGraph(5)
            sage: G = RightAngledArtinGroup(Gamma)
            sage: G.as_permutation_group()
            Traceback (most recent call last):
            ...
            ValueError: the group is infinite

            sage: A = ArtinGroup(['D',4], 'g')                                          # needs sage.rings.number_field
            sage: A.as_permutation_group()                                              # needs sage.rings.number_field
            Traceback (most recent call last):
            ...
            ValueError: the group is infinite
        """
        raise ValueError("the group is infinite")

    def coxeter_type(self):
        """
        Return the Coxeter type of ``self``.

        EXAMPLES::

            sage: A = ArtinGroup(['D',4])                                               # needs sage.rings.number_field
            sage: A.coxeter_type()                                                      # needs sage.rings.number_field
            Coxeter type of ['D', 4]
        """
        return self._coxeter_group.coxeter_type()

    def coxeter_matrix(self):
        """
        Return the Coxeter matrix of ``self``.

        EXAMPLES::

            sage: A = ArtinGroup(['B',3])                                               # needs sage.rings.number_field
            sage: A.coxeter_matrix()                                                    # needs sage.rings.number_field
            [1 3 2]
            [3 1 4]
            [2 4 1]
        """
        return self._coxeter_group.coxeter_matrix()

    def coxeter_group(self):
        """
        Return the Coxeter group of ``self``.

        EXAMPLES::

            sage: A = ArtinGroup(['D',4])                                               # needs sage.rings.number_field
            sage: A.coxeter_group()                                                     # needs sage.rings.number_field
            Finite Coxeter group over Integer Ring with Coxeter matrix:
            [1 3 2 2]
            [3 1 3 3]
            [2 3 1 2]
            [2 3 2 1]
        """
        return self._coxeter_group

    def index_set(self):
        """
        Return the index set of ``self``.

        OUTPUT: tuple

        EXAMPLES::

            sage: A = ArtinGroup(['E',7])                                               # needs sage.rings.number_field
            sage: A.index_set()                                                         # needs sage.rings.number_field
            (1, 2, 3, 4, 5, 6, 7)
        """
        return self._coxeter_group.index_set()

    def _element_constructor_(self, x):
        """
        TESTS::

            sage: A = ArtinGroup(['B',3])                                               # needs sage.rings.number_field
            sage: A([2,1,-2,3,3,3,1])                                                   # needs sage.rings.number_field
            s2*s1*s2^-1*s3^3*s1
        """
        if x in self._coxeter_group:
            return self._standard_lift(x)
        return self.element_class(self, x)

    @cached_method
    def an_element(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: A = ArtinGroup(['B',2])                                               # needs sage.rings.number_field
            sage: A.an_element()                                                        # needs sage.rings.number_field
            s1
        """
        return self.gen(0)

    def some_elements(self):
        """
        Return a list of some elements of ``self``.

        EXAMPLES::

            sage: A = ArtinGroup(['B',3])                                               # needs sage.rings.number_field
            sage: A.some_elements()                                                     # needs sage.rings.number_field
            [s1, s1*s2*s3, (s1*s2*s3)^3]
        """
        rank = self.coxeter_matrix().rank()
        elements_list = [self.gen(0)]
        elements_list.append(self.prod(self.gens()))
        elements_list.append(elements_list[-1] ** min(rank, 3))
        return elements_list

    def _standard_lift_Tietze(self, w):
        """
        Return a Tietze word representing the Coxeter element ``w``
        under the natural section.

        INPUT:

        - ``w`` -- an element of the Coxeter group of ``self``

        EXAMPLES::

            sage: A = ArtinGroup(['B',3])                                               # needs sage.rings.number_field
            sage: A._standard_lift_Tietze(A.coxeter_group().long_element())             # needs sage.rings.number_field
            [3, 2, 3, 1, 2, 3, 1, 2, 1]
        """
        return w.reduced_word()

    @cached_method
    def _standard_lift(self, w):
        """
        Return the element of ``self`` that corresponds to the given
        Coxeter element ``w`` under the natural section.

        INPUT:

        - ``w`` -- an element of the Coxeter group of ``self``

        EXAMPLES::

            sage: A = ArtinGroup(['B',3])                                               # needs sage.rings.number_field
            sage: A._standard_lift(A.coxeter_group().long_element())                    # needs sage.rings.number_field
            s3*(s2*s3*s1)^2*s2*s1

            sage: B = BraidGroup(5)
            sage: P = Permutation([5, 3, 1, 2, 4])
            sage: B._standard_lift(P)
            s0*s1*s2*s3*s0*s1
        """
        return self(self._standard_lift_Tietze(w))

    Element = ArtinGroupElement


class FiniteTypeArtinGroup(ArtinGroup):
    r"""
    A finite-type Artin group.

    An Artin group is *finite-type* or *spherical* if the corresponding
    Coxeter group is finite. Finite type Artin groups are known to be
    torsion free, have a Garside structure given by `\Delta` (see
    :meth:`delta`) and have a center generated by `\Delta`.

    .. SEEALSO::

        :class:`ArtinGroup`

    EXAMPLES::

        sage: ArtinGroup(['E',7])                                                       # needs sage.rings.number_field
        Artin group of type ['E', 7]

    Since the word problem for finite-type Artin groups is solvable, their
    Cayley graph can be locally obtained as follows (see :issue:`16059`)::

        sage: def ball(group, radius):
        ....:     ret = set()
        ....:     ret.add(group.one())
        ....:     for length in range(1, radius):
        ....:         for w in Words(alphabet=group.gens(), length=length):
        ....:              ret.add(prod(w))
        ....:     return ret
        sage: A = ArtinGroup(['B',3])                                                   # needs sage.rings.number_field
        sage: GA = A.cayley_graph(elements=ball(A, 4), generators=A.gens()); GA         # needs sage.rings.number_field
        Digraph on 32 vertices

    Since the Artin group has nontrivial relations, this graph contains less
    vertices than the one associated to the free group (which is a tree)::

        sage: F = FreeGroup(3)
        sage: GF = F.cayley_graph(elements=ball(F, 4), generators=F.gens()); GF         # needs sage.combinat
        Digraph on 40 vertices
    """
    def delta(self):
        r"""
        Return the `\Delta` element of ``self``.

        EXAMPLES::

            sage: A = ArtinGroup(['B',3])                                               # needs sage.rings.number_field
            sage: A.delta()                                                             # needs sage.rings.number_field
            s3*(s2*s3*s1)^2*s2*s1

            sage: A = ArtinGroup(['G',2])                                               # needs sage.rings.number_field
            sage: A.delta()                                                             # needs sage.rings.number_field
            (s2*s1)^3

            sage: B = BraidGroup(5)
            sage: B.delta()
            s0*s1*s2*s3*s0*s1*s2*s0*s1*s0
        """
        return self._standard_lift(self._coxeter_group.long_element())

    Element = FiniteTypeArtinGroupElement
