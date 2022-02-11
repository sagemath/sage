# -*- coding: utf-8 -*-
r"""
Examples of simplicial sets.

These are accessible via ``simplicial_sets.Sphere(3)``,
``simplicial_sets.Torus()``, etc.  Type ``simplicial_sets.[TAB]`` to
see a complete list.

AUTHORS:

- John H. Palmieri (2016-07)
"""
#*****************************************************************************
#  Copyright (C) 2016 John H. Palmieri <palmieri at math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#
#*****************************************************************************

import re
import os
from pyparsing import OneOrMore, nestedExpr

from sage.env import SAGE_ENV
from sage.groups.abelian_gps.abelian_group import AbelianGroup
from sage.misc.cachefunc import cached_method, cached_function
from sage.misc.latex import latex
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from sage.structure.parent import Parent

from .delta_complex import delta_complexes
from .simplicial_set import AbstractSimplex, \
    SimplicialSet_arbitrary, SimplicialSet_finite

import sage.topology.simplicial_complex_catalog as simplicial_complexes

from sage.misc.lazy_import import lazy_import
lazy_import('sage.categories.simplicial_sets', 'SimplicialSets')

########################################################################
# The nerve of a finite monoid, used in sage.categories.finite_monoid.

class Nerve(SimplicialSet_arbitrary):
    def __init__(self, monoid):
        """
        The nerve of a multiplicative monoid.

        INPUT:

        - ``monoid`` -- a multiplicative monoid

        See
        :meth:`sage.categories.finite_monoids.FiniteMonoids.ParentMethods.nerve`
        for full documentation.

        EXAMPLES::

            sage: M = FiniteMonoids().example()
            sage: M
            An example of a finite multiplicative monoid: the integers modulo 12
            sage: X = M.nerve()
            sage: list(M)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
            sage: X.n_cells(0)
            [1]
            sage: X.n_cells(1)
            [0, 10, 11, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        category = SimplicialSets().Pointed()
        Parent.__init__(self, category=category)
        self.rename("Nerve of {}".format(str(monoid)))
        self.rename_latex("B{}".format(latex(monoid)))

        e = AbstractSimplex(0, name=str(monoid.one()),
                            latex_name=latex(monoid.one()))
        self._basepoint = e
        vertex = SimplicialSet_finite({e: None}, base_point=e)
        # self._n_skeleton: cache the highest dimensional skeleton
        # calculated so far for this simplicial set, along with its
        # dimension.
        self._n_skeleton = (0, vertex)
        self._monoid = monoid
        # self._simplex_data: a tuple whose elements are pairs (simplex, list
        # of monoid elements). Omit the base point.
        self._simplex_data = ()

    def __eq__(self, other):
        """
        Return ``True`` if ``self`` and ``other`` are equal.

        This checks that the underlying monoids and the underlying
        base points are the same. Because the base points will be
        different each time the nerve is constructed, different
        instances will not be equal.

        EXAMPLES::

            sage: C3 = groups.misc.MultiplicativeAbelian([3])
            sage: C3.nerve() == C3.nerve()
            False
            sage: BC3 = C3.nerve()
            sage: BC3 == BC3
            True
        """
        return (isinstance(other, Nerve)
                and self._monoid == other._monoid
                and self.base_point() == other.base_point())

    def __ne__(self, other):
        """
        Return the negation of `__eq__`.

        EXAMPLES::

            sage: C3 = groups.misc.MultiplicativeAbelian([3])
            sage: G3 = groups.permutation.Cyclic(3)
            sage: C3.nerve() != G3.nerve()
            True
            sage: C3.nerve() != C3.nerve()
            True
        """
        return not self == other

    @cached_method
    def __hash__(self):
        """
        The hash is formed from the monoid and the base point.

        EXAMPLES::

            sage: G3 = groups.permutation.Cyclic(3)
            sage: hash(G3.nerve()) # random
            17

        Different instances yield different base points, hence different hashes::

            sage: X = G3.nerve()
            sage: Y = G3.nerve()
            sage: X.base_point() != Y.base_point()
            True
            sage: hash(X) != hash(Y)
            True
        """
        return hash(self._monoid) ^ hash(self.base_point())

    def n_skeleton(self, n):
        """
        Return the `n`-skeleton of this simplicial set.

        That is, the simplicial set generated by all nondegenerate
        simplices of dimension at most `n`.

        INPUT:

        - ``n`` -- the dimension

        EXAMPLES::

            sage: K4 = groups.misc.MultiplicativeAbelian([2,2])
            sage: BK4 = simplicial_sets.ClassifyingSpace(K4)
            sage: BK4.n_skeleton(3)
            Simplicial set with 40 non-degenerate simplices
            sage: BK4.n_cells(1) == BK4.n_skeleton(3).n_cells(1)
            True
            sage: BK4.n_cells(3) == BK4.n_skeleton(1).n_cells(3)
            False
        """
        from .simplicial_set_constructions import SubSimplicialSet
        monoid = self._monoid
        one = monoid.one()
        # Build up chains of elements inductively, from dimension d-1
        # to dimension d. We start with the cached
        # self._n_skeleton. If only the 0-skeleton has been
        # constructed, we construct the 1-cells by hand.
        start, skel = self._n_skeleton
        if start == n:
            return skel
        elif start > n:
            return skel.n_skeleton(n)

        # There is a single vertex. Name it after the identity
        # element of the monoid.
        e = skel.n_cells(0)[0]
        # Build the dictionary simplices, to be used for
        # constructing the simplicial set.
        simplices = skel.face_data()

        # face_dict: dictionary of simplices: keys are
        # composites of monoid elements (as tuples), values are
        # the corresponding simplices.
        face_dict = dict(self._simplex_data)

        if start == 0:
            for g in monoid:
                if g != one:
                    x = AbstractSimplex(1, name=str(g), latex_name=latex(g))
                    simplices[x] = (e, e)
                    face_dict[(g,)] = x
            start = 1

        for d in range(start+1, n+1):
            for g in monoid:
                if g == one:
                    continue
                new_faces = {}
                for t in face_dict.keys():
                    if len(t) != d-1:
                        continue
                    # chain: chain of group elements to multiply,
                    # as a tuple.
                    chain = t + (g,)
                    # bdries: the face maps applied to chain, in a
                    # format suitable for passing to the DeltaComplex
                    # constructor.
                    x = AbstractSimplex(d,
                                        name=' * '.join(str(_) for _ in chain),
                                        latex_name = ' * '.join(latex(_) for _ in chain))
                    new_faces[chain] = x

                    # Compute faces of x.
                    faces = [face_dict[chain[1:]]]
                    for i in range(d-1):
                        product = chain[i] * chain[i+1]
                        if product == one:
                            # Degenerate.
                            if d == 2:
                                face = e.apply_degeneracies(i)
                            else:
                                face = (face_dict[chain[:i]
                                         + chain[i+2:]].apply_degeneracies(i))
                        else:
                            # Non-degenerate.
                            face = (face_dict[chain[:i]
                                              + (product,) + chain[i+2:]])
                        faces.append(face)
                    faces.append(face_dict[chain[:-1]])
                    simplices[x] = faces
                face_dict.update(new_faces)

        K = SubSimplicialSet(simplices, self)
        self._n_skeleton = (n, K)
        self._simplex_data = face_dict.items()
        return K


########################################################################
# Catalog of examples. These are accessed via simplicial_set_catalog.py.

def Sphere(n):
    r"""
    Return the `n`-sphere as a simplicial set.

    It is constructed with two non-degenerate simplices: a vertex
    `v_0` (which is the base point) and an `n`-simplex `\sigma_n`.

    INPUT:

    - ``n`` -- integer

    EXAMPLES::

        sage: S0 = simplicial_sets.Sphere(0)
        sage: S0
        S^0
        sage: S0.nondegenerate_simplices()
        [v_0, w_0]
        sage: S0.is_pointed()
        True
        sage: simplicial_sets.Sphere(4)
        S^4
        sage: latex(simplicial_sets.Sphere(4))
        S^{4}
        sage: simplicial_sets.Sphere(4).nondegenerate_simplices()
        [v_0, sigma_4]
    """
    v_0 = AbstractSimplex(0, name='v_0')
    if n == 0:
        w_0 = AbstractSimplex(0, name='w_0')
        return SimplicialSet_finite({v_0: None, w_0: None}, base_point=v_0,
                             name='S^0')
    degens = range(n-2, -1, -1)
    degen_v = v_0.apply_degeneracies(*degens)
    sigma = AbstractSimplex(n, name='sigma_{}'.format(n),
                            latex_name='\\sigma_{}'.format(n))
    return SimplicialSet_finite({sigma: [degen_v] * (n+1)}, base_point=v_0,
                                name='S^{}'.format(n),
                                latex_name='S^{{{}}}'.format(n))


def ClassifyingSpace(group):
    r"""
    Return the classifying space of ``group``, as a simplicial set.

    INPUT:

    - ``group`` -- a finite group or finite monoid

    See
    :meth:`sage.categories.finite_monoids.FiniteMonoids.ParentMethods.nerve`
    for more details and more examples.

    EXAMPLES::

        sage: C2 = groups.misc.MultiplicativeAbelian([2])
        sage: BC2 = simplicial_sets.ClassifyingSpace(C2)
        sage: H = BC2.homology(range(9), base_ring=GF(2))
        sage: [H[i].dimension() for i in range(9)]
        [0, 1, 1, 1, 1, 1, 1, 1, 1]

        sage: Klein4 = groups.misc.MultiplicativeAbelian([2, 2])
        sage: BK = simplicial_sets.ClassifyingSpace(Klein4)
        sage: BK
        Classifying space of Multiplicative Abelian group isomorphic to C2 x C2
        sage: BK.homology(range(5), base_ring=GF(2))  # long time (1 second)
        {0: Vector space of dimension 0 over Finite Field of size 2,
         1: Vector space of dimension 2 over Finite Field of size 2,
         2: Vector space of dimension 3 over Finite Field of size 2,
         3: Vector space of dimension 4 over Finite Field of size 2,
         4: Vector space of dimension 5 over Finite Field of size 2}
    """
    X = group.nerve()
    X.rename('Classifying space of {}'.format(group))
    return X


def RealProjectiveSpace(n):
    r"""
    Return real `n`-dimensional projective space, as a simplicial set.

    This is constructed as the `n`-skeleton of the nerve of the group
    of order 2, and therefore has a single non-degenerate simplex in
    each dimension up to `n`.

    EXAMPLES::

        sage: simplicial_sets.RealProjectiveSpace(7)
        RP^7
        sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
        sage: RP5.homology()
        {0: 0, 1: C2, 2: 0, 3: C2, 4: 0, 5: Z}
        sage: RP5
        RP^5
        sage: latex(RP5)
        RP^{5}

        sage: BC2 = simplicial_sets.RealProjectiveSpace(Infinity)
        sage: latex(BC2)
        RP^{\infty}
    """
    if n == Infinity:
        X = AbelianGroup([2]).nerve()
        X.rename('RP^oo')
        X.rename_latex('RP^{\\infty}')
    else:
        X = RealProjectiveSpace(Infinity).n_skeleton(n)
        X.rename('RP^{}'.format(n))
        X.rename_latex('RP^{{{}}}'.format(n))
    return X


def KleinBottle():
    r"""
    Return the Klein bottle as a simplicial set.

    This converts the `\Delta`-complex version to a simplicial set. It
    has one 0-simplex, three 1-simplices, and two 2-simplices.

    EXAMPLES::

        sage: K = simplicial_sets.KleinBottle()
        sage: K.f_vector()
        [1, 3, 2]
        sage: K.homology(reduced=False)
        {0: Z, 1: Z x C2, 2: 0}
        sage: K
        Klein bottle
    """
    temp = SimplicialSet_finite(delta_complexes.KleinBottle())
    pt = temp.n_cells(0)[0]
    return SimplicialSet_finite(temp.face_data(), base_point=pt,
                         name='Klein bottle')


def Torus():
    r"""
    Return the torus as a simplicial set.

    This computes the product of the circle with itself, where the
    circle is represented using a single 0-simplex and a single
    1-simplex. Thus it has one 0-simplex, three 1-simplices, and two
    2-simplices.

    EXAMPLES::

        sage: T = simplicial_sets.Torus()
        sage: T.f_vector()
        [1, 3, 2]
        sage: T.homology(reduced=False)
        {0: Z, 1: Z x Z, 2: Z}
    """
    S1 = Sphere(1)
    T = S1.product(S1)
    T.rename('Torus')
    return T


def Simplex(n):
    r"""
    Return the `n`-simplex as a simplicial set.

    EXAMPLES::

        sage: K = simplicial_sets.Simplex(2)
        sage: K
        2-simplex
        sage: latex(K)
        \Delta^{2}
        sage: K.n_cells(0)
        [(0,), (1,), (2,)]
        sage: K.n_cells(1)
        [(0, 1), (0, 2), (1, 2)]
        sage: K.n_cells(2)
        [(0, 1, 2)]
    """
    return SimplicialSet_finite(simplicial_complexes.Simplex(n),
                                name='{}-simplex'.format(n),
                                latex_name='\\Delta^{{{}}}'.format(n))


@cached_function
def Empty():
    """
    Return the empty simplicial set.

    This should return the same simplicial set each time it is called.

    EXAMPLES::

        sage: from sage.topology.simplicial_set_examples import Empty
        sage: E = Empty()
        sage: E
        Empty simplicial set
        sage: E.nondegenerate_simplices()
        []
        sage: E is Empty()
        True
    """
    return SimplicialSet_finite({}, name='Empty simplicial set')


@cached_function
def Point():
    """
    Return a single point called "*" as a simplicial set.

    This should return the same simplicial set each time it is called.

    EXAMPLES::

        sage: P = simplicial_sets.Point()
        sage: P.is_pointed()
        True
        sage: P.nondegenerate_simplices()
        [*]

        sage: Q = simplicial_sets.Point()
        sage: P is Q
        True
        sage: P == Q
        True
    """
    star = AbstractSimplex(0, name='*')
    return SimplicialSet_finite({star: None}, base_point=star,
                                name='Point',
                                latex_name='*')


def Horn(n, k):
    r"""
    Return the horn $\Lambda^n_k$.

    This is the subsimplicial set of the $n$-simplex obtained by
    removing its $k$-th face.

    EXAMPLES::

        sage: L = simplicial_sets.Horn(3, 0)
        sage: L
        (3, 0)-Horn
        sage: L.n_cells(3)
        []
        sage: L.n_cells(2)
        [(0, 1, 2), (0, 1, 3), (0, 2, 3)]

        sage: L20 = simplicial_sets.Horn(2, 0)
        sage: latex(L20)
        \Lambda^{2}_{0}
        sage: L20.inclusion_map()
        Simplicial set morphism:
          From: (2, 0)-Horn
          To:   2-simplex
          Defn: [(0,), (1,), (2,), (0, 1), (0, 2)] --> [(0,), (1,), (2,), (0, 1), (0, 2)]
    """
    K = Simplex(n)
    sigma = K.n_cells(n)[0]
    L = K.subsimplicial_set(K.faces(sigma)[:k] + K.faces(sigma)[k+1:])
    L.rename('({}, {})-Horn'.format(n, k))
    L.rename_latex('\\Lambda^{{{}}}_{{{}}}'.format(n, k))
    return L


def ComplexProjectiveSpace(n):
    r"""
    Return complex `n`-dimensional projective space, as a simplicial set.

    This is only defined when `n` is at most 4. It is constructed
    using the simplicial set decomposition provided by Kenzo, as
    described by Sergeraert [Ser2010]_

    EXAMPLES::

        sage: simplicial_sets.ComplexProjectiveSpace(2).homology(reduced=False)
        {0: Z, 1: 0, 2: Z, 3: 0, 4: Z}
        sage: CP3 = simplicial_sets.ComplexProjectiveSpace(3)
        sage: CP3
        CP^3
        sage: latex(CP3)
        CP^{3}
        sage: CP3.f_vector()
        [1, 0, 3, 10, 25, 30, 15]

        sage: K = CP3.suspension() # long time (1 second)
        sage: R = K.cohomology_ring(GF(2)) # long time
        sage: R.gens()        # long time
        (h^{0,0}, h^{3,0}, h^{5,0}, h^{7,0})
        sage: x = R.gens()[1] # long time
        sage: x.Sq(2)         # long time
        h^{5,0}

        sage: simplicial_sets.ComplexProjectiveSpace(4).f_vector()
        [1, 0, 4, 22, 97, 255, 390, 315, 105]

        sage: simplicial_sets.ComplexProjectiveSpace(5)
        Traceback (most recent call last):
        ...
        ValueError: complex projective spaces are only available in dimensions between 0 and 4
    """
    if n < 0 or n > 4:
        raise ValueError('complex projective spaces are only available in dimensions between 0 and 4')
    if n == 0:
        return Point()
    if n == 1:
        return Sphere(2)
    if n == 2:
        # v: Kenzo name <<GBar>>
        v = AbstractSimplex(0, name='v')
        # f_2_i: Kenzo name <<GBar<- (i)><- NIL>>> for i=1,2
        f2_1 = AbstractSimplex(2, name='rho_0')
        f2_2 = AbstractSimplex(2, name='rho_1')
        # f3_110: Kenzo name <<GBar<- (1 1)><0 NIL><- NIL>>>
        # f3_011: Kenzo name <<GBar<0 (1)><- (1)><- NIL>>>
        # f3_111: Kenzo name <<GBar<1 (1)><- (1)><- NIL>>>
        f3_110 = AbstractSimplex(3, name='sigma_0', latex_name='\\sigma_0')
        f3_011 = AbstractSimplex(3, name='sigma_1', latex_name='\\sigma_1')
        f3_111 = AbstractSimplex(3, name='sigma_2', latex_name='\\sigma_2')
        # f4_101101: Kenzo name <<GBar<1-0 (1)><1-0 NIL><- (1)><- NIL>>>
        # f4_201110: Kenzo name <<GBar<2-0 (1)><1 (1)><0 NIL><- NIL>>>
        # f4_211010: Kenzo name <<GBar<2-1 (1)><0 (1)><0 NIL><- NIL>>>
        f4_101101 = AbstractSimplex(4, name='tau_0', latex_name='\\tau_0')
        f4_201110 = AbstractSimplex(4, name='tau_1', latex_name='\\tau_1')
        f4_211010 = AbstractSimplex(4, name='tau_2', latex_name='\\tau_2')
        K = SimplicialSet_finite({f2_1: (v.apply_degeneracies(0),
                                  v.apply_degeneracies(0),
                                  v.apply_degeneracies(0)),
                           f2_2: (v.apply_degeneracies(0),
                                  v.apply_degeneracies(0),
                                  v.apply_degeneracies(0)),
                           f3_110: (f2_1, f2_2, f2_1, v.apply_degeneracies(1, 0)),
                           f3_011: (f2_1, f2_1, f2_1, f2_1),
                           f3_111: (v.apply_degeneracies(1, 0), f2_1, f2_2, f2_1),
                           f4_101101: (f2_1.apply_degeneracies(0),
                                       f2_1.apply_degeneracies(0),
                                       f3_011,
                                       f2_1.apply_degeneracies(2),
                                       f2_1.apply_degeneracies(2)),
                           f4_201110: (f2_1.apply_degeneracies(1),
                                       f3_111,
                                       f3_011,
                                       f3_110,
                                       f2_1.apply_degeneracies(1)),
                           f4_211010: (f2_1.apply_degeneracies(2),
                                       f3_111,
                                       f2_1.apply_degeneracies(1),
                                       f3_110,
                                       f2_1.apply_degeneracies(0))},
                           base_point=v, name='CP^2',
                           latex_name='CP^{2}')
        return K
    if n == 3:
        file = os.path.join(SAGE_ENV['SAGE_EXTCODE'], 'kenzo', 'CP3.txt')
        data = simplicial_data_from_kenzo_output(file)
        v = [_ for _ in data.keys() if _.dimension() == 0][0]
        K = SimplicialSet_finite(data, base_point=v, name='CP^3',
                                 latex_name='CP^{3}')
        return K
    if n == 4:
        file = os.path.join(SAGE_ENV['SAGE_EXTCODE'], 'kenzo', 'CP4.txt')
        data = simplicial_data_from_kenzo_output(file)
        v = [_ for _ in data.keys() if _.dimension() == 0][0]
        K = SimplicialSet_finite(data, base_point=v, name='CP^4',
                                 latex_name='CP^{4}')
        return K


def simplicial_data_from_kenzo_output(filename):
    """
    Return data to construct a simplicial set, given Kenzo output.

    INPUT:

    - ``filename`` -- name of file containing the output from Kenzo's
      :func:`show-structure` function

    OUTPUT: data to construct a simplicial set from the Kenzo output

    Several files with Kenzo output are in the directory
    :file:`SAGE_EXTCODE/kenzo/`.

    EXAMPLES::

        sage: from sage.topology.simplicial_set_examples import simplicial_data_from_kenzo_output
        sage: from sage.topology.simplicial_set import SimplicialSet
        sage: sphere = os.path.join(SAGE_ENV['SAGE_EXTCODE'], 'kenzo', 'S4.txt')
        sage: S4 = SimplicialSet(simplicial_data_from_kenzo_output(sphere))
        sage: S4.homology(reduced=False)
        {0: Z, 1: 0, 2: 0, 3: 0, 4: Z}
    """
    with open(filename, 'r') as f:
        data = f.read()
    dim = 0
    start = 0
    # simplex_data: data for constructing the simplicial set.
    simplex_data = {}
    # simplex_names: simplices indexed by their names
    simplex_names = {}
    dim_idx = data.find('Dimension = {}:'.format(dim), start)
    while dim_idx != -1:
        start = dim_idx + len('Dimension = {}:'.format(dim))
        new_dim_idx = data.find('Dimension = {}:'.format(dim+1), start)
        if new_dim_idx == -1:
            end = len(data)
        else:
            end = new_dim_idx
        if dim == 0:
            simplex_string = data[data.find('Vertices :') + len('Vertices :'):end]
            vertices = OneOrMore(nestedExpr()).parseString(simplex_string).asList()[0]
            for v in vertices:
                vertex = AbstractSimplex(0, name=v)
                simplex_data[vertex] = None
                simplex_names[v] = vertex
        else:
            simplex_string = data[start:end].strip()

            for s in [_.strip() for _ in simplex_string.split('Simplex : ')]:
                if s:
                    name, face_str = [_.strip() for _ in s.split('Faces : ')]
                    face_str = face_str.strip('()')
                    face_str = face_str.split('<AbSm ')
                    faces = []
                    for f in face_str[1:]:
                        # f has the form 'DEGENS NAME>', possibly with a trailing space.
                        # DEGENS is a hyphen-separated list, like
                        # '3-2-1-0' or '0' or '-'.
                        m = re.match('[-[0-9]+', f)
                        degen_str = m.group(0)
                        if degen_str.find('-') != -1:
                            if degen_str == '-':
                                degens = []
                            else:
                                degens = [Integer(_)
                                          for _ in degen_str.split('-')]
                        else:
                            degens = [Integer(degen_str)]

                        face_name = f[m.end(0):].strip()[:-1]
                        nondegen = simplex_names[face_name]
                        faces.append(nondegen.apply_degeneracies(*degens))

                    simplex = AbstractSimplex(dim, name=name)
                    simplex_data[simplex] = faces
                    simplex_names[name] = simplex
        dim += 1
        dim_idx = new_dim_idx
    return simplex_data


def HopfMap():
    r"""
    Return a simplicial model of the Hopf map `S^3 \to S^2`

    This is taken from Exemple II.1.19 in the thesis of Clemens Berger
    [Ber1991]_.

    The Hopf map is a fibration `S^3 \to S^2`. If it is viewed as
    attaching a 4-cell to the 2-sphere, the resulting adjunction space
    is 2-dimensional complex projective space. The resulting model is
    a bit larger than the one obtained from
    ``simplicial_sets.ComplexProjectiveSpace(2)``.

    EXAMPLES::

        sage: g = simplicial_sets.HopfMap()
        sage: g.domain()
        Simplicial set with 20 non-degenerate simplices
        sage: g.codomain()
        S^2

    Using the Hopf map to attach a cell::

        sage: X = g.mapping_cone()
        sage: CP2 = simplicial_sets.ComplexProjectiveSpace(2)
        sage: X.homology() == CP2.homology()
        True

        sage: X.f_vector()
        [1, 0, 5, 9, 6]
        sage: CP2.f_vector()
        [1, 0, 2, 3, 3]
    """
    # The 2-sphere and its simplices.
    S2 = Sphere(2)
    sigma = S2.n_cells(2)[0]
    s0_sigma = sigma.apply_degeneracies(0)
    s1_sigma = sigma.apply_degeneracies(1)
    s2_sigma = sigma.apply_degeneracies(2)
    # The 3-sphere and its simplices.
    w_0 = AbstractSimplex(0, name='w')
    w_1 = w_0.apply_degeneracies(0)
    w_2 = w_0.apply_degeneracies(0, 0)
    beta_11 = AbstractSimplex(1, name='beta_11', latex_name='\\beta_{11}')
    beta_22 = AbstractSimplex(1, name='beta_22', latex_name='\\beta_{22}')
    beta_23 = AbstractSimplex(1, name='beta_23', latex_name='\\beta_{23}')
    beta_44 = AbstractSimplex(1, name='beta_44', latex_name='\\beta_{44}')
    beta_1 = AbstractSimplex(2, name='beta_1', latex_name='\\beta_1')
    beta_2 = AbstractSimplex(2, name='beta_2', latex_name='\\beta_2')
    beta_3 = AbstractSimplex(2, name='beta_3', latex_name='\\beta_3')
    beta_4 = AbstractSimplex(2, name='beta_4', latex_name='\\beta_4')
    alpha_12 = AbstractSimplex(2, name='alpha_12', latex_name='\\alpha_{12}')
    alpha_23 = AbstractSimplex(2, name='alpha_23', latex_name='\\alpha_{23}')
    alpha_34 = AbstractSimplex(2, name='alpha_34', latex_name='\\alpha_{34}')
    alpha_45 = AbstractSimplex(2, name='alpha_45', latex_name='\\alpha_{45}')
    alpha_56 = AbstractSimplex(2, name='alpha_56', latex_name='\\alpha_{56}')
    alpha_1 = AbstractSimplex(3, name='alpha_1', latex_name='\\alpha_1')
    alpha_2 = AbstractSimplex(3, name='alpha_2', latex_name='\\alpha_2')
    alpha_3 = AbstractSimplex(3, name='alpha_3', latex_name='\\alpha_3')
    alpha_4 = AbstractSimplex(3, name='alpha_4', latex_name='\\alpha_4')
    alpha_5 = AbstractSimplex(3, name='alpha_5', latex_name='\\alpha_5')
    alpha_6 = AbstractSimplex(3, name='alpha_6', latex_name='\\alpha_6')
    S3 = SimplicialSet_finite({beta_11: (w_0, w_0), beta_22: (w_0, w_0),
                        beta_23: (w_0, w_0), beta_44: (w_0, w_0),
                        beta_1: (w_1, beta_11, w_1),
                        beta_2: (w_1, beta_22, beta_23),
                        beta_3: (w_1, beta_23, w_1),
                        beta_4: (w_1, beta_44, w_1),
                        alpha_12: (beta_11, beta_23, w_1),
                        alpha_23: (beta_11, beta_22, w_1),
                        alpha_34: (beta_11, beta_22, beta_44),
                        alpha_45: (w_1, beta_23, beta_44),
                        alpha_56: (w_1, beta_23, w_1),
                        alpha_1: (beta_1, beta_3, alpha_12, w_2),
                        alpha_2: (beta_11.apply_degeneracies(1), beta_2,
                                  alpha_23, alpha_12),
                        alpha_3: (beta_11.apply_degeneracies(0), alpha_34,
                                  alpha_23, beta_4),
                        alpha_4: (beta_1, beta_2, alpha_34, alpha_45),
                        alpha_5: (w_2, alpha_45, alpha_56, beta_4),
                        alpha_6: (w_2, beta_3, alpha_56, w_2)},
                       base_point=w_0)
    return S3.Hom(S2)({alpha_1:s0_sigma, alpha_2:s1_sigma,
                       alpha_3:s2_sigma, alpha_4:s0_sigma,
                       alpha_5:s2_sigma, alpha_6:s1_sigma})

