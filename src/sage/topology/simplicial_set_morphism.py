r"""
Morphisms and homsets for simplicial sets

.. NOTE::

    Morphisms with infinite domain are not implemented in general:
    only constant maps and identity maps are currently implemented.

AUTHORS:

- John H. Palmieri (2016-07)

This module implements morphisms and homsets of simplicial sets.
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

import itertools

from sage.categories.homset import Hom, Homset
from sage.categories.morphism import Morphism
from sage.categories.simplicial_sets import SimplicialSets
from sage.matrix.constructor import matrix, zero_matrix
from sage.misc.latex import latex
from sage.rings.integer_ring import ZZ

from sage.homology.chain_complex_morphism import ChainComplexMorphism
from sage.homology.homology_morphism import InducedHomologyMorphism
from .simplicial_set import SimplicialSet_arbitrary

class SimplicialSetHomset(Homset):
    r"""
    A set of morphisms between simplicial sets.

    Once a homset has been constructed in Sage, typically via
    ``Hom(X,Y)`` or ``X.Hom(Y)``, one can use it to construct a
    morphism `f` by specifying a dictionary, the keys of which are the
    nondegenerate simplices in the domain, and the value corresponding
    to `\sigma` is the simplex `f(\sigma)` in the codomain.

    EXAMPLES::

        sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
        sage: v = AbstractSimplex(0, name='v')
        sage: w = AbstractSimplex(0, name='w')
        sage: e = AbstractSimplex(1, name='e')
        sage: f = AbstractSimplex(1, name='f')
        sage: X = SimplicialSet({e: (v, w), f: (w, v)})
        sage: Y = SimplicialSet({e: (v, v)})

    Define the homset::

        sage: H = Hom(X, Y)

    Now define a morphism by specifying a dictionary::

        sage: H({v: v, w: v, e: e, f: e})
        Simplicial set morphism:
          From: Simplicial set with 4 non-degenerate simplices
          To:   Simplicial set with 2 non-degenerate simplices
          Defn: [v, w, e, f] --> [v, v, e, e]
    """
    def __call__(self, f, check=True):
        r"""
        INPUT:

        - ``f`` -- a dictionary with keys the simplices of the domain
          and values simplices of the codomain

        - ``check`` -- optional, default ``True``. Pass this to the
          morphism constructor.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: v0 = S1.n_cells(0)[0]
            sage: e = S1.n_cells(1)[0]
            sage: f = {v0: v0, e: v0.apply_degeneracies(0)} # constant map
            sage: Hom(S1, S1)(f)
            Simplicial set endomorphism of S^1
              Defn: Constant map at v_0
        """
        return SimplicialSetMorphism(f, self.domain(), self.codomain(), check=check)

    def diagonal_morphism(self):
        r"""
        Return the diagonal morphism in `\operatorname{Hom}(S, S \times S)`.

        EXAMPLES::

            sage: RP2 = simplicial_sets.RealProjectiveSpace(2)
            sage: Hom(RP2, RP2.product(RP2)).diagonal_morphism()
            Simplicial set morphism:
              From: RP^2
              To:   RP^2 x RP^2
              Defn: [1, f, f * f] --> [(1, 1), (f, f), (f * f, f * f)]
        """
        domain = self.domain()
        codomain = self.codomain()
        if not hasattr(codomain, 'factors'):
            raise ValueError('diagonal morphism is only defined for Hom(X, XxX)')
        factors = codomain.factors()
        if len(factors) != 2 or factors[0] != domain or factors[1] != domain:
            raise ValueError('diagonal morphism is only defined for Hom(X, XxX)')
        f = {}
        for i in range(domain.dimension()+1):
            for s in domain.n_cells(i):
                f[s] = dict(codomain._translation)[((s, ()), (s, ()))]
        return self(f)

    def identity(self):
        r"""
        Return the identity morphism in `\operatorname{Hom}(S, S)`.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: Hom(S1, S1).identity()
            Simplicial set endomorphism of S^1
              Defn: Identity map
            sage: T = simplicial_sets.Torus()
            sage: Hom(S1, T).identity()
            Traceback (most recent call last):
            ...
            TypeError: identity map is only defined for endomorphism sets
        """
        return SimplicialSetMorphism(domain=self.domain(),
                                     codomain=self.codomain(),
                                     identity=True)

    def constant_map(self, point=None):
        r"""
        Return the constant map in this homset.

        INPUT:

        - ``point`` -- optional, default ``None``. If specified, it
          must be a 0-simplex in the codomain, and it will be the
          target of the constant map.

        If ``point`` is specified, it is the target of the constant
        map. Otherwise, if the codomain is pointed, the target is its
        base point. If the codomain is not pointed and ``point`` is
        not specified, raise an error.

        EXAMPLES::

            sage: S3 = simplicial_sets.Sphere(3)
            sage: T = simplicial_sets.Torus()
            sage: T.n_cells(0)[0].rename('w')
            sage: Hom(S3,T).constant_map()
            Simplicial set morphism:
              From: S^3
              To:   Torus
              Defn: Constant map at w

            sage: S0 = simplicial_sets.Sphere(0)
            sage: v, w = S0.n_cells(0)
            sage: Hom(S3, S0).constant_map(v)
            Simplicial set morphism:
              From: S^3
              To:   S^0
              Defn: Constant map at v_0
            sage: Hom(S3, S0).constant_map(w)
            Simplicial set morphism:
              From: S^3
              To:   S^0
              Defn: Constant map at w_0

        This constant map is not pointed, since it doesn't send the
        base point of `S^3` to the base point of `S^0`::

            sage: Hom(S3, S0).constant_map(w).is_pointed()
            False

        TESTS::

            sage: S0 = S0.unset_base_point()
            sage: Hom(S3, S0).constant_map()
            Traceback (most recent call last):
            ...
            ValueError: codomain is not pointed, so specify a target for the constant map
        """
        codomain = self.codomain()
        if point is None:
            if codomain.is_pointed():
                point = codomain.base_point()
            else:
                raise ValueError('codomain is not pointed, so specify a '
                                 'target for the constant map')
        return SimplicialSetMorphism(domain=self.domain(),
                                     codomain=self.codomain(),
                                     constant=point)

    def an_element(self):
        """
        Return an element of this homset: a constant map.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: S2 = simplicial_sets.Sphere(2)
            sage: Hom(S2, S1).an_element()
            Simplicial set morphism:
              From: S^2
              To:   S^1
              Defn: Constant map at v_0

            sage: K = simplicial_sets.Simplex(3)
            sage: L = simplicial_sets.Simplex(4)
            sage: d = {K.n_cells(3)[0]: L.n_cells(0)[0].apply_degeneracies(2, 1, 0)}
            sage: Hom(K,L)(d) == Hom(K,L).an_element()
            True
        """
        codomain = self.codomain()
        if codomain.is_pointed():
            target = codomain.base_point()
        else:
            target = codomain.n_cells(0)[0]
        return self.constant_map(target)

    def __iter__(self):
        """
        Iterate through all morphisms in this homset.

        This is very slow: it tries all possible targets for the
        maximal nondegenerate simplices and yields those which are
        valid morphisms of simplicial sets. ("Maximal" means
        nondegenerate simplices which are not the faces of other
        nondegenerate simplices.) So if either the domain or the
        codomain has many simplices, the number of possibilities may
        be quite large.

        This is only implemented when the domain is finite.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: T = simplicial_sets.Torus()
            sage: H = Hom(S1, T)
            sage: list(H)
            [Simplicial set morphism:
               From: S^1
               To:   Torus
               Defn: [v_0, sigma_1] --> [(v_0, v_0), (s_0 v_0, sigma_1)],
             Simplicial set morphism:
               From: S^1
               To:   Torus
               Defn: [v_0, sigma_1] --> [(v_0, v_0), (sigma_1, s_0 v_0)],
             Simplicial set morphism:
               From: S^1
               To:   Torus
               Defn: [v_0, sigma_1] --> [(v_0, v_0), (sigma_1, sigma_1)],
             Simplicial set morphism:
               From: S^1
               To:   Torus
               Defn: Constant map at (v_0, v_0)]
            sage: [f.induced_homology_morphism().to_matrix() for f in H]
            [
            [ 1| 0]  [1|0]  [1|0]  [1|0]
            [--+--]  [-+-]  [-+-]  [-+-]
            [ 0|-1]  [0|1]  [0|0]  [0|0]
            [ 0| 1]  [0|0]  [0|1]  [0|0]
            [--+--]  [-+-]  [-+-]  [-+-]
            [ 0| 0], [0|0], [0|0], [0|0]
            ]
        """
        if not self.domain().is_finite():
            raise NotImplementedError('domain must be finite to iterate '
                                      'through all morphisms')
        codomain = self.codomain()
        facets = self.domain()._facets_()
        dims = [f.dimension() for f in facets]
        # Record all of the n-simplices in the codomain once for each
        # relevant dimension.
        all_n_simplices = {d: codomain.all_n_simplices(d) for d in set(dims)}
        for target in itertools.product(*[all_n_simplices[d] for d in dims]):
            try:
                yield self({sigma: tau for (sigma, tau) in zip(facets, target)})
            except ValueError:
                # Not a valid morphism.
                pass

    def _latex_(self):
        r"""
        LaTeX representation

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: T = simplicial_sets.Torus()
            sage: H = Hom(S1, T)
            sage: latex(H)
            \operatorname{Map} (S^{1}, S^{1} \times S^{1})
        """
        return '\\operatorname{{Map}} ({}, {})'.format(latex(self.domain()), latex(self.codomain()))


class SimplicialSetMorphism(Morphism):
    def __init__(self, data=None, domain=None, codomain=None,
                 constant=None, identity=False, check=True):
        r"""
        Return a morphism of simplicial sets.

        INPUT:

        - ``data`` -- optional. Dictionary defining the map.
        - ``domain`` -- simplicial set
        - ``codomain`` -- simplicial set
        - ``constant`` -- optional: if not ``None``, then this should
          be a vertex in the codomain, in which case return the
          constant map with this vertex as the target.
        - ``identity`` -- optional: if ``True``, return the identity
          morphism.
        - ``check`` -- optional, default ``True``. If ``True``, check
          that this is actually a morphism: it commutes with the face
          maps.

        So to define a map, you must specify ``domain`` and
        ``codomain``. If the map is constant, specify the target (a
        vertex in the codomain) as ``constant``. If the map is the
        identity map, specify ``identity=True``. Otherwise, pass a
        dictionary, ``data``.  The keys of the dictionary are the
        nondegenerate simplices of the domain, the corresponding
        values are simplices in the codomain.

        In fact, the keys in ``data`` do not need to include all of
        the nondegenerate simplices, only those which are not faces of
        other nondegenerate simplices: if `\sigma` is a face of
        `\tau`, then the image of `\sigma` need not be specified.

        EXAMPLES::

            sage: from sage.topology.simplicial_set_morphism import SimplicialSetMorphism
            sage: K = simplicial_sets.Simplex(1)
            sage: S1 = simplicial_sets.Sphere(1)
            sage: v0 = K.n_cells(0)[0]
            sage: v1 = K.n_cells(0)[1]
            sage: e01 = K.n_cells(1)[0]
            sage: w = S1.n_cells(0)[0]
            sage: sigma = S1.n_cells(1)[0]

            sage: f = {v0: w, v1: w, e01: sigma}
            sage: SimplicialSetMorphism(f, K, S1)
            Simplicial set morphism:
              From: 1-simplex
              To:   S^1
              Defn: [(0,), (1,), (0, 1)] --> [v_0, v_0, sigma_1]

        The same map can be defined as follows::

            sage: H = Hom(K, S1)
            sage: H(f)
            Simplicial set morphism:
              From: 1-simplex
              To:   S^1
              Defn: [(0,), (1,), (0, 1)] --> [v_0, v_0, sigma_1]

        Also, this map can be defined by specifying where the
        1-simplex goes; the vertices then go where they have to, to
        satisfy the condition `d_i \circ f = f \circ d_i`::

            sage: H = Hom(K, S1)
            sage: H({e01: sigma})
            Simplicial set morphism:
              From: 1-simplex
              To:   S^1
              Defn: [(0,), (1,), (0, 1)] --> [v_0, v_0, sigma_1]

        A constant map::

            sage: g = {e01: w.apply_degeneracies(0)}
            sage: SimplicialSetMorphism(g, K, S1)
            Simplicial set morphism:
              From: 1-simplex
              To:   S^1
              Defn: Constant map at v_0

        The same constant map::

            sage: SimplicialSetMorphism(domain=K, codomain=S1, constant=w)
            Simplicial set morphism:
              From: 1-simplex
              To:   S^1
              Defn: Constant map at v_0

        An identity map::

            sage: SimplicialSetMorphism(domain=K, codomain=K, identity=True)
            Simplicial set endomorphism of 1-simplex
              Defn: Identity map

        Defining a map by specifying it on only some of the simplices
        in the domain::

            sage: S5 = simplicial_sets.Sphere(5)
            sage: s = S5.n_cells(5)[0]
            sage: one = S5.Hom(S5)({s: s})
            sage: one
            Simplicial set endomorphism of S^5
              Defn: Identity map

        TESTS:

        A non-map::

            sage: h = {w: v0, sigma: e01}
            sage: SimplicialSetMorphism(h, S1, K)
            Traceback (most recent call last):
            ...
            ValueError: the dictionary does not define a map of simplicial sets

        Another non-map::

            sage: h = {w: v0, v0: w, sigma: e01}
            sage: SimplicialSetMorphism(h, S1, K)
            Traceback (most recent call last):
            ...
            ValueError: at least one simplex in the defining dictionary is not in the domain

        A non-identity map::

            sage: SimplicialSetMorphism(domain=K, codomain=S1, identity=True)
            Traceback (most recent call last):
            ...
            TypeError: identity map is only defined for endomorphism sets

        An improperly partially defined map::

            sage: h = {w: v0}
            sage: SimplicialSetMorphism(h, S1, K)
            Traceback (most recent call last):
            ...
            ValueError: the image of at least one simplex in the domain is not defined
        """
        self._is_identity = False
        if not domain.is_finite():
            if identity:
                if codomain is None:
                    codomain = domain
                elif domain is not codomain:
                    raise TypeError("identity map is only defined for endomorphism sets")
                self._is_identity = True
                Morphism.__init__(self, Hom(domain, codomain, SimplicialSets()))
                return
            if constant is not None:
                # If self._constant is set, it should be a vertex in
                # the codomain, the target of the constant map.
                self._constant = constant
                Morphism.__init__(self, Hom(domain, codomain, SimplicialSets()))
                return
            raise NotImplementedError('morphisms with infinite domain '
                                      'are not implemented in general')
        else:
            if identity:
                self._is_identity = True
                check = False
                if domain is not codomain:
                    raise TypeError("identity map is only defined for endomorphism sets")
                data = {}
                for i in range(domain.dimension() + 1):
                    for s in domain.n_cells(i):
                        data[s] = s
            if constant is not None:
                self._constant = constant
                check = False
                data = {sigma: constant.apply_degeneracies(*range(sigma.dimension()-1,-1,-1))
                        for sigma in domain.nondegenerate_simplices()}

        if (not isinstance(domain, SimplicialSet_arbitrary)
            or not isinstance(codomain, SimplicialSet_arbitrary)):
            raise TypeError('the domain and codomain must be simplicial sets')
        if any(x.nondegenerate() not in
               domain.nondegenerate_simplices() for x in data.keys()):
            raise ValueError('at least one simplex in the defining '
                             'dictionary is not in the domain')
        # Remove degenerate simplices from the domain specification.
        d = {sigma:data[sigma] for sigma in data if sigma.is_nondegenerate()}
        # For each simplex in d.keys(), add its faces, and the faces
        # of its faces, etc., to d.
        for simplex in list(d):
            faces = domain.faces(simplex)
            add = []
            if faces:
                for (i, sigma) in enumerate(faces):
                    nondegen = sigma.nondegenerate()
                    if nondegen not in d:
                        add.append((sigma, i, simplex))
            while add:
                (sigma, i, tau) = add.pop()
                # sigma is the ith face of tau.
                face_f = codomain.face(d[tau], i)
                degens = sigma.degeneracies()
                x = face_f
                for j in degens:
                    x = codomain.face(x, j)
                d[sigma.nondegenerate()] = x
                faces = domain.faces(sigma.nondegenerate())
                if faces:
                    for (i,rho) in enumerate(faces):
                        nondegen = rho.nondegenerate()
                        if nondegen not in d:
                            add.append((rho,i,sigma))
        # Now check that the proposed map commutes with the face
        # maps. (The degeneracy maps should work automatically.)
        if check:
            for simplex in d:
                # Compare d[d_i (simplex)] to d_i d[simplex]. Since
                # d_i(simplex) may be degenerate, we have to be careful
                # when applying f to it. We can skip vertices and start
                # with 1-simplices.
                bad = False
                for i in range(simplex.dimension()+1):
                    face_f = codomain.face(d[simplex], i)
                    face = domain.face(simplex, i)
                    if face is None:
                        f_face = None
                    elif face.is_nondegenerate():
                        f_face = d[face]
                    else:
                        nondegen = face.nondegenerate()
                        f_face = d[nondegen].apply_degeneracies(*face.degeneracies())
                    if face_f != f_face:
                        bad = True
                        break
                if bad:
                    raise ValueError('the dictionary does not define a map of simplicial sets')
        if any(x not in d.keys() for x in domain.nondegenerate_simplices()):
            raise ValueError('the image of at least one simplex in '
                             'the domain is not defined')
        self._dictionary = d
        Morphism.__init__(self, Hom(domain, codomain, SimplicialSets()))

    def __eq__(self, other):
        """
        Two morphisms are equal iff their domains are the same, their
        codomains are the same, and their defining dictionaries are
        the same.

        EXAMPLES::

            sage: S = simplicial_sets.Sphere(1)
            sage: T = simplicial_sets.Torus()
            sage: T_c = T.constant_map() * T.base_point_map()
            sage: S_c = S.constant_map() * S.base_point_map()
            sage: T_c == S_c
            True
            sage: T.constant_map() == S.constant_map()
            False
            sage: K = simplicial_sets.Sphere(1)
            sage: K.constant_map() == S.constant_map()
            False

            sage: Point = simplicial_sets.Point()
            sage: f = Point._map_from_empty_set()
            sage: Empty = f.domain()
            sage: g = Empty.constant_map()
            sage: f == g
            True
        """
        if self.domain().is_finite() and other.domain().is_finite():
            return (self.domain() == other.domain()
                    and self.codomain() == other.codomain()
                    and self._dictionary == other._dictionary)
        else:
            return False

    def __ne__(self, other):
        """
        The negation of ``__eq__``.

        EXAMPLES::

            sage: S0 = simplicial_sets.Sphere(0)
            sage: v,w = S0.n_cells(0)
            sage: H = Hom(S0, S0)
            sage: H({v:v, w:w}) != H({v:w, w:v})
            True
            sage: H({v:v, w:w}) != H({w:w, v:v})
            False
        """
        return not self == other

    def __call__(self, x):
        """
        INPUT: a simplex of the domain.

        Return its image under this morphism.

        EXAMPLES::

            sage: K = simplicial_sets.Simplex(1)
            sage: S1 = simplicial_sets.Sphere(1)
            sage: v0 = K.n_cells(0)[0]
            sage: v1 = K.n_cells(0)[1]
            sage: e01 = K.n_cells(1)[0]
            sage: w = S1.n_cells(0)[0]
            sage: sigma = S1.n_cells(1)[0]
            sage: d = {v0: w, v1: w, e01: sigma}
            sage: f = Hom(K, S1)(d)
            sage: f(e01) # indirect doctest
            sigma_1

            sage: one = Hom(S1, S1).identity()
            sage: e = S1.n_cells(1)[0]
            sage: one(e) == e
            True

            sage: B = AbelianGroup([2]).nerve()
            sage: c = B.constant_map()
            sage: c(B.n_cells(2)[0])
            s_1 s_0 *
        """
        if x not in self.domain():
            raise ValueError('element is not a simplex in the domain')
        if self.is_constant():
            target = self._constant
            return target.apply_degeneracies(*range(x.dimension()-1, -1, -1))
        if self._is_identity:
            return x
        return self._dictionary[x.nondegenerate()].apply_degeneracies(*x.degeneracies())

    def _composition_(self, right, homset):
        """
        Return the composition of two morphisms.

        INPUT:

        - ``self``, ``right`` -- maps
        - ``homset`` -- a homset

        ASSUMPTION:

        The codomain of ``right`` is contained in the domain of
        ``self``.  This assumption should be verified by the
        ``Map.__mul__`` method in ``categories/map.pyx``, so we don't
        need to check it here.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: f = S1.Hom(S1).identity()
            sage: f * f # indirect doctest
            Simplicial set endomorphism of S^1
              Defn: Identity map
            sage: T = S1.product(S1)
            sage: K = T.factor(0, as_subset=True)
            sage: g = S1.Hom(T)({S1.n_cells(0)[0]:K.n_cells(0)[0], S1.n_cells(1)[0]:K.n_cells(1)[0]})
            sage: g
            Simplicial set morphism:
              From: S^1
              To:   S^1 x S^1
              Defn: [v_0, sigma_1] --> [(v_0, v_0), (sigma_1, s_0 v_0)]
            sage: (g*f).image()
            Simplicial set with 2 non-degenerate simplices
            sage: f.image().homology()
            {0: 0, 1: Z}
        """
        if self.is_identity():
            return right
        if right.is_identity():
            return self
        d = {}
        for sigma in right._dictionary:
            d[sigma] = self(right(sigma))
        return homset(d)

    def image(self):
        """
        Return the image of this morphism as a subsimplicial set of the
        codomain.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: T = S1.product(S1)
            sage: K = T.factor(0, as_subset=True)
            sage: f = S1.Hom(T)({S1.n_cells(0)[0]:K.n_cells(0)[0], S1.n_cells(1)[0]:K.n_cells(1)[0]})
            sage: f
            Simplicial set morphism:
              From: S^1
              To:   S^1 x S^1
              Defn: [v_0, sigma_1] --> [(v_0, v_0), (sigma_1, s_0 v_0)]
            sage: f.image()
            Simplicial set with 2 non-degenerate simplices
            sage: f.image().homology()
            {0: 0, 1: Z}

            sage: B = simplicial_sets.ClassifyingSpace(groups.misc.MultiplicativeAbelian([2]))
            sage: B.constant_map().image()
            Point
            sage: Hom(B,B).identity().image() == B
            True
        """
        if self._is_identity:
            return self.codomain()
        if self.is_constant():
            return self.codomain().subsimplicial_set([self._constant])
        simplices = self._dictionary.values()
        if set(simplices) == set(self.codomain().nondegenerate_simplices()):
            return self.codomain()
        return self.codomain().subsimplicial_set(simplices)

    def is_identity(self):
        """
        Return ``True`` if this morphism is an identity map.

        EXAMPLES::

            sage: K = simplicial_sets.Simplex(1)
            sage: v0 = K.n_cells(0)[0]
            sage: v1 = K.n_cells(0)[1]
            sage: e01 = K.n_cells(1)[0]
            sage: L = simplicial_sets.Simplex(2).n_skeleton(1)
            sage: w0 = L.n_cells(0)[0]
            sage: w1 = L.n_cells(0)[1]
            sage: w2 = L.n_cells(0)[2]
            sage: f01 = L.n_cells(1)[0]
            sage: f02 = L.n_cells(1)[1]
            sage: f12 = L.n_cells(1)[2]

            sage: d = {v0:w0, v1:w1, e01:f01}
            sage: f = K.Hom(L)(d)
            sage: f.is_identity()
            False
            sage: d = {w0:v0, w1:v1, w2:v1, f01:e01, f02:e01, f12: v1.apply_degeneracies(0,)}
            sage: g = L.Hom(K)(d)
            sage: (g*f).is_identity()
            True
            sage: (f*g).is_identity()
            False
            sage: (f*g).induced_homology_morphism().to_matrix(1)
            [0]

            sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
            sage: RP5.n_skeleton(2).inclusion_map().is_identity()
            False
            sage: RP5.n_skeleton(5).inclusion_map().is_identity()
            True

            sage: B = simplicial_sets.ClassifyingSpace(groups.misc.MultiplicativeAbelian([2]))
            sage: Hom(B,B).identity().is_identity()
            True
            sage: Hom(B,B).constant_map().is_identity()
            False
        """
        ans = (self._is_identity or
                (self.domain() == self.codomain()
                 and self.domain().is_finite()
                 and all(a == b for a,b in self._dictionary.items())))
        self._is_identity = ans
        return ans

    def is_surjective(self):
        """
        Return ``True`` if this map is surjective.

        EXAMPLES::

            sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
            sage: RP2 = RP5.n_skeleton(2)
            sage: RP2.inclusion_map().is_surjective()
            False

            sage: RP5_2 = RP5.quotient(RP2)
            sage: RP5_2.quotient_map().is_surjective()
            True

            sage: K = RP5_2.pullback(RP5_2.quotient_map(), RP5_2.base_point_map())
            sage: f = K.universal_property(RP2.inclusion_map(), RP2.constant_map())
            sage: f.is_surjective()
            True
        """
        return self._is_identity or self.image() == self.codomain()

    def is_injective(self):
        """
        Return ``True`` if this map is injective.

        EXAMPLES::

            sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
            sage: RP2 = RP5.n_skeleton(2)
            sage: RP2.inclusion_map().is_injective()
            True

            sage: RP5_2 = RP5.quotient(RP2)
            sage: RP5_2.quotient_map().is_injective()
            False

            sage: K = RP5_2.pullback(RP5_2.quotient_map(), RP5_2.base_point_map())
            sage: f = K.universal_property(RP2.inclusion_map(), RP2.constant_map())
            sage: f.is_injective()
            True
        """
        if self._is_identity:
            return True
        domain = self.domain()
        for n in range(domain.dimension()+1):
            input = domain.n_cells(n)
            output = set([self(sigma) for sigma in input if self(sigma).is_nondegenerate()])
            if len(input) > len(output):
                return False
        return True

    def is_bijective(self):
        """
        Return ``True`` if this map is bijective.

        EXAMPLES::

            sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
            sage: RP2 = RP5.n_skeleton(2)
            sage: RP2.inclusion_map().is_bijective()
            False

            sage: RP5_2 = RP5.quotient(RP2)
            sage: RP5_2.quotient_map().is_bijective()
            False

            sage: K = RP5_2.pullback(RP5_2.quotient_map(), RP5_2.base_point_map())
            sage: f = K.universal_property(RP2.inclusion_map(), RP2.constant_map())
            sage: f.is_bijective()
            True
        """
        return self.is_injective() and self.is_surjective()

    def is_pointed(self):
        """
        Return ``True`` if this is a pointed map.

        That is, return ``True`` if the domain and codomain are
        pointed and this morphism preserves the base point.

        EXAMPLES::

            sage: S0 = simplicial_sets.Sphere(0)
            sage: f = Hom(S0,S0).identity()
            sage: f.is_pointed()
            True
            sage: v = S0.n_cells(0)[0]
            sage: w = S0.n_cells(0)[1]
            sage: g = Hom(S0,S0)({v:v, w:v})
            sage: g.is_pointed()
            True
            sage: t = Hom(S0,S0)({v:w, w:v})
            sage: t.is_pointed()
            False
        """
        return (self.domain().is_pointed() and self.codomain().is_pointed()
                and self(self.domain().base_point()) == self.codomain().base_point())

    def is_constant(self):
        """
        Return ``True`` if this morphism is a constant map.

        EXAMPLES::

            sage: K = simplicial_sets.KleinBottle()
            sage: S4 = simplicial_sets.Sphere(4)
            sage: c = Hom(K, S4).constant_map()
            sage: c.is_constant()
            True
            sage: X = S4.n_skeleton(3) # a point
            sage: X.inclusion_map().is_constant()
            True
            sage: eta = simplicial_sets.HopfMap()
            sage: eta.is_constant()
            False
        """
        try:
            return self._constant is not None
        except AttributeError:
            pass
        if not self.domain().is_finite():
            # The domain is infinite, so there is no safe way to
            # determine if the map is constant.
            return False
        targets = [tau.nondegenerate() for tau in self._dictionary.values()]
        if len(set(targets)) == 1:
            # It's constant, so save the target.
            self._constant = targets[0]
            return True
        return False

    def pushout(self, *others):
        """
        Return the pushout of this morphism along with ``others``.

        INPUT:

        - ``others`` -- morphisms of simplicial sets, the domains of
          which must all equal that of ``self``.

        This returns the pushout as a simplicial set. See
        :class:`sage.topology.simplicial_set_constructions.PushoutOfSimplicialSets`
        for more documentation and examples.

        EXAMPLES::

            sage: T = simplicial_sets.Torus()
            sage: K = simplicial_sets.KleinBottle()
            sage: init_T = T._map_from_empty_set()
            sage: init_K = K._map_from_empty_set()
            sage: D = init_T.pushout(init_K) # the disjoint union as a pushout
            sage: D
            Pushout of maps:
              Simplicial set morphism:
                From: Empty simplicial set
                To:   Torus
                Defn: [] --> []
              Simplicial set morphism:
                From: Empty simplicial set
                To:   Klein bottle
                Defn: [] --> []
        """
        domain = self.domain()
        if any(domain != f.domain() for f in others):
            raise ValueError('the domains of the maps must be equal')
        return self.domain().pushout(*(self,) + others)

    def pullback(self, *others):
        """
        Return the pullback of this morphism along with ``others``.

        INPUT:

        - ``others`` -- morphisms of simplicial sets, the codomains of
          which must all equal that of ``self``.

        This returns the pullback as a simplicial set. See
        :class:`sage.topology.simplicial_set_constructions.PullbackOfSimplicialSets`
        for more documentation and examples.

        EXAMPLES::

            sage: T = simplicial_sets.Torus()
            sage: K = simplicial_sets.KleinBottle()
            sage: term_T = T.constant_map()
            sage: term_K = K.constant_map()
            sage: P = term_T.pullback(term_K) # the product as a pullback
            sage: P
            Pullback of maps:
              Simplicial set morphism:
                From: Torus
                To:   Point
                Defn: Constant map at *
              Simplicial set morphism:
                From: Klein bottle
                To:   Point
                Defn: Constant map at *
        """
        codomain = self.codomain()
        if any(codomain != f.codomain() for f in others):
            raise ValueError('the codomains of the maps must be equal')
        return self.codomain().pullback(*(self,) + others)

    def equalizer(self, other):
        r"""
        Return the equalizer of this map with ``other``.

        INPUT:

        - ``other`` -- a morphism with the same domain and codomain as this map

        If the two maps are `f, g: X \to Y`, then the equalizer `P` is
        constructed as the pullback ::

            P ----> X
            |       |
            V       V
            X --> X x Y

        where the two maps `X \to X \times Y` are `(1,f)` and `(1,g)`.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: w = AbstractSimplex(0, name='w')
            sage: x = AbstractSimplex(0, name='x')
            sage: evw = AbstractSimplex(1, name='vw')
            sage: evx = AbstractSimplex(1, name='vx')
            sage: ewx = AbstractSimplex(1, name='wx')
            sage: X = SimplicialSet({evw: (w, v), evx: (x, v)})
            sage: Y = SimplicialSet({evw: (w, v), evx: (x, v), ewx: (x, w)})

        Here `X` is a wedge of two 1-simplices (a horn, that is), and
        `Y` is the boundary of a 2-simplex. The map `f` includes the
        two 1-simplices into `Y`, while the map `g` maps both
        1-simplices to the same edge in `Y`. ::

            sage: f = Hom(X, Y)({v:v, w:w, x:x, evw:evw, evx:evx})
            sage: g = Hom(X, Y)({v:v, w:x, x:x, evw:evx, evx:evx})
            sage: P = f.equalizer(g)
            sage: P
            Pullback of maps:
              Simplicial set morphism:
                From: Simplicial set with 5 non-degenerate simplices
                To:   Simplicial set with 5 non-degenerate simplices x Simplicial set with 6 non-degenerate simplices
                Defn: [v, w, x, vw, vx] --> [(v, v), (w, w), (x, x), (vw, vw), (vx, vx)]
              Simplicial set morphism:
                From: Simplicial set with 5 non-degenerate simplices
                To:   Simplicial set with 5 non-degenerate simplices x Simplicial set with 6 non-degenerate simplices
                Defn: [v, w, x, vw, vx] --> [(v, v), (w, x), (x, x), (vw, vx), (vx, vx)]
        """
        domain = self.domain()
        codomain = self.codomain()
        if domain != other.domain() or codomain != other.codomain():
            raise ValueError('the maps must have the same domain and the same codomain')
        prod = domain.product(codomain)
        one = domain.Hom(domain).identity()
        f = prod.universal_property(one, self)
        g = prod.universal_property(one, other)
        return f.pullback(g)

    def coequalizer(self, other):
        r"""
        Return the coequalizer of this map with ``other``.

        INPUT:

        - ``other`` -- a morphism with the same domain and codomain as this map

        If the two maps are `f, g: X \to Y`, then the coequalizer `P` is
        constructed as the pushout ::

            X v Y --> Y
              |       |
              V       V
              Y ----> P

        where the upper left corner is the coproduct of `X` and `Y`
        (the wedge if they are pointed, the disjoint union otherwise),
        and the two maps `X \amalg Y \to Y` are `f \amalg 1` and `g
        \amalg 1`.

        EXAMPLES::

            sage: L = simplicial_sets.Simplex(2)
            sage: pt = L.n_cells(0)[0]
            sage: e = L.n_cells(1)[0]
            sage: K = L.subsimplicial_set([e])
            sage: f = K.inclusion_map()
            sage: v,w = K.n_cells(0)
            sage: g = Hom(K,L)({v:pt, w:pt, e:pt.apply_degeneracies(0)})
            sage: P = f.coequalizer(g)
            sage: P
            Pushout of maps:
              Simplicial set morphism:
                From: Disjoint union: (Simplicial set with 3 non-degenerate simplices u 2-simplex)
                To:   2-simplex
                Defn: ...
              Simplicial set morphism:
                From: Disjoint union: (Simplicial set with 3 non-degenerate simplices u 2-simplex)
                To:   2-simplex
                Defn: ...
        """
        domain = self.domain()
        codomain = self.codomain()
        if domain != other.domain() or codomain != other.codomain():
            raise ValueError('the maps must have the same domain and the same codomain')
        coprod = domain.coproduct(codomain)
        one = codomain.Hom(codomain).identity()
        f = coprod.universal_property(self, one)
        g = coprod.universal_property(other, one)
        return f.pushout(g)

    def mapping_cone(self):
        r"""
        Return the mapping cone defined by this map.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: v_0, sigma_1 = S1.nondegenerate_simplices()
            sage: K = simplicial_sets.Simplex(2).n_skeleton(1)

        The mapping cone will be a little smaller if we use only
        pointed simplicial sets. `S^1` is already pointed, but not
        `K`. ::

            sage: L = K.set_base_point(K.n_cells(0)[0])
            sage: u,v,w = L.n_cells(0)
            sage: e,f,g = L.n_cells(1)
            sage: h = L.Hom(S1)({u:v_0, v:v_0, w:v_0, e:sigma_1, f:v_0.apply_degeneracies(0), g:sigma_1})
            sage: h
            Simplicial set morphism:
              From: Simplicial set with 6 non-degenerate simplices
              To:   S^1
              Defn: [(0,), (1,), (2,), (0, 1), (0, 2), (1, 2)] --> [v_0, v_0, v_0, sigma_1, s_0 v_0, sigma_1]
            sage: h.induced_homology_morphism().to_matrix()
            [1|0]
            [-+-]
            [0|2]
            sage: X = h.mapping_cone()
            sage: X.homology() == simplicial_sets.RealProjectiveSpace(2).homology()
            True
        """
        dom = self.domain()
        cone = dom.cone()
        i = cone.map_from_base()
        return self.pushout(i)

    def product(self, *others):
        r"""
        Return the product of this map with ``others``.

        - ``others`` -- morphisms of simplicial sets.

        If the relevant maps are `f_i: X_i \to Y_i`, this returns the
        natural map `\prod X_i \to \prod Y_i`.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: f = Hom(S1,S1).identity()
            sage: f.product(f).is_bijective()
            True
            sage: g = S1.constant_map(S1)
            sage: g.product(g).is_bijective()
            False
        """
        domain = self.domain().product(*[g.domain() for g in others])
        codomain = self.codomain().product(*[g.codomain() for g in others])
        factors = []
        for (i,f) in enumerate([self] + list(others)):
            factors.append(f * domain.projection_map(i))
        return codomain.universal_property(*factors)

    def coproduct(self, *others):
        r"""
        Return the coproduct of this map with ``others``.

        - ``others`` -- morphisms of simplicial sets.

        If the relevant maps are `f_i: X_i \to Y_i`, this returns the
        natural map `\amalg X_i \to \amalg Y_i`.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: f = Hom(S1,S1).identity()
            sage: f.coproduct(f).is_bijective()
            True
            sage: g = S1.constant_map(S1)
            sage: g.coproduct(g).is_bijective()
            False
        """
        codomain = self.codomain().coproduct(*[g.codomain() for g in others])
        factors = []
        for i, f in enumerate([self] + list(others)):
            factors.append(codomain.inclusion_map(i) * f)
        return codomain.universal_property(*factors)

    def suspension(self, n=1):
        """
        Return the `n`-th suspension of this morphism of simplicial sets.

        INPUT:

        - ``n`` (optional) -- non-negative integer, default 1

        EXAMPLES::

            sage: eta = simplicial_sets.HopfMap()
            sage: susp_eta = eta.suspension()
            sage: susp_eta.mapping_cone().homology() == eta.mapping_cone().suspension().homology()
            True

        This uses reduced suspensions if the original morphism is
        pointed, unreduced otherwise. So for example, if a constant
        map is not pointed, its suspension is not a constant map::

            sage: L = simplicial_sets.Simplex(1)
            sage: L.constant_map().is_pointed()
            False
            sage: f = L.constant_map().suspension()
            sage: f.is_constant()
            False

            sage: K = simplicial_sets.Sphere(3)
            sage: K.constant_map().is_pointed()
            True
            sage: g = K.constant_map().suspension()
            sage: g.is_constant()
            True

            sage: h = K.identity().suspension()
            sage: h.is_identity()
            True
        """
        domain = self.domain()
        codomain = self.codomain()
        if not self.is_pointed():
            # Make sure to use unreduced suspensions for both domain
            # and codomain.
            if domain.is_pointed():
                domain = domain.unset_base_point()
            if codomain.is_pointed():
                codomain = codomain.unset_base_point()
        f = self
        for i in range(n):
            new_dom = domain.suspension()
            new_cod = codomain.suspension()
            data = {new_dom.base_point(): new_cod.base_point()}
            for sigma in f._dictionary:
                target = f(sigma)
                underlying = target.nondegenerate()
                degens = target.degeneracies()
                data[new_dom._suspensions[sigma]] = new_cod._suspensions[underlying].apply_degeneracies(*degens)
            f = new_dom.Hom(new_cod)(data)
            domain = f.domain()
            codomain = f.codomain()
        return f

    def n_skeleton(self, n, domain=None, codomain=None):
        """
        Return the restriction of this morphism to the n-skeleta of the
        domain and codomain

        INPUT:

        - ``n`` -- the dimension

        - ``domain`` -- optional, the domain. Specify this to
          explicitly specify the domain; otherwise, Sage will attempt
          to compute it. Specifying this can be useful if the domain
          is built as a pushout or pullback, so trying to compute it
          may lead to computing the `n`-skeleton of a map, causing an
          infinite recursion. (Users should not have to specify this,
          but it may be useful for developers.)

        - ``codomain`` -- optional, the codomain.

        EXAMPLES::

            sage: B = simplicial_sets.ClassifyingSpace(groups.misc.MultiplicativeAbelian([2]))
            sage: one = Hom(B,B).identity()
            sage: one.n_skeleton(3)
            Simplicial set endomorphism of Simplicial set with 4 non-degenerate simplices
              Defn: Identity map
            sage: c = Hom(B,B).constant_map()
            sage: c.n_skeleton(3)
            Simplicial set endomorphism of Simplicial set with 4 non-degenerate simplices
              Defn: Constant map at 1

            sage: K = simplicial_sets.Simplex(2)
            sage: L = K.subsimplicial_set(K.n_cells(0)[:2])
            sage: L.nondegenerate_simplices()
            [(0,), (1,)]
            sage: L.inclusion_map()
            Simplicial set morphism:
              From: Simplicial set with 2 non-degenerate simplices
              To:   2-simplex
              Defn: [(0,), (1,)] --> [(0,), (1,)]
            sage: L.inclusion_map().n_skeleton(1)
            Simplicial set morphism:
              From: Simplicial set with 2 non-degenerate simplices
              To:   Simplicial set with 6 non-degenerate simplices
              Defn: [(0,), (1,)] --> [(0,), (1,)]
        """
        if domain is None:
            domain = self.domain().n_skeleton(n)
        if codomain is None:
            codomain = self.codomain().n_skeleton(n)
        if self.is_constant():
            return Hom(domain, codomain).constant_map(self._constant)
        if self.is_identity():
            return Hom(domain, domain).identity()
        old = self._dictionary
        new = {d: old[d] for d in old if d.dimension() <= n}
        return Hom(domain, codomain)(new)

    def associated_chain_complex_morphism(self, base_ring=ZZ,
                                          augmented=False, cochain=False):
        """
        Return the associated chain complex morphism of ``self``.

        INPUT:

        - ``base_ring`` -- default ``ZZ``
        - ``augmented`` -- boolean, default ``False``. If ``True``,
          return the augmented complex.
        - ``cochain`` -- boolean, default ``False``. If ``True``,
          return the cochain complex.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: v0 = S1.n_cells(0)[0]
            sage: e = S1.n_cells(1)[0]
            sage: f = {v0: v0, e: v0.apply_degeneracies(0)} # constant map
            sage: g = Hom(S1, S1)(f)
            sage: g.associated_chain_complex_morphism().to_matrix()
            [1|0]
            [-+-]
            [0|0]
        """
        # One or the other chain complex is trivial between these
        # dimensions:
        max_dim = max(self.domain().dimension(), self.codomain().dimension())
        min_dim = min(self.domain().dimension(), self.codomain().dimension())
        matrices = {}
        if augmented is True:
            m = matrix(base_ring,1,1,1)
            if not cochain:
                matrices[-1] = m
            else:
                matrices[-1] = m.transpose()
        for dim in range(min_dim+1):
            X_faces = list(self.domain().n_cells(dim))
            Y_faces = list(self.codomain().n_cells(dim))
            num_faces_X = len(X_faces)
            num_faces_Y = len(Y_faces)
            mval = [0 for _ in range(num_faces_X * num_faces_Y)]
            for idx,x in enumerate(X_faces):
                y = self(x)
                if y.is_nondegenerate():
                    mval[idx + (Y_faces.index(y) * num_faces_X)] = 1
            m = matrix(base_ring, num_faces_Y, num_faces_X, mval, sparse=True)
            if not cochain:
                matrices[dim] = m
            else:
                matrices[dim] = m.transpose()
        for dim in range(min_dim+1,max_dim+1):
            try:
                l1 = len(self.codomain().n_cells(dim))
            except KeyError:
                l1 = 0
            try:
                l2 = len(self.domain().n_cells(dim))
            except KeyError:
                l2 = 0
            m = zero_matrix(base_ring,l1,l2,sparse=True)
            if not cochain:
                matrices[dim] = m
            else:
                matrices[dim] = m.transpose()
        if not cochain:
            return ChainComplexMorphism(matrices,
                    self.domain().chain_complex(base_ring=base_ring, augmented=augmented, cochain=False),
                    self.codomain().chain_complex(base_ring=base_ring, augmented=augmented, cochain=False))
        else:
            return ChainComplexMorphism(matrices,
                    self.codomain().chain_complex(base_ring=base_ring, augmented=augmented, cochain=True),
                    self.domain().chain_complex(base_ring=base_ring, augmented=augmented, cochain=True))

    def induced_homology_morphism(self, base_ring=None, cohomology=False):
        """
        Return the map in (co)homology induced by this map

        INPUT:

        - ``base_ring`` -- must be a field (optional, default ``QQ``)

        - ``cohomology`` -- boolean (optional, default ``False``). If
          ``True``, the map induced in cohomology rather than homology.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: w = AbstractSimplex(0, name='w')
            sage: e = AbstractSimplex(1, name='e')
            sage: f = AbstractSimplex(1, name='f')
            sage: X = SimplicialSet({e: (v, w), f: (w, v)})
            sage: Y = SimplicialSet({e: (v, v)})
            sage: H = Hom(X, Y)
            sage: f = H({v: v, w: v, e: e, f: e})
            sage: g = f.induced_homology_morphism()
            sage: g.to_matrix()
            [1|0]
            [-+-]
            [0|2]
            sage: g3 = f.induced_homology_morphism(base_ring=GF(3), cohomology=True)
            sage: g3.to_matrix()
            [1|0]
            [-+-]
            [0|2]
        """
        return InducedHomologyMorphism(self, base_ring, cohomology)

    def _repr_type(self):
        """
        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: f = Hom(S1,S1).identity()
            sage: f._repr_type()
            'Simplicial set'
        """
        return "Simplicial set"

    def _repr_defn(self):
        """
        EXAMPLES::

            sage: K1 = simplicial_sets.Simplex(1)
            sage: v = K1.n_cells(0)[0]
            sage: e = K1.n_cells(1)[0]
            sage: f = Hom(K1,K1)({e:v.apply_degeneracies(0)})
            sage: f._repr_defn()
            'Constant map at (0,)'

            sage: K2 = simplicial_sets.Simplex(2)
            sage: tau = K2.n_cells(1)[0]
            sage: Hom(K1, K2)({e:tau})._repr_defn()
            '[(0,), (1,), (0, 1)] --> [(0,), (1,), (0, 1)]'

            sage: S1 = simplicial_sets.Sphere(1)
            sage: Hom(S1,S1).identity()._repr_defn()
            'Identity map'
        """
        if self.is_identity():
            return 'Identity map'
        if self.is_constant():
            return 'Constant map at {}'.format(self._constant)
        d = self._dictionary
        keys = sorted(d.keys())
        return "{} --> {}".format(keys, [d[x] for x in keys])

    def _latex_(self):
        """
        LaTeX representation.

        EXAMPLES::

            sage: eta = simplicial_sets.HopfMap()
            sage: eta.domain().rename_latex('S^{3}')
            sage: latex(eta)
            S^{3} \to S^{2}
        """
        return '{} \\to {}'.format(latex(self.domain()), latex(self.codomain()))
