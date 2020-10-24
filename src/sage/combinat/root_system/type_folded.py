r"""
Root system data for folded Cartan types

AUTHORS:

- Travis Scrimshaw (2013-01-12) - Initial version
"""
#*****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.sets.family import Family
from sage.combinat.root_system.cartan_type import CartanType


class CartanTypeFolded(UniqueRepresentation, SageObject):
    r"""
    A Cartan type realized from a (Dynkin) diagram folding.

    Given a Cartan type `X`, we say `\hat{X}` is a folded Cartan
    type of `X` if there exists a diagram folding of the Dynkin
    diagram of `\hat{X}` onto `X`.

    A folding of a simply-laced Dynkin diagram `D` with index set `I` is an
    automorphism `\sigma` of `D` where all nodes any orbit of `\sigma` are not
    connected. The resulting Dynkin diagram `\hat{D}` is induced by
    `I / \sigma` where we identify edges in `\hat{D}` which are not incident
    and add a `k`-edge if we identify `k` incident edges and the arrow is
    pointing towards the indicent note. We denote the index set of `\hat{D}`
    by `\hat{I}`, and by abuse of notation, we denote the folding by `\sigma`.

    We also have scaling factors `\gamma_i` for `i \in \hat{I}` and defined
    as the unique numbers such that the map
    `\Lambda_j \mapsto \gamma_j \sum_{i \in \sigma^{-1}(j)} \Lambda_i`
    is the smallest proper embedding of the weight lattice of `X` to `\hat{X}`.

    If the Cartan type is simply laced, the default folding is the one
    induced from the identity map on `D`.

    If `X` is affine type, the default embeddings we consider here are:

    .. MATH::

        \begin{array}{ccl}
        C_n^{(1)}, A_{2n}^{(2)}, A_{2n}^{(2)\dagger}, D_{n+1}^{(2)}
        & \hookrightarrow & A_{2n-1}^{(1)}, \\
        A_{2n-1}^{(2)}, B_n^{(1)} & \hookrightarrow & D_{n+1}^{(1)}, \\
        E_6^{(2)}, F_4^{(1)} & \hookrightarrow & E_6^{(1)}, \\
        D_4^{(3)}, G_2^{(1)} & \hookrightarrow & D_4^{(1)},
        \end{array}

    and were chosen based on virtual crystals. In particular, the diagram
    foldings extend to crystal morphisms and gives a realization of
    Kirillov-Reshetikhin crystals for non-simply-laced types as simply-laced
    types. See [OSShimo03]_ and [FOS2009]_ for more details. Here we can compute
    `\gamma_i = \max(c) / c_i` where `(c_i)_i` are the translation factors
    of the root system. In a more type-dependent way, we can define `\gamma_i`
    as follows:

    1. There exists a unique arrow (multiple bond) in `X`.

       a. Suppose the arrow points towards 0. Then `\gamma_i = 1` for all
          `i \in I`.
       b. Otherwise `\gamma_i` is the order of `\sigma` for all `i` in the
          connected component of 0 after removing the arrow, else
          `\gamma_i = 1`.

    2. There is not a unique arrow. Thus `\hat{X} = A_{2n-1}^{(1)}` and
       `\gamma_i = 1` for all `1 \leq i \leq n-1`. If `i \in \{0, n\}`, then
       `\gamma_i = 2` if the arrow incident to `i` points away and is `1`
       otherwise.

    We note that `\gamma_i` only depends upon `X`.

    If the Cartan type is finite, then we consider the classical
    foldings/embeddings induced by the above affine foldings/embeddings:

    .. MATH::

        \begin{aligned}
        C_n & \hookrightarrow A_{2n-1}, \\
        B_n & \hookrightarrow D_{n+1}, \\
        F_4 & \hookrightarrow E_6, \\
        G_2 & \hookrightarrow D_4.
        \end{aligned}

    For more information on Cartan types, see
    :mod:`sage.combinat.root_system.cartan_type`.

    Other foldings may be constructed by passing in an optional
    ``folding_of`` second argument. See below.

    INPUT:

    - ``cartan_type`` -- the Cartan type `X` to create the folded type

    - ``folding_of`` -- the Cartan type `\hat{X}` which `X` is a folding of

    - ``orbit`` -- the orbit of the Dynkin diagram automorphism `\sigma`
      given as a list of lists where the `a`-th list corresponds to the `a`-th
      entry in `I` or a dictionary with keys in `I` and values as lists

    .. NOTE::

        If `X` is an affine type, we assume the special node is fixed
        under `\sigma`.

    EXAMPLES::

        sage: fct = CartanType(['C',4,1]).as_folding(); fct
        ['C', 4, 1] as a folding of ['A', 7, 1]
        sage: fct.scaling_factors()
        Finite family {0: 2, 1: 1, 2: 1, 3: 1, 4: 2}
        sage: fct.folding_orbit()
        Finite family {0: (0,), 1: (1, 7), 2: (2, 6), 3: (3, 5), 4: (4,)}

    A simply laced Cartan type can be considered as a virtual type of
    itself::

        sage: fct = CartanType(['A',4,1]).as_folding(); fct
        ['A', 4, 1] as a folding of ['A', 4, 1]
        sage: fct.scaling_factors()
        Finite family {0: 1, 1: 1, 2: 1, 3: 1, 4: 1}
        sage: fct.folding_orbit()
        Finite family {0: (0,), 1: (1,), 2: (2,), 3: (3,), 4: (4,)}

    Finite types::

        sage: fct = CartanType(['C',4]).as_folding(); fct
        ['C', 4] as a folding of ['A', 7]
        sage: fct.scaling_factors()
        Finite family {1: 1, 2: 1, 3: 1, 4: 2}
        sage: fct.folding_orbit()
        Finite family {1: (1, 7), 2: (2, 6), 3: (3, 5), 4: (4,)}

        sage: fct = CartanType(['F',4]).dual().as_folding(); fct
        ['F', 4] relabelled by {1: 4, 2: 3, 3: 2, 4: 1} as a folding of ['E', 6]
        sage: fct.scaling_factors()
        Finite family {1: 1, 2: 1, 3: 2, 4: 2}
        sage: fct.folding_orbit()
        Finite family {1: (1, 6), 2: (3, 5), 3: (4,), 4: (2,)}

    REFERENCES:

    - :wikipedia:`Dynkin_diagram#Folding`

    .. [OSShimo03] \M. Okado, A. Schilling, M. Shimozono.
       "Virtual crystals and fermionic formulas for type `D_{n+1}^{(2)}`,
       `A_{2n}^{(2)}`, and `C_n^{(1)}`". Representation Theory. **7** (2003).
       101-163. :doi:`10.1.1.192.2095`, :arxiv:`0810.5067`.
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type, virtual, orbit):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: from sage.combinat.root_system.type_folded import CartanTypeFolded
            sage: sigma_list = [[0], [1,5], [2,4], [3]]
            sage: fct1 = CartanTypeFolded(['C',3,1], ['A',5,1], sigma_list)
            sage: sigma_tuple = tuple(map(tuple, sigma_list))
            sage: fct2 = CartanTypeFolded(CartanType(['C',3,1]), CartanType(['A',5,1]), sigma_tuple)
            sage: fct3 = CartanTypeFolded('C3~', 'A5~', {0:[0], 2:[2,4], 1:[1,5], 3:[3]})
            sage: fct1 is fct2 and fct2 is fct3
            True
        """
        if isinstance(cartan_type, CartanTypeFolded):
            return cartan_type
        cartan_type = CartanType(cartan_type)
        virtual = CartanType(virtual)
        if isinstance(orbit, dict):
            i_set = cartan_type.index_set()
            orb = [None]*len(i_set)
            for k,v in orbit.items():
                orb[i_set.index(k)] = tuple(v)
            orbit = tuple(orb)
        else:
            orbit = tuple(map(tuple, orbit))
        return super(CartanTypeFolded, cls).__classcall__(cls, cartan_type, virtual, orbit)

    def __init__(self, cartan_type, folding_of, orbit):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: fct = CartanType(['C',4,1]).as_folding()
            sage: TestSuite(fct).run()
            sage: hash(fct)  # random
            42
        """
        self._cartan_type = cartan_type
        self._folding = folding_of
        self._orbit = orbit

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: CartanType(['C',4,1]).as_folding()
            ['C', 4, 1] as a folding of ['A', 7, 1]
        """
        return "{} as a folding of {}".format(self._cartan_type, self._folding)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: fct = CartanType(['C', 4, 1]).as_folding()
            sage: latex(fct)
            C_{4}^{(1)} \hookrightarrow A_{7}^{(1)}
        """
        return self._cartan_type._latex_() + " \\hookrightarrow " + self._folding._latex_()

    def cartan_type(self):
        """
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: fct = CartanType(['C', 4, 1]).as_folding()
            sage: fct.cartan_type()
            ['C', 4, 1]
        """
        return self._cartan_type

    def folding_of(self):
        """
        Return the Cartan type of the virtual space.

        EXAMPLES::

            sage: fct = CartanType(['C', 4, 1]).as_folding()
            sage: fct.folding_of()
            ['A', 7, 1]
        """
        return self._folding

    @cached_method
    def folding_orbit(self):
        r"""
        Return the orbits under the automorphism `\sigma` as a
        dictionary (of tuples).

        EXAMPLES::

            sage: fct = CartanType(['C', 4, 1]).as_folding()
            sage: fct.folding_orbit()
            Finite family {0: (0,), 1: (1, 7), 2: (2, 6), 3: (3, 5), 4: (4,)}
        """
        return Family({i:tuple(self._orbit[pos])
                       for pos,i in enumerate(self._cartan_type.index_set())})

    @cached_method
    def scaling_factors(self):
        """
        Return the scaling factors of ``self``.

        EXAMPLES::

            sage: fct = CartanType(['C', 4, 1]).as_folding()
            sage: fct.scaling_factors()
            Finite family {0: 2, 1: 1, 2: 1, 3: 1, 4: 2}
            sage: fct = CartanType(['BC', 4, 2]).as_folding()
            sage: fct.scaling_factors()
            Finite family {0: 1, 1: 1, 2: 1, 3: 1, 4: 2}
            sage: fct = CartanType(['BC', 4, 2]).dual().as_folding()
            sage: fct.scaling_factors()
            Finite family {0: 2, 1: 1, 2: 1, 3: 1, 4: 1}
            sage: CartanType(['BC', 4, 2]).relabel({0:4, 1:3, 2:2, 3:1, 4:0}).as_folding().scaling_factors()
            Finite family {0: 2, 1: 1, 2: 1, 3: 1, 4: 1}
        """
        if self._cartan_type.is_finite():
            L = self._cartan_type.root_system().ambient_space()
            def f(i):
                root = L.simple_root(i)
                coroot = L.simple_coroot(i)
                return root.leading_coefficient() / coroot.leading_coefficient()
            index_set = self._cartan_type.index_set()
            min_f = min(f(j) for j in index_set)
            return Family(dict( (i, int(f(i) / min_f)) for i in index_set ))
        elif self._cartan_type.is_affine():
            c = self._cartan_type.translation_factors()
            cmax = max(c)
            return Family(dict( (i, int(cmax / c[i]))
                                for i in self._cartan_type.index_set() ))

