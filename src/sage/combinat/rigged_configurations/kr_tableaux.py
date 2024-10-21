# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
r"""
Kirillov-Reshetikhin Tableaux

Kirillov-Reshetikhin tableaux are rectangular tableaux with `r` rows and
`s` columns that naturally arise under the bijection between rigged
configurations and tableaux [RigConBijection]_. They are in bijection with
the elements of the Kirillov-Reshetikhin crystal `B^{r,s}` under the (inverse)
filling map [OSS13]_ [SS2015]_. They do not have to satisfy the semistandard row or column
restrictions. These tensor products are the result from the bijection from
rigged configurations [RigConBijection]_.

For more information, see :class:`~sage.combinat.rigged_configurations.kr_tableaux.KirillovReshetikhinTableaux`
and :class:`~sage.combinat.rigged_configurations.tensor_product_kr_tableaux.TensorProductOfKirillovReshetikhinTableaux`.

AUTHORS:

- Travis Scrimshaw (2012-01-03): Initial version
- Travis Scrimshaw (2012-11-14): Added bijection to KR crystals

REFERENCES:

.. [OSS13] Masato Okado, Reiho Sakamoto, and Anne Schilling.
   *Affine crystal structure on rigged configurations of type* `D_n^{(1)}`.
   J. Algebraic Combinatorics, **37** (2013). 571-599. :arxiv:`1109.3523`.
"""

#*****************************************************************************
#       Copyright (C) 2012 Travis Scrimshaw <tscrim@ucdavis.edu>
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# This contains both the parent and element classes. These should be split if
#   the classes grow larger.

from sage.misc.cachefunc import cached_method
from sage.misc.abstract_method import abstract_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.flatten import flatten

from sage.structure.parent import Parent

from sage.categories.loop_crystals import KirillovReshetikhinCrystals

from sage.combinat.crystals.letters import CrystalOfLetters, EmptyLetter
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.crystals.tensor_product import CrystalOfWords
from sage.combinat.crystals.tensor_product import TensorProductOfRegularCrystalsElement
from sage.combinat.crystals.kirillov_reshetikhin import horizontal_dominoes_removed, \
  KashiwaraNakashimaTableaux, KirillovReshetikhinGenericCrystalElement, \
  partitions_in_box, vertical_dominoes_removed
from sage.combinat.partition import Partition
from sage.combinat.tableau import Tableau


class KirillovReshetikhinTableaux(CrystalOfWords):
    r"""
    Kirillov-Reshetikhin tableaux.

    Kirillov-Reshetikhin tableaux are rectangular tableaux with `r` rows and
    `s` columns that naturally arise under the bijection between rigged
    configurations and tableaux [RigConBijection]_. They are in bijection with
    the elements of the Kirillov-Reshetikhin crystal `B^{r,s}` under the
    (inverse) filling map.

    Whenever `B^{r,s} \cong B(s\Lambda_r)` as a classical crystal (which is
    the case for `B^{r,s}` in type `A_n^{(1)}`, `B^{n,s}` in type `C_n^{(1)}`
    and `D_{n+1}^{(2)}`, `B^{n,s}` and `B^{n-1,s}` in type `D_n^{(1)}`) then
    the filling map is trivial.

    For `B^{r,s}` in:

    - type `D_n^{(1)}` when `r \leq n-2`,
    - type `B_n^{(1)}` when `r < n`,
    - type `A_{2n-1}^{(2)}` for all `r`,

    the filling map is defined in [OSS2011]_.

    For the spinor cases in type `D_n^{(1)}`, the crystal `B^{k,s}` where
    `k = n-1, n`,  is isomorphic as a classical crystal to `B(s\Lambda_k)`,
    and here we consider the Kirillov-Reshetikhin tableaux as living in
    `B(2s \Lambda_k)` under the natural doubling map. In this case, the
    crystal operators `e_i` and `f_i` act as `e_i^2` and `f_i^2` respectively.
    See [BijectionDn]_.

    For the spinor case in type `B_n^{(1)}`, the crystal `B^{n,s}`, we
    consider the images under the natural doubling map into `B^{n,2s}`.
    The classical components of this crystal are now given by
    removing `2 \times 2` boxes. The filling map is the same as below
    (see the non-spin type `C_n^{(1)}`).

    For `B^{r,s}` in:

    - type `C_n^{(1)}` when `r < n`,
    - type `A_{2n}^{(2)\dagger}` for all `r`,

    the filling map is given as follows. Suppose we are considering the
    (classically) highest weight element in the classical component
    `B(\lambda)`. Then we fill it in with the horizontal dominoes
    `[\bar{\imath}, i]` in the `i`-th row from the top (in English notation)
    and reordering the columns so that they are increasing. Recall from above
    that `B^{n,s} \cong B(s\Lambda_n)` in type `C^{(1)}_n`.

    For `B^{r,s}` in:

    - type `A_{2n}^{(2)}` for all `r`,
    - type `D_{n+1}^{(2)}` when `r < n`,
    - type `D_4^{(3)}` when `r = 1`,

    the filling map is the same as given in [OSS2011]_ except for
    the rightmost column which is given by the column `[1, 2, \ldots, k,
    \emptyset, \ldots \emptyset]` where `k = (r+x-1)/2` in Step 3 of
    [OSS2011]_.

    For the spinor case in type `D_{n+1}^{(2)}`, the crystal `B^{n,s}`, we
    define the filling map in the same way as in type `D_n^{(1)}`.

    .. NOTE::

        The filling map and classical decompositions in non-spinor cases can
        be classified by how the special node `0` connects with the
        corresponding classical diagram.

    The classical crystal structure is given by the usual Kashiwara-Nakashima
    tableaux rules. That is to embed this into `B(\Lambda_1)^{\otimes n s}`
    by using the reading word and then applying the classical crystal
    operator. The affine crystal structure is given by converting to
    the corresponding KR crystal element, performing the affine crystal
    operator, and pulling back to a KR tableau.

    For more information about the bijection between rigged configurations
    and tensor products of Kirillov-Reshetikhin tableaux, see
    :class:`~sage.combinat.rigged_configurations.tensor_product_kr_tableaux.TensorProductOfKirillovReshetikhinTableaux`.

    .. NOTE::

        The tableaux for all non-simply-laced types are provably correct if the
        bijection with :class:`rigged configurations
        <sage.combinat.rigged_configurations.rigged_configurations.RiggedConfigurations>`
        holds. Therefore this is currently only proven for `B^{r,1}` or
        `B^{1,s}` and in general for types `A_n^{(1)}` and `D_n^{(1)}`.

    INPUT:

    - ``cartan_type`` -- the Cartan type

    - ``r`` -- the Dynkin diagram index (typically the number of rows)

    - ``s`` -- the number of columns

    EXAMPLES::

        sage: KRT = crystals.KirillovReshetikhin(['A', 4, 1], 2, 1, model='KR')
        sage: elt = KRT(4, 3); elt
        [[3], [4]]

        sage: KRT = crystals.KirillovReshetikhin(['D', 4, 1], 2, 1, model='KR')
        sage: elt = KRT(-1, 1); elt
        [[1], [-1]]

    We can create highest weight crystals from a given shape or weight::

        sage: KRT = crystals.KirillovReshetikhin(['D', 4, 1], 2, 2, model='KR')
        sage: KRT.module_generator(shape=[1,1])
        [[1, 1], [2, -1]]
        sage: KRT.module_generator(column_shape=[2])
        [[1, 1], [2, -1]]
        sage: WS = RootSystem(['D',4,1]).weight_space()
        sage: KRT.module_generator(weight=WS.sum_of_terms([[0,-2],[2,1]]))
        [[1, 1], [2, -1]]
        sage: WSC = RootSystem(['D',4]).weight_space()
        sage: KRT.module_generator(classical_weight=WSC.fundamental_weight(2))
        [[1, 1], [2, -1]]

    We can go between
    :func:`~sage.combinat.crystals.kirillov_reshetikhin.KashiwaraNakashimaTableaux`
    and
    :class:`~sage.combinat.rigged_configurations.kr_tableaux.KirillovReshetikhinTableaux`
    elements::

        sage: KRCrys = crystals.KirillovReshetikhin(['D', 4, 1], 2, 2, model='KN')
        sage: KRTab = crystals.KirillovReshetikhin(['D', 4, 1], 2, 2, model='KR')
        sage: elt = KRCrys(3, 2); elt
        [[2], [3]]
        sage: k = KRTab(elt); k
        [[2, 1], [3, -1]]
        sage: KRCrys(k)
        [[2], [3]]

    We check that the classical weights in the classical decompositions
    agree in a few different type::

        sage: KRCrys = crystals.KirillovReshetikhin(['D', 4, 1], 2, 2, model='KN')
        sage: KRTab = crystals.KirillovReshetikhin(['D', 4, 1], 2, 2, model='KR')
        sage: all(t.classical_weight() == KRCrys(t).classical_weight() for t in KRTab)
        True
        sage: KRCrys = crystals.KirillovReshetikhin(['B', 3, 1], 2, 2, model='KN')
        sage: KRTab = crystals.KirillovReshetikhin(['B', 3, 1], 2, 2, model='KR')
        sage: all(t.classical_weight() == KRCrys(t).classical_weight() for t in KRTab)
        True
        sage: KRCrys = crystals.KirillovReshetikhin(['C', 3, 1], 2, 2, model='KN')
        sage: KRTab = crystals.KirillovReshetikhin(['C', 3, 1], 2, 2, model='KR')
        sage: all(t.classical_weight() == KRCrys(t).classical_weight() for t in KRTab)
        True
        sage: KRCrys = crystals.KirillovReshetikhin(['D', 4, 2], 2, 2, model='KN')
        sage: KRTab = crystals.KirillovReshetikhin(['D', 4, 2], 2, 2, model='KR')
        sage: all(t.classical_weight() == KRCrys(t).classical_weight() for t in KRTab)
        True
        sage: KRCrys = crystals.KirillovReshetikhin(['A', 4, 2], 2, 2, model='KN')
        sage: KRTab = crystals.KirillovReshetikhin(['A', 4, 2], 2, 2, model='KR')
        sage: all(t.classical_weight() == KRCrys(t).classical_weight() for t in KRTab)
        True
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type, r, s):
        """
        Normalize the input arguments to ensure unique representation.

        EXAMPLES::

            sage: KRT1 = crystals.KirillovReshetikhin(CartanType(['A',3,1]), 2, 3, model='KR')
            sage: KRT2 = crystals.KirillovReshetikhin(['A',3,1], 2, 3, model='KR')
            sage: KRT1 is KRT2
            True
        """
        ct = CartanType(cartan_type)
        if not ct.is_affine():
            raise ValueError("The Cartan type must be affine")

        typ = ct.type()
        if ct.is_untwisted_affine():
            if typ == 'A':
                return KRTableauxRectangle(ct, r, s)
            if typ == 'B':
                if r == ct.classical().rank():
                    return KRTableauxBn(ct, r, s)
                return KRTableauxTypeVertical(ct, r, s)
            if typ == 'C':
                if r == ct.classical().rank():
                    return KRTableauxRectangle(ct, r, s)
                return KRTableauxTypeHorizonal(ct, r, s)
            if typ == 'D':
                if r == ct.classical().rank() or r == ct.classical().rank() - 1:
                    return KRTableauxSpin(ct, r, s)
                return KRTableauxTypeVertical(ct, r, s)
            if typ == 'E':
                return KRTableauxTypeFromRC(ct, r, s)
        else:
            if typ == 'BC': # A_{2n}^{(2)}
                return KRTableauxTypeBox(ct, r, s)
            typ = ct.dual().type()
            if typ == 'BC': # A_{2n}^{(2)\dagger}
                return KRTableauxTypeHorizonal(ct, r, s)
            if typ == 'B': # A_{2n-1}^{(2)}
                return KRTableauxTypeVertical(ct, r, s)
            if typ == 'C': # D_{n+1}^{(2)}
                if r == ct.dual().classical().rank():
                    return KRTableauxDTwistedSpin(ct, r, s)
                return KRTableauxTypeBox(ct, r, s)
            #if typ == 'F': # E_6^{(2)}
            if typ == 'G': # D_4^{(3)}
                if r == 1:
                    return KRTableauxTypeBox(ct, r, s)
                return KRTableauxTypeFromRC(ct, r, s)

        raise NotImplementedError
        #return super(KirillovReshetikhinTableaux, cls).__classcall__(cls, ct, r, s)

    def __init__(self, cartan_type, r, s):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['A', 4, 1], 2, 2, model='KR')
            sage: TestSuite(KRT).run()
            sage: KRT = crystals.KirillovReshetikhin(['D', 4, 1], 2, 2, model='KR')
            sage: TestSuite(KRT).run()  # long time
            sage: KRT = crystals.KirillovReshetikhin(['D', 4, 1], 4, 1, model='KR'); KRT
            Kirillov-Reshetikhin tableaux of type ['D', 4, 1] and shape (4, 1)
            sage: TestSuite(KRT).run()
        """
        self._r = r
        self._s = s
        self._cartan_type = cartan_type

        Parent.__init__(self, category=KirillovReshetikhinCrystals())

        self.letters = CrystalOfLetters(cartan_type.classical())
        self.module_generators = self._build_module_generators()

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: crystals.KirillovReshetikhin(['A', 4, 1], 2, 3, model='KR')
            Kirillov-Reshetikhin tableaux of type ['A', 4, 1] and shape (2, 3)
        """
        return "Kirillov-Reshetikhin tableaux of type {} and shape ({}, {})".format(
                self._cartan_type, self._r, self._s)

    def __iter__(self):
        """
        Return the iterator of ``self``.

        EXAMPLES::

            sage: KR = crystals.KirillovReshetikhin(['A', 5, 2], 2, 1, model='KR')
            sage: L = [x for x in KR]
            sage: len(L)
            15
        """
        index_set = self._cartan_type.classical().index_set()
        from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet
        return RecursivelyEnumeratedSet(self.module_generators,
                    lambda x: [x.f(i) for i in index_set],
                    structure='graded').breadth_first_search_iterator()

    def module_generator(self, i=None, **options):
        r"""
        Return the specified module generator.

        INPUT:

        - ``i`` -- the index of the module generator

        We can also get a module generator by using one of the following
        optional arguments:

        - ``shape`` -- the associated shape
        - ``column_shape`` -- the shape given as columns (a column of length
          `k` correspond to a classical weight `\omega_k`)
        - ``weight`` -- the weight
        - ``classical_weight`` -- the classical weight

        If no arguments are specified, then return the unique module generator
        of classical weight `s \Lambda_r`.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D', 4, 1], 2, 2, model='KR')
            sage: KRT.module_generator(1)
            [[1, 1], [2, -1]]
            sage: KRT.module_generator(shape=[1,1])
            [[1, 1], [2, -1]]
            sage: KRT.module_generator(column_shape=[2])
            [[1, 1], [2, -1]]
            sage: WS = RootSystem(['D',4,1]).weight_space()
            sage: KRT.module_generator(weight=WS.sum_of_terms([[0,-2],[2,1]]))
            [[1, 1], [2, -1]]
            sage: WSC = RootSystem(['D',4]).weight_space()
            sage: KRT.module_generator(classical_weight=WSC.fundamental_weight(2))
            [[1, 1], [2, -1]]
            sage: KRT.module_generator()
            [[1, 1], [2, 2]]

            sage: KRT = crystals.KirillovReshetikhin(['A', 3, 1], 2, 2, model='KR')
            sage: KRT.module_generator()
            [[1, 1], [2, 2]]
        """
        if i is not None:
            return self.module_generators[i]
        n = self._cartan_type.classical().rank()

        if "shape" in options:
            shape = list(options["shape"])
            # Make sure the shape is the correct length
            if len(shape) < n:
                shape.extend( [0]*(n - len(shape)) )
            for mg in self.module_generators:
                if list(mg.classical_weight().to_vector()) == shape:
                    return mg
            return None

        if "column_shape" in options:
            shape = list(Partition(options["column_shape"]).conjugate())
            if len(shape) < n:
                shape.extend( [0]*(n - len(shape)) )
            for mg in self.module_generators:
                if list(mg.classical_weight().to_vector()) == shape:
                    return mg
            return None

        if "weight" in options:
            wt = options["weight"]
            for mg in self.module_generators:
                if mg.weight() == wt:
                    return mg
            return None

        if "classical_weight" in options:
            wt = options["classical_weight"]
            for mg in self.module_generators:
                if mg.classical_weight() == wt:
                    return mg
            return None

        # Otherwise return the unique module generator of classical weight `s \Lambda_r`
        R = self.weight_lattice_realization()
        Lambda = R.fundamental_weights()
        r = self.r()
        s = self.s()
        weight = s*Lambda[r] - s*Lambda[0] * Lambda[r].level() / Lambda[0].level()
        for b in self.module_generators:
            if b.weight() == weight:
                return b
        assert False

    @abstract_method
    def _build_module_generators(self):
        """
        Build the module generators.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['A', 4, 1], 2, 3, model='KR')
            sage: KRT._build_module_generators()
            ([[1, 1, 1], [2, 2, 2]],)
        """

    @abstract_method(optional=True)
    def from_kirillov_reshetikhin_crystal(self, krc):
        """
        Construct an element of ``self`` from the Kirillov-Reshetikhin
        crystal element ``krc``.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['A', 4, 1], 2, 1, model='KR')
            sage: C = crystals.KirillovReshetikhin(['A',4,1], 2, 1, model='KN')
            sage: krc = C(4,3); krc
            [[3], [4]]
            sage: KRT.from_kirillov_reshetikhin_crystal(krc)
            [[3], [4]]
        """

    def _element_constructor_(self, *lst, **options):
        """
        Construct a
        :class:`~sage.combinat.rigged_configurations.kr_tableaux.KirillovReshetikhinTableauxElement`.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['A', 4, 1], 2, 1, model='KR')
            sage: KRT(3, 4) # indirect doctest
            [[4], [3]]
            sage: KRT(4, 3)
            [[3], [4]]
        """
        if isinstance(lst[0], KirillovReshetikhinGenericCrystalElement):
            # Check to make sure it can be converted
            if lst[0].cartan_type() != self.cartan_type() \
              or lst[0].parent().r() != self._r or lst[0].parent().s() != self._s:
                raise ValueError("the Kirillov-Reshetikhin crystal must have the same Cartan type and (r,s)")
            return self.from_kirillov_reshetikhin_crystal(lst[0])

        return self.element_class(self, list(lst), **options)

    def r(self):
        """
        Return the value `r` for this tableaux class which corresponds to the
        number of rows.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['A', 4, 1], 2, 1, model='KR')
            sage: KRT.r()
            2
        """
        return self._r

    def s(self):
        """
        Return the value `s` for this tableaux class which corresponds to the
        number of columns.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['A', 4, 1], 2, 1, model='KR')
            sage: KRT.s()
            1
        """
        return self._s

    @cached_method
    def kirillov_reshetikhin_crystal(self):
        """
        Return the corresponding KR crystal in the
        :func:`Kashiwara-Nakashima model
        <sage.combinat.crystals.kirillov_reshetikhin.KashiwaraNakashimaTableaux>`.

        EXAMPLES::

            sage: crystals.KirillovReshetikhin(['A', 4, 1], 2, 1, model='KR').kirillov_reshetikhin_crystal()
            Kirillov-Reshetikhin crystal of type ['A', 4, 1] with (r,s)=(2,1)
        """
        return KashiwaraNakashimaTableaux(self._cartan_type, self._r, self._s)

    def classical_decomposition(self):
        """
        Return the classical crystal decomposition of ``self``.

        EXAMPLES::

            sage: crystals.KirillovReshetikhin(['D', 4, 1], 2, 2, model='KR').classical_decomposition()
            The crystal of tableaux of type ['D', 4] and shape(s) [[], [1, 1], [2, 2]]
        """
        return self.kirillov_reshetikhin_crystal().classical_decomposition()

    def tensor(self, *crystals, **options):
        """
        Return the tensor product of ``self`` with ``crystals``.

        If ``crystals`` is a list of (a tensor product of) KR tableaux, this
        returns a
        :class:`~sage.combinat.rigged_configurations.tensor_product_kr_tableaux.TensorProductOfKirillovReshetikhinTableaux`.

        EXAMPLES::

            sage: K = crystals.KirillovReshetikhin(['A', 3, 1], 2, 2, model='KR')
            sage: TP = crystals.TensorProductOfKirillovReshetikhinTableaux(['A', 3, 1], [[1,3],[3,1]])
            sage: K.tensor(TP, K)
            Tensor product of Kirillov-Reshetikhin tableaux of type ['A', 3, 1]
             and factor(s) ((2, 2), (1, 3), (3, 1), (2, 2))

            sage: C = crystals.KirillovReshetikhin(['A',3,1], 3, 1, model='KN')
            sage: K.tensor(K, C)
            Full tensor product of the crystals
             [Kirillov-Reshetikhin tableaux of type ['A', 3, 1] and shape (2, 2),
              Kirillov-Reshetikhin tableaux of type ['A', 3, 1] and shape (2, 2),
              Kirillov-Reshetikhin crystal of type ['A', 3, 1] with (r,s)=(3,1)]
        """
        ct = self._cartan_type
        from sage.combinat.rigged_configurations.tensor_product_kr_tableaux \
                import TensorProductOfKirillovReshetikhinTableaux
        if all(isinstance(B, (KirillovReshetikhinTableaux, TensorProductOfKirillovReshetikhinTableaux))
               and B.cartan_type() == ct for B in crystals):
            dims = [[self._r, self._s]]
            for B in crystals:
                if isinstance(B, TensorProductOfKirillovReshetikhinTableaux):
                    dims += B.dims
                elif isinstance(B, KirillovReshetikhinTableaux):
                    dims.append([B._r, B._s])
            return TensorProductOfKirillovReshetikhinTableaux(ct, dims)
        return super().tensor(*crystals, **options)

    @lazy_attribute
    def _tableau_height(self):
        """
        The height of the tableaux in ``self``.

        EXAMPLES::

            sage: K = crystals.KirillovReshetikhin(['A', 3, 1], 3, 2, model='KR')
            sage: K._tableau_height
            3
        """
        return self._r


class KRTableauxRectangle(KirillovReshetikhinTableaux):
    r"""
    Kirillov-Reshetkhin tableaux `B^{r,s}` whose module generator is a single
    `r \times s` rectangle.

    These are Kirillov-Reshetkhin tableaux `B^{r,s}` of type:

    - `A_n^{(1)}` for all `1 \leq r \leq n`,
    - `C_n^{(1)}` when `r = n`.

    TESTS::

        sage: KRT = crystals.KirillovReshetikhin(['A', 3, 1], 2, 2, model='KR')
        sage: TestSuite(KRT).run()
        sage: KRT = crystals.KirillovReshetikhin(['C', 3, 1], 3, 2, model='KR')
        sage: TestSuite(KRT).run() # long time
    """

    def _build_module_generators(self):
        r"""
        Build the module generators.

        There is only one module generator which corresponds to a single
        `r \times s` rectangle.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['A', 4, 1], 2, 3, model='KR')
            sage: KRT._build_module_generators()
            ([[1, 1, 1], [2, 2, 2]],)
        """
        tableau = []
        for i in range(self._s):
            tableau.append( [self._r - j for j in range(self._r)] )

        return (self.element_class(self, [self.letters(x) for x in flatten(tableau)]),)

    def from_kirillov_reshetikhin_crystal(self, krc):
        """
        Construct a
        :class:`~sage.combinat.rigged_configurations.kr_tableaux.KirillovReshetikhinTableauxElement`.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['A', 4, 1], 2, 1, model='KR')
            sage: C = crystals.KirillovReshetikhin(['A',4,1], 2, 1, model='KN')
            sage: krc = C(4,3); krc
            [[3], [4]]
            sage: KRT.from_kirillov_reshetikhin_crystal(krc)
            [[3], [4]]
        """
        # To build a KR tableau from a KR crystal:
        # 1 - start with the highest weight KR tableau
        # 2 - determine a path from the KR crystal to its highest weight
        # 3 - apply the inverse path to the highest weight KR tableau
        f_str = reversed(krc.lift().to_highest_weight()[1])
        return self.module_generators[0].f_string(f_str)


class KRTableauxTypeVertical(KirillovReshetikhinTableaux):
    r"""
    Kirillov-Reshetkihn tableaux `B^{r,s}` of type:

    - `D_n^{(1)}` for all `1 \leq r < n-1`,
    - `B_n^{(1)}` for all `1 \leq r < n`,
    - `A_{2n-1}^{(2)}` for all `1 \leq r \leq n`.

    TESTS::

        sage: KRT = crystals.KirillovReshetikhin(['D', 4, 1], 1, 1, model='KR')
        sage: TestSuite(KRT).run()
        sage: KRT = crystals.KirillovReshetikhin(['B', 3, 1], 2, 2, model='KR')
        sage: TestSuite(KRT).run() # long time
        sage: KRT = crystals.KirillovReshetikhin(['A', 5, 2], 2, 2, model='KR')
        sage: TestSuite(KRT).run() # long time
    """

    def _fill(self, weight):
        r"""
        Return the highest weight KR tableau of weight ``weight``.

        INPUT:

        - ``weight`` -- the weight of the highest weight KR tableau (the
          conjugate of the shape of the KR crystal's tableau)

        OUTPUT: a `r \times s` tableau

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D', 4, 1], 2, 1, model='KR')
            sage: KRT._fill([])
            [[1], [-1]]
            sage: KRT = crystals.KirillovReshetikhin(['D', 14, 1], 12, 7, model='KR')
            sage: KRT._fill([10,10,8,2,2,2])
            [[1, 1, 1, 1, 1, 7, 1], [2, 2, 2, 2, 2, 8, 2], [3, 3, 7, 9, 7, 9, 3], [4, 4, 8, 10, 8, 10, 4], [5, 5, 9, 11, 9, 11, 5], [6, 6, 10, 12, 10, 12, 6], [7, 7, 11, -12, 11, -12, 7], [8, 8, 12, -11, 12, -11, 8], [9, 9, -12, -10, -12, -10, 9], [10, 10, -11, -9, -11, -9, -9], [-12, 11, -10, -8, -10, -8, -8], [-11, 12, -9, -7, -9, -7, -7]]
            sage: KRT._fill([10,10,6,2,2,2])
            [[1, 1, 1, 1, 1, 5, 1], [2, 2, 2, 2, 2, 6, 2], [3, 3, 9, 7, 9, 7, 3], [4, 4, 10, 8, 10, 8, 4], [5, 5, 11, 9, 11, 9, 5], [6, 6, 12, 10, 12, 10, 6], [7, 7, -12, 11, -12, 11, 7], [8, 8, -11, 12, -11, 12, 8], [9, 9, -10, -12, -10, -12, -8], [10, 10, -9, -11, -9, -11, -7], [-12, 11, -8, -10, -8, -10, -6], [-11, 12, -7, -9, -7, -9, -5]]
        """
        # Add zeros until the shape has length s
        weight_list = list(weight) # Make sure we have a list
        while len(weight_list) != self._s:
            weight_list.append(0)

        tableau = []
        i = 0
        # Step 0 - Fill first columns of height r
        while i < self._s and weight_list[i] == self._r:
            tableau.append( [self._r - j for j in range(self._r)] )
            i += 1

        # Step 1 - Add the alternating columns until we hit an odd number of columns
        c = -1
        while i < self._s:
            # If it is an odd number of columns
            if i == self._s - 1 or weight_list[i] != weight_list[i+1]:
                c = weight_list[i]
                i += 1
                break
            temp_list = [-(weight_list[i] + j + 1) for j in range(self._r - weight_list[i])]
            for j in range(weight_list[i]):
                temp_list.append(weight_list[i] - j)
            tableau.append(temp_list)
            tableau.append( [self._r - j for j in range(self._r)] )
            i += 2

        # Step 2 - Add the x dependent columns
        x = c + 1
        while i < self._s:
            temp_list = [-x - j for j in range(self._r - x + 1)] # +1 for indexing
            for j in range(x - weight_list[i] - 1): # +1 for indexing
                temp_list.append(self._r - j)
            x = temp_list[-1] # This is the h+1 entry of the column
            for j in range(weight_list[i]):
                temp_list.append(weight_list[i] - j)

            tableau.append(temp_list)
            i += 1

        # Step 3 - Add the final column
        if c > -1:
            val = (self._r + x - 1) // 2
            temp_list = [-x - j for j in range(self._r - val)]
            for j in range(val):
                temp_list.append(val - j)
            tableau.append(temp_list)

        return self.element_class(self, [self.letters(x) for x in flatten(tableau)])

    def _build_module_generators(self):
        """
        Build the module generators.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D',4,1], 2, 3, model='KR')
            sage: KRT._build_module_generators()
            ([[-2, 1, 1], [-1, 2, -1]], [[1, -2, 1], [2, -1, 2]],
             [[1, 1, 1], [2, 2, -1]], [[1, 1, 1], [2, 2, 2]])
        """
        return tuple(self._fill(weight) for weight in
                     horizontal_dominoes_removed(self._s, self._r))

    def from_kirillov_reshetikhin_crystal(self, krc):
        """
        Construct an element of ``self`` from the Kirillov-Reshetikhin
        crystal element ``krc``.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D',4,1], 2, 3, model='KR')
            sage: C = crystals.KirillovReshetikhin(['D',4,1], 2, 3, model='KN')
            sage: krc = C(4,3); krc
            [[3], [4]]
            sage: KRT.from_kirillov_reshetikhin_crystal(krc)
            [[3, -2, 1], [4, -1, 2]]
        """
        # To build a KR tableau from a KR crystal:
        # 1 - start with a highest weight KR tableau generated from the
        #  shape of the KR crystal
        # 2 - determine a path from the KR crystal to its highest weight
        # 3 - apply the inverse path to the highest weight KR tableau
        lifted = krc.lift()
        weight = lifted.to_tableau().shape().conjugate()
        f_str = reversed(lifted.to_highest_weight()[1])
        return self._fill(weight).f_string(f_str)


class KRTableauxTypeHorizonal(KirillovReshetikhinTableaux):
    r"""
    Kirillov-Reshetikhin tableaux `B^{r,s}` of type:

    - `C_n^{(1)}` for `1 \leq r < n`,
    - `A_{2n}^{(2)\dagger}` for `1 \leq r \leq n`.

    TESTS::

        sage: KRT = crystals.KirillovReshetikhin(['C', 3, 1], 2, 2, model='KR')
        sage: TestSuite(KRT).run() # long time
        sage: KRT = crystals.KirillovReshetikhin(CartanType(['A', 4, 2]).dual(), 2, 2, model='KR')
        sage: TestSuite(KRT).run()
    """

    def _fill(self, shape):
        r"""
        Return the highest weight KR tableau of weight ``shape``.

        INPUT:

        - ``shape`` -- the shape of the KR crystal's tableau

        OUTPUT: a `r \times s` tableau

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['C', 5, 1], 3, 5, model='KR')
            sage: KRT._fill([3,3,1])
            [[1, 1, 1, -3, 1], [2, 2, 2, -2, 2], [3, -3, 3, -1, 3]]
            sage: KRT = crystals.KirillovReshetikhin(['C', 10, 1], 5, 6, model='KR')
            sage: KRT._fill([6,4,2,2])
            [[1, 1, 1, 1, 1, 1], [2, 2, 2, 2, -5, 2], [3, 3, -5, 3, -4, 3], [4, 4, -4, 4, -3, 4], [-5, 5, -3, 5, -2, 5]]
            sage: KRT._fill([6,4])
            [[1, 1, 1, 1, 1, 1], [2, 2, 2, 2, -5, 2], [-5, 3, -5, 3, -4, 3], [-4, 4, -4, 4, -3, 4], [-3, 5, -3, 5, -2, 5]]
        """
        # Add zeros until the shape has length s
        shape_list = list(shape) # Make sure we have a list
        while len(shape_list) != self._r:
            shape_list.append(0)

        lst = []
        for col in range(1, self._s+1):
            if (self._s - col) % 2 == 0:
                lst.extend( [self.letters(self._r - x) for x in range(self._r)] )
            else:
                m = self._r
                for j, val in enumerate(shape_list):
                    if col >= val:
                        m = j
                        break
                lst.extend([self.letters(-x) for x in range(m+1, self._r+1)])
                lst.extend([self.letters(m - x) for x in range(m)])

        return self.element_class(self, lst)

    def _build_module_generators(self):
        """
        Build the module generators.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['C',4,1], 2, 3, model='KR')
            sage: KRT._build_module_generators()
            ([[1, -2, 1], [2, -1, 2]], [[1, 1, 1], [2, -2, 2]], [[1, 1, 1], [2, 2, 2]])
        """
        return tuple(self._fill(shape) for shape in horizontal_dominoes_removed(self._r, self._s))

    def from_kirillov_reshetikhin_crystal(self, krc):
        """
        Construct an element of ``self`` from the Kirillov-Reshetikhin
        crystal element ``krc``.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['C',4,1], 2, 3, model='KR')
            sage: C = crystals.KirillovReshetikhin(['C',4,1], 2, 3, model='KN')
            sage: krc = C(4,3); krc
            [[3], [4]]
            sage: KRT.from_kirillov_reshetikhin_crystal(krc)
            [[3, -2, 1], [4, -1, 2]]
        """
        # To build a KR tableau from a KR crystal:
        # 1 - start with a highest weight KR tableau generated from the
        #  shape of the KR crystal
        # 2 - determine a path from the KR crystal to its highest weight
        # 3 - apply the inverse path to the highest weight KR tableau
        lifted = krc.lift()
        shape = lifted.to_tableau().shape()
        f_str = reversed(lifted.to_highest_weight()[1])
        return self._fill(shape).f_string(f_str)


class KRTableauxTypeBox(KRTableauxTypeVertical):
    r"""
    Kirillov-Reshetikhin tableaux `B^{r,s}` of type:

    - `A_{2n}^{(2)}` for all `r \leq n`,
    - `D_{n+1}^{(2)}` for all `r < n`,
    - `D_4^{(3)}` for `r = 1`.

    TESTS::

        sage: KRT = crystals.KirillovReshetikhin(['A', 4, 2], 2, 2, model='KR')
        sage: TestSuite(KRT).run()
        sage: KRT = crystals.KirillovReshetikhin(['D', 4, 2], 2, 2, model='KR')
        sage: TestSuite(KRT).run() # long time
        sage: KRT = crystals.KirillovReshetikhin(['D', 4, 3], 1, 2, model='KR')
        sage: TestSuite(KRT).run() # long time
    """

    def _fill(self, weight):
        r"""
        Return the highest weight KR tableau of weight ``weight``.

        INPUT:

        - ``weight`` -- the weight of the highest weight KR tableau (the
          conjugate of the shape of the KR crystal's tableau)

        OUTPUT: a `r \times s` tableau

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D', 4, 1], 2, 1, model='KR')
            sage: KRT._fill([])
            [[1], [-1]]
            sage: KRT = crystals.KirillovReshetikhin(['D', 14, 1], 12, 7, model='KR')
            sage: KRT._fill([10,10,8,2,2,2])
            [[1, 1, 1, 1, 1, 7, 1], [2, 2, 2, 2, 2, 8, 2], [3, 3, 7, 9, 7, 9, 3], [4, 4, 8, 10, 8, 10, 4], [5, 5, 9, 11, 9, 11, 5], [6, 6, 10, 12, 10, 12, 6], [7, 7, 11, -12, 11, -12, 7], [8, 8, 12, -11, 12, -11, 8], [9, 9, -12, -10, -12, -10, 9], [10, 10, -11, -9, -11, -9, -9], [-12, 11, -10, -8, -10, -8, -8], [-11, 12, -9, -7, -9, -7, -7]]
            sage: KRT._fill([10,10,6,2,2,2])
            [[1, 1, 1, 1, 1, 5, 1], [2, 2, 2, 2, 2, 6, 2], [3, 3, 9, 7, 9, 7, 3], [4, 4, 10, 8, 10, 8, 4], [5, 5, 11, 9, 11, 9, 5], [6, 6, 12, 10, 12, 10, 6], [7, 7, -12, 11, -12, 11, 7], [8, 8, -11, 12, -11, 12, 8], [9, 9, -10, -12, -10, -12, -8], [10, 10, -9, -11, -9, -11, -7], [-12, 11, -8, -10, -8, -10, -6], [-11, 12, -7, -9, -7, -9, -5]]
        """
        # Add zeros until the shape has length s
        weight_list = list(weight) # Make sure we have a list
        while len(weight_list) != self._s:
            weight_list.append(0)

        tableau = []
        i = 0
        # Step 0 - Fill first columns of height r
        while i < self._s and weight_list[i] == self._r:
            tableau.append( [self._r - j for j in range(self._r)] )
            i += 1

        # Step 1 - Add the alternating columns until we hit an odd number of columns
        c = -1
        while i < self._s:
            # If it is an odd number of columns
            if i == self._s - 1 or weight_list[i] != weight_list[i+1]:
                c = weight_list[i]
                i += 1
                break
            temp_list = [-(weight_list[i] + j + 1) for j in range(self._r - weight_list[i])]
            for j in range(weight_list[i]):
                temp_list.append(weight_list[i] - j)
            tableau.append(temp_list)
            tableau.append( [self._r - j for j in range(self._r)] )
            i += 2

        # Step 2 - Add the x dependent columns
        x = c + 1
        while i < self._s:
            temp_list = [-x - j for j in range(self._r - x + 1)] # +1 for indexing
            for j in range(x - weight_list[i] - 1): # +1 for indexing
                temp_list.append(self._r - j)
            x = temp_list[-1] # This is the h+1 entry of the column
            for j in range(weight_list[i]):
                temp_list.append(weight_list[i] - j)

            tableau.append(temp_list)
            i += 1

        # Step 3 - Add the final column
        if c > -1:
            val = x - 1
            temp_list = ['E' for j in range(self._r - val)]
            for j in range(val):
                temp_list.append(val - j)
            tableau.append(temp_list)

        return self.element_class(self, [self.letters(x) for x in flatten(tableau)])

    def _build_module_generators(self):
        """
        Build the module generators.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['A',4,2], 2, 2, model='KR')
            sage: KRT._build_module_generators()
            ([[-2, 1], [-1, 2]], [[2, 1], [-2, E]], [[1, E], [2, E]],
             [[1, 1], [-2, 2]], [[1, 1], [2, E]], [[1, 1], [2, 2]])
        """
        return tuple(self._fill(weight) for weight in partitions_in_box(self._s, self._r))


class KRTableauxSpin(KRTableauxRectangle):
    r"""
    Kirillov-Reshetikhin tableaux `B^{r,s}` of type `D_n^{(1)}` with
    `r = n, n-1`.

    TESTS::

        sage: KRT = crystals.KirillovReshetikhin(['D', 4, 1], 3, 2, model='KR')
        sage: TestSuite(KRT).run()
        sage: KRT = crystals.KirillovReshetikhin(['D', 4, 1], 4, 2, model='KR')
        sage: TestSuite(KRT).run()
    """

    def _build_module_generators(self):
        r"""
        Build the module generators.

        There is only one module generator which corresponds to a single
        `n \times s` rectangle.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D', 4, 1], 3, 3, model='KR')
            sage: KRT._build_module_generators()
            ([[1, 1, 1], [2, 2, 2], [3, 3, 3], [-4, -4, -4]],)
            sage: KRT = crystals.KirillovReshetikhin(['D', 4, 1], 4, 3, model='KR')
            sage: KRT._build_module_generators()
            ([[1, 1, 1], [2, 2, 2], [3, 3, 3], [4, 4, 4]],)
        """
        n = self.cartan_type().classical().rank()
        if self._r == n:
            return KRTableauxRectangle._build_module_generators(self)

        tableau = []
        for i in range(self._s):
            tableau.append( [-n] + [self._r - j for j in range(self._r)] )

        return (self.element_class(self, [self.letters(x) for x in flatten(tableau)]),)


class KRTableauxBn(KRTableauxTypeHorizonal):
    """
    Kirillov-Reshetkhin tableaux `B^{n,s}` of type `B_n^{(1)}`.

    TESTS::

        sage: KRT = crystals.KirillovReshetikhin(['B', 2, 1], 2, 3, model='KR')
        sage: TestSuite(KRT).run()
    """

    def _build_module_generators(self):
        """
        Build the module generators.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['B', 2, 1], 2, 2, model='KR')
            sage: KRT._build_module_generators()
            ([[-2, 1], [-1, 2]], [[1, 1], [2, 2]])
        """
        odd = int(self._s % 2)
        shapes = [[int(x * 2 + odd) for x in sh] for sh
                   in vertical_dominoes_removed(self._r, self._s // 2)]
        return tuple(self._fill(sh) for sh in shapes)

    def from_kirillov_reshetikhin_crystal(self, krc):
        """
        Construct an element of ``self`` from the Kirillov-Reshetikhin
        crystal element ``krc``.

        EXAMPLES::

            sage: KR = crystals.KirillovReshetikhin(['B',3,1], 3, 3, model='KR')
            sage: C = crystals.KirillovReshetikhin(['B',3,1], 3, 3, model='KN')
            sage: krc = C.module_generators[1].f_string([3,2,3,1,3,3]); krc
            [++-, [[2], [0], [-3]]]
            sage: KR.from_kirillov_reshetikhin_crystal(krc)
            [[1, 1, 2], [2, 2, -3], [-3, -3, -1]]
        """
        # To build a KR tableau from a type B_n spinor KR crystal:
        # 1 - determine a path from the KR crystal to its highest weight
        # 2 - find the corresponding highest weight KR tableau
        # 3 - apply the inverse path to the highest weight KR tableau
        lifted = krc.lift()
        to_hw = lifted.to_highest_weight()
        f_str = reversed(to_hw[1])
        wt = to_hw[0].weight()
        for x in self.module_generators:
            if x.classical_weight() == wt:
                return x.f_string(f_str)
        raise ValueError("no matching highest weight element found")


class KirillovReshetikhinTableauxElement(TensorProductOfRegularCrystalsElement):
    r"""
    A Kirillov-Reshetikhin tableau.

    For more information, see
    :class:`~sage.combinat.rigged_configurations.kr_tableaux.KirillovReshetikhinTableaux`
    and
    :class:`~sage.combinat.rigged_configurations.tensor_product_kr_tableaux.TensorProductOfKirillovReshetikhinTableaux`.
    """

    def __init__(self, parent, list, **options):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['A', 4, 1], 2, 1, model='KR')
            sage: elt = KRT(4, 3); elt
            [[3], [4]]
            sage: TestSuite(elt).run()
        """
        # Make sure we are a list of letters
        if list and not isinstance(list[0], (parent.letters.element_class, EmptyLetter)):
            list = [parent.letters(x) for x in list]
        TensorProductOfRegularCrystalsElement.__init__(self, parent, list)

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['A', 4, 1], 2, 1, model='KR')
            sage: KRT(3,2)
            [[2], [3]]
        """
        return repr(self.to_array())

    def _repr_diagram(self):
        """
        Return a string representation of ``self`` as a diagram.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['A',4,1], 2, 2, model='KR')
            sage: elt = KRT(2,1,4,3)
            sage: print(elt._repr_diagram())
              1  3
              2  4
        """
        return self.to_tableau()._repr_diagram()

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['A', 4, 1], 2, 3, model='KR')
            sage: latex(KRT(3,2,4,2,4,3))
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{3}c}\cline{1-3}
            \lr{2}&\lr{2}&\lr{3}\\\cline{1-3}
            \lr{3}&\lr{4}&\lr{4}\\\cline{1-3}
            \end{array}$}
            }
        """
        from sage.combinat.output import tex_from_array
        return tex_from_array([[val._latex_() for val in row] for row in self.to_array()])

    def _ascii_art_(self):
        r"""
        Return an ASCII art representation of ``self``.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['A',4,1], 2, 2, model='KR')
            sage: ascii_art(KRT(2,1,4,3))
              1  3
              2  4
        """
        from sage.typeset.ascii_art import AsciiArt
        return AsciiArt(self._repr_diagram().splitlines())

    def _unicode_art_(self):
        r"""
        Return a unicode art representation of ``self``.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D',4,1], 2, 2, model='KR')
            sage: unicode_art(KRT(2,1,-4,3))
            ┌───┬───┐
            │ 1 │ 3 │
            ├───┼───┤
            │ 2 │ 4̄ │
            └───┴───┘
        """
        return self.to_tableau()._unicode_art_()

    def to_kirillov_reshetikhin_crystal(self):
        r"""
        Construct a
        :func:`~sage.combinat.crystals.kirillov_reshetihkin.KashiwaraNakashimaTableaux`
        element from ``self``.

        We construct the Kirillov-Reshetikhin crystal element as follows:

        1. Determine the shape `\lambda` of the KR crystal from the weight.
        2. Determine a path `e_{i_1} e_{i_2} \cdots e_{i_k}` to the highest
           weight.
        3. Apply `f_{i_k} \cdots f_{i_2} f_{i_1}` to a highest weight KR
           crystal of shape `\lambda`.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D',4,1], 2, 2, model='KR')
            sage: elt = KRT(3,2,-1,1); elt
            [[2, 1], [3, -1]]
            sage: elt.to_kirillov_reshetikhin_crystal()
            [[2], [3]]

        TESTS:

        Spinor tests::

            sage: KRT = crystals.KirillovReshetikhin(['D',4,1], 4, 3, model='KR')
            sage: KRC = crystals.KirillovReshetikhin(['D',4,1], 4, 3, model='KN')
            sage: elt = KRT(-3,-4,2,1,-3,-4,2,1,-2,-4,3,1); elt
            [[1, 1, 1], [2, 2, 3], [-4, -4, -4], [-3, -3, -2]]
            sage: ret = elt.to_kirillov_reshetikhin_crystal(); ret
            [++--, [[1], [3], [-4], [-3]]]
            sage: test = KRT(ret); test
            [[1, 1, 1], [2, 2, 3], [-4, -4, -4], [-3, -3, -2]]
            sage: test == elt
            True
        """
        return self.parent().kirillov_reshetikhin_crystal()(self)

    @cached_method
    def to_array(self, rows=True):
        r"""
        Return a 2-dimensional array representation of this
        Kirillov-Reshetikhin element.

        If the output is in rows, then it outputs the top row first (in the
        English convention) from left to right.

        For example: if the reading word is `[2, 1, 4, 3]`, so as a
        `2 \times 2` tableau::

            1 3
            2 4

        we output ``[[1, 3], [2, 4]]``.

        If the output is in columns, then it outputs the leftmost column first
        with the bottom element first. In other words this parses the reading
        word into its columns.

        Continuing with the previous example, the output would be
        ``[[2, 1], [4, 3]]``.

        INPUT:

        - ``rows`` -- boolean (default: ``True``); set to ``True`` if the
          resulting array is by row, otherwise it is by column

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['A', 4, 1], 2, 2, model='KR')
            sage: elt = KRT(2, 1, 4, 3)
            sage: elt.to_array()
            [[1, 3], [2, 4]]
            sage: elt.to_array(False)
            [[2, 1], [4, 3]]
        """
        ret_list = []
        h = self.parent()._tableau_height
        s = self.parent()._s
        if rows:
            for i in reversed(range(h)):
                row = []
                for j in range(s):
                    row.append(self[j * h + i])
                ret_list.append(row)
        else:
            for j in range(s):
                col = []
                for i in range(h):
                    col.append(self[j * h + i])
                ret_list.append(col)

        return ret_list

    @cached_method
    def to_tableau(self):
        """
        Return a :class:`Tableau` object of ``self``.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['A', 4, 1], 2, 2, model='KR')
            sage: elt = KRT(2, 1, 4, 3); elt
            [[1, 3], [2, 4]]
            sage: t = elt.to_tableau(); t
            [[1, 3], [2, 4]]
            sage: type(t)
            <class 'sage.combinat.tableau.Tableaux_all_with_category.element_class'>
        """
        return Tableau(self.to_array())

    def pp(self):
        """
        Pretty print ``self``.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['A', 4, 1], 2, 2, model='KR')
            sage: elt = KRT(2, 1, 4, 3); elt
            [[1, 3], [2, 4]]
            sage: elt.pp()
            1  3
            2  4
        """
        self.to_tableau().pp()

    def to_classical_highest_weight(self, index_set=None):
        r"""
        Return the classical highest weight element corresponding to ``self``.

        INPUT:

        - ``index_set`` -- (default: ``None``) return the highest weight
          with respect to the index set; if ``None`` is passed in, then this
          uses the classical index set

        OUTPUT:

        A pair ``[H, f_str]`` where ``H`` is the highest weight element and
        ``f_str`` is a list of `a_i` of `f_{a_i}` needed to reach ``H``.

        EXAMPLES::

            sage: KRTab = crystals.KirillovReshetikhin(['D',4,1], 2, 2, model='KR')
            sage: elt = KRTab(3,2,-1,1); elt
            [[2, 1], [3, -1]]
            sage: elt.to_classical_highest_weight()
            [[[1, 1], [2, -1]], [1, 2]]
        """
        if index_set is None:
            index_set = self.parent()._cartan_type.classical().index_set()
        for i in index_set:
            next = self.e(i)
            if next is not None:
                hw = next.to_classical_highest_weight(index_set=index_set)
                return [hw[0], [i] + hw[1]]
        return [self, []]

    def weight(self):
        """
        Return the weight of ``self``.

        EXAMPLES::

            sage: KR = crystals.KirillovReshetikhin(['D',4,1], 2, 2, model='KR')
            sage: KR.module_generators[1].weight()
            -2*Lambda[0] + Lambda[2]
        """
        return self.Phi() - self.Epsilon()

    @cached_method
    def classical_weight(self):
        r"""
        Return the classical weight of ``self``.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D',4,1], 2, 2, model='KR')
            sage: elt = KRT(3,2,-1,1); elt
            [[2, 1], [3, -1]]
            sage: elt.classical_weight()
            (0, 1, 1, 0)
        """
        F = self.cartan_type().classical().root_system()
        if F.ambient_space() is None:
            WLR = F.weight_lattice()
        else:
            WLR = F.ambient_space()
        return sum((self[j].weight() for j in range(len(self))), WLR.zero())

    def e(self, i):
        """
        Perform the action of `e_i` on ``self``.

        .. TODO::

            Implement a direct action of `e_0` without moving to KR crystals.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D',4,1], 2, 2, model='KR')
            sage: KRT.module_generators[0].e(0)
            [[-2, 1], [-1, -1]]
        """
        if i == self.parent()._cartan_type.special_node():
            ret = self.to_kirillov_reshetikhin_crystal().e0()
            if ret is None:
                return None
            return ret.to_kirillov_reshetikhin_tableau()
        return TensorProductOfRegularCrystalsElement.e(self, i)

    def f(self, i):
        """
        Perform the action of `f_i` on ``self``.

        .. TODO::

            Implement a direct action of `f_0` without moving to KR crystals.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D',4,1], 2, 2, model='KR')
            sage: KRT.module_generators[0].f(0)
            [[1, 1], [2, -1]]
        """
        if i == self.parent()._cartan_type.special_node():
            ret = self.to_kirillov_reshetikhin_crystal().f0()
            if ret is None:
                return None
            return ret.to_kirillov_reshetikhin_tableau()
        return TensorProductOfRegularCrystalsElement.f(self, i)

    def epsilon(self, i):
        r"""
        Compute `\varepsilon_i` of ``self``.

        .. TODO::

            Implement a direct action of `\varepsilon_0` without moving to
            KR crystals.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D',4,1], 2, 2, model='KR')
            sage: KRT.module_generators[0].epsilon(0)
            2
        """
        if i == self.parent()._cartan_type.special_node():
            return self.to_kirillov_reshetikhin_crystal().epsilon0()
        return TensorProductOfRegularCrystalsElement.epsilon(self, i)

    def phi(self, i):
        r"""
        Compute `\varphi_i` of ``self``.

        .. TODO::

            Compute `\varphi_0` without moving to KR crystals.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D',4,1], 2, 2, model='KR')
            sage: KRT.module_generators[0].phi(0)
            2
        """
        if i == self.parent()._cartan_type.special_node():
            return self.to_kirillov_reshetikhin_crystal().phi0()
        return TensorProductOfRegularCrystalsElement.phi(self, i)

    def left_split(self):
        r"""
        Return the image of ``self`` under the left column splitting map.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D',4,1], 2, 3, model='KR')
            sage: mg = KRT.module_generators[1]; mg.pp()
              1 -2  1
              2 -1  2
            sage: ls = mg.left_split(); ls.pp()
              1 (X)  -2  1
              2      -1  2
            sage: ls.parent()
            Tensor product of Kirillov-Reshetikhin tableaux of type ['D', 4, 1] and factor(s) ((2, 1), (2, 2))
        """
        P = self.parent()
        if P._s == 1:
            raise ValueError("cannot split a single column")
        from sage.combinat.rigged_configurations.tensor_product_kr_tableaux import \
                TensorProductOfKirillovReshetikhinTableaux
        r = P._r
        TP = TensorProductOfKirillovReshetikhinTableaux(P._cartan_type, [[r, 1], [r, P._s-1]])
        lf = TP.crystals[0](*(self[:r]))
        rf = TP.crystals[1](*(self[r:]))
        return TP(lf, rf)

    def right_split(self):
        r"""
        Return the image of ``self`` under the right column splitting map.

        Let `\ast` denote the :meth:`Lusztig involution<lusztig_involution>`,
        and `\mathrm{ls}` as the :meth:`left splitting map<left_split>`.
        The right splitting map is defined as
        `\mathrm{rs} := \ast \circ \mathrm{ls} \circ \ast`.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D',4,1], 2, 3, model='KR')
            sage: mg = KRT.module_generators[1]; mg.pp()
              1 -2  1
              2 -1  2
            sage: ls = mg.right_split(); ls.pp()
             -2  1 (X)   1
             -1  2       2
            sage: ls.parent()
            Tensor product of Kirillov-Reshetikhin tableaux of type ['D', 4, 1] and factor(s) ((2, 2), (2, 1))
        """
        return self.lusztig_involution().left_split().lusztig_involution()


KirillovReshetikhinTableaux.Element = KirillovReshetikhinTableauxElement


class KRTableauxSpinElement(KirillovReshetikhinTableauxElement):
    r"""
    Kirillov-Reshetikhin tableau for spinors.

    Here we are in the embedding `B(\Lambda_n) \hookrightarrow
    B(2 \Lambda_n)`, so `e_i` and `f_i` act by `e_i^2` and `f_i^2`
    respectively for all `i \neq 0`. We do this so our columns are full
    width (as opposed to half width and/or uses a `\pm` representation).
    """

    def e(self, i):
        r"""
        Calculate the action of `e_i` on ``self``.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D',4,1], 4, 1, model='KR')
            sage: KRT(-1, -4, 3, 2).e(1)
            [[1], [3], [-4], [-2]]
            sage: KRT(-1, -4, 3, 2).e(3)
        """
        if i == self.parent()._cartan_type.special_node():
            # Only need to do it once since we pull to the KR crystal
            return KirillovReshetikhinTableauxElement.e(self, i)

        half = KirillovReshetikhinTableauxElement.e(self, i)
        if half is None:
            return None
        return KirillovReshetikhinTableauxElement.e(half, i)

    def f(self, i):
        r"""
        Calculate the action of `f_i` on ``self``.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D',4,1], 4, 1, model='KR')
            sage: KRT(-1, -4, 3, 2).f(1)
            sage: KRT(-1, -4, 3, 2).f(3)
            [[2], [4], [-3], [-1]]
        """
        if i == self.parent()._cartan_type.special_node():
            # Only need to do it once since we pull to the KR crystal
            return KirillovReshetikhinTableauxElement.f(self, i)

        half = KirillovReshetikhinTableauxElement.f(self, i)
        if half is None:
            return None

        return KirillovReshetikhinTableauxElement.f(half, i)

    def epsilon(self, i):
        r"""
        Compute `\varepsilon_i` of ``self``.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D',4,1], 4, 1, model='KR')
            sage: KRT(-1, -4, 3, 2).epsilon(1)
            1
            sage: KRT(-1, -4, 3, 2).epsilon(3)
            0
        """
        if i == self.parent()._cartan_type.special_node():
            # Don't need to half it since we pull to the KR crystal
            return KirillovReshetikhinTableauxElement.epsilon(self, i)
        return KirillovReshetikhinTableauxElement.epsilon(self, i) // 2

    def phi(self, i):
        r"""
        Compute `\varphi_i` of ``self``.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D',4,1], 4, 1, model='KR')
            sage: KRT(-1, -4, 3, 2).phi(1)
            0
            sage: KRT(-1, -4, 3, 2).phi(3)
            1
        """
        if i == self.parent()._cartan_type.special_node():
            # Don't need to half it since we pull to the KR crystal
            return KirillovReshetikhinTableauxElement.phi(self, i)
        return KirillovReshetikhinTableauxElement.phi(self, i) // 2

    @cached_method
    def to_array(self, rows=True):
        r"""
        Return a 2-dimensional array representation of this
        Kirillov-Reshetikhin element.

        If the output is in rows, then it outputs the top row first (in the
        English convention) from left to right.

        For example: if the reading word is `[2, 1, 4, 3]`, so as a
        `2 \times 2` tableau::

            1 3
            2 4

        we output ``[[1, 3], [2, 4]]``.

        If the output is in columns, then it outputs the leftmost column first
        with the bottom element first. In other words this parses the reading
        word into its columns.

        Continuing with the previous example, the output would be
        ``[[2, 1], [4, 3]]``.

        INPUT:

        - ``rows`` -- boolean (default: ``True``); set to ``True`` if the
          resulting array is by row, otherwise it is by column

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D', 4, 1], 4, 3, model='KR')
            sage: elt = KRT(-3,-4,2,1,-3,-4,2,1,-2,-4,3,1)
            sage: elt.to_array()
            [[1, 1, 1], [2, 2, 3], [-4, -4, -4], [-3, -3, -2]]
            sage: elt.to_array(False)
            [[-3, -4, 2, 1], [-3, -4, 2, 1], [-2, -4, 3, 1]]
        """
        ret_list = []
        h = self.parent()._cartan_type.classical().rank()
        s = self.parent()._s
        if rows:
            for i in reversed(range(h)):
                row = []
                for j in range(s):
                    row.append(self[j * h + i])
                ret_list.append(row)
        else:
            for j in range(s):
                col = []
                for i in range(h):
                    col.append(self[j * h + i])
                ret_list.append(col)

        return ret_list

    def left_split(self):
        """
        Return the image of ``self`` under the left column splitting map.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D', 4, 1], 4, 3, model='KR')
            sage: elt = KRT(-3,-4,2,1,-3,-4,2,1,-2,-4,3,1); elt.pp()
              1  1  1
              2  2  3
             -4 -4 -4
             -3 -3 -2
            sage: elt.left_split().pp()
              1 (X)   1  1
              2       2  3
             -4      -4 -4
             -3      -3 -2
        """
        P = self.parent()
        if P._s == 1:
            raise ValueError("cannot split a single column")
        from sage.combinat.rigged_configurations.tensor_product_kr_tableaux import \
                TensorProductOfKirillovReshetikhinTableaux
        h = P._cartan_type.classical().rank()
        TP = TensorProductOfKirillovReshetikhinTableaux(P._cartan_type, [[P._r, 1], [P._r, P._s-1]])
        lf = TP.crystals[0](*(self[:h]))
        rf = TP.crystals[1](*(self[h:]))
        return TP(lf, rf)

    # FIXME: This is a copy of the above classical weight, and cached_method
    #   overwrites this method if it is called via super.
    @cached_method
    def classical_weight(self):
        r"""
        Return the classical weight of ``self``.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D', 4, 1], 4, 1, model='KR')
            sage: KRT.module_generators[0].classical_weight()
            (1/2, 1/2, 1/2, 1/2)
        """
        F = self.cartan_type().classical().root_system()
        if F.ambient_space() is None:
            WLR = F.weight_lattice()
        else:
            WLR = F.ambient_space()
        return sum((self[j].weight() for j in range(len(self))), WLR.zero()) / 2


KRTableauxBn.Element = KRTableauxSpinElement
KRTableauxSpin.Element = KRTableauxSpinElement


class KRTableauxDTwistedSpin(KRTableauxRectangle):
    r"""
    Kirillov-Reshetikhin tableaux `B^{r,s}` of type `D_n^{(2)}` with `r = n`.

    EXAMPLES::

        sage: KRT = crystals.KirillovReshetikhin(['D', 4, 2], 1, 1, model='KR')
        sage: KRT.cardinality()
        8
        sage: KRC = crystals.KirillovReshetikhin(['D', 4, 2], 1, 1, model='KN')
        sage: KRT.cardinality() == KRC.cardinality()
        True
    """
    Element = KRTableauxSpinElement


class KRTableauxTypeFromRCElement(KirillovReshetikhinTableauxElement):
    r"""
    A Kirillov-Reshetikhin tableau constructed from rigged configurations
    under the bijection `\Phi`.
    """

    def e(self, i):
        """
        Perform the action of `e_i` on ``self``.

        .. TODO::

            Implement a direct action of `e_0` without moving to
            rigged configurations.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D',4,3], 2, 1, model='KR')
            sage: KRT.module_generators[0].e(0)
            [[2], [E]]
        """
        if i == self.parent().cartan_type().special_node():
            P = self.parent()
            from sage.combinat.rigged_configurations.tensor_product_kr_tableaux import TensorProductOfKirillovReshetikhinTableaux
            K = TensorProductOfKirillovReshetikhinTableaux(P.cartan_type(), [[2, P.s()]])
            ret = K(self).to_rigged_configuration()
            RC = ret.parent()
            ret = ret.to_virtual_configuration().e(0)
            if ret is None:
                return None
            ret = RC.from_virtual(ret)
            return ret.to_tensor_product_of_kirillov_reshetikhin_tableaux()[0]
        return TensorProductOfRegularCrystalsElement.e(self, i)

    def f(self, i):
        """
        Perform the action of `f_i` on ``self``.

        .. TODO::

            Implement a direct action of `f_0` without moving to
            rigged configurations.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D',4,3], 2, 1, model='KR')
            sage: KRT.module_generators[0].f(0)
            sage: KRT.module_generators[3].f(0)
            [[1], [0]]
        """
        if i == self.parent().cartan_type().special_node():
            P = self.parent()
            from sage.combinat.rigged_configurations.tensor_product_kr_tableaux import TensorProductOfKirillovReshetikhinTableaux
            K = TensorProductOfKirillovReshetikhinTableaux(P.cartan_type(), [[2, P.s()]])
            ret = K(self).to_rigged_configuration()
            RC = ret.parent()
            ret = ret.to_virtual_configuration().f(0)
            if ret is None:
                return None
            ret = RC.from_virtual(ret)
            return ret.to_tensor_product_of_kirillov_reshetikhin_tableaux()[0]
        return TensorProductOfRegularCrystalsElement.f(self, i)

    def epsilon(self, i):
        r"""
        Compute `\varepsilon_i` of ``self``.

        .. TODO::

            Implement a direct action of `\epsilon_0` without moving to
            KR crystals.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D',4,3], 2, 2, model='KR')
            sage: KRT.module_generators[0].epsilon(0)
            6
        """
        if i == self.parent().cartan_type().special_node():
            P = self.parent()
            from sage.combinat.rigged_configurations.tensor_product_kr_tableaux import TensorProductOfKirillovReshetikhinTableaux
            K = TensorProductOfKirillovReshetikhinTableaux(P.cartan_type(), [[2, P.s()]])
            rc = K(self).to_rigged_configuration().to_virtual_configuration()
            return rc.epsilon(0)
        return TensorProductOfRegularCrystalsElement.epsilon(self, i)

    def phi(self, i):
        r"""
        Compute `\varphi_i` of ``self``.

        .. TODO::

            Compute `\phi_0` without moving to KR crystals.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D',4,3], 2, 2, model='KR')
            sage: KRT.module_generators[0].phi(0)
            0
        """
        if i == self.parent().cartan_type().special_node():
            P = self.parent()
            from sage.combinat.rigged_configurations.tensor_product_kr_tableaux import TensorProductOfKirillovReshetikhinTableaux
            K = TensorProductOfKirillovReshetikhinTableaux(P.cartan_type(), [[2, P.s()]])
            rc = K(self).to_rigged_configuration().to_virtual_configuration()
            return rc.phi(0)
        return TensorProductOfRegularCrystalsElement.phi(self, i)


class KRTableauxTypeFromRC(KirillovReshetikhinTableaux):
    r"""
    Kirillov-Reshetikhin tableaux `B^{r,s}` constructed from rigged
    configurations under the bijection `\Phi`.

    .. WARNING::

        The Kashiwara-Nakashima version is not implemented due to the
        non-trivial multiplicities of classical components, so
        :meth:`classical_decomposition` does not work.
    """

    def __init__(self, cartan_type, r, s):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D', 4, 3], 2, 1, model='KR')
            sage: TestSuite(KRT).run() # long time
        """
        # We must modify the constructor of KirillovReshetikhinTableaux
        self._r = r
        self._s = s
        self._cartan_type = cartan_type
        Parent.__init__(self, category=KirillovReshetikhinCrystals())
        self.letters = CrystalOfLetters(cartan_type.classical())

    @lazy_attribute
    def module_generators(self):
        """
        The module generators of ``self``.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D',4,3], 2, 1, model='KR')
            sage: KRT.module_generators
            ([[1], [2]], [[1], [0]], [[1], [E]], [[E], [E]])
        """
        return self._build_module_generators()

    def _build_module_generators(self):
        r"""
        Return the module generators of ``self``.

        EXAMPLES::

            sage: KRT = crystals.KirillovReshetikhin(['D',4,3], 2, 1, model='KR')
            sage: KRT._build_module_generators()
            ([[1], [2]], [[1], [0]], [[1], [E]], [[E], [E]])
        """
        from sage.combinat.rigged_configurations.rigged_configurations import RiggedConfigurations
        RC = RiggedConfigurations(self._cartan_type, [[self._r, self._s]])
        return tuple(mg.to_tensor_product_of_kirillov_reshetikhin_tableaux()[0]
                     for mg in RC.module_generators)

    @lazy_attribute
    def _tableau_height(self):
        """
        The height of the tableaux in ``self``.

        EXAMPLES::

            sage: ct = CartanType(['E',6,1])
            sage: [crystals.KirillovReshetikhin(ct, r, 1, model='KR')._tableau_height
            ....:  for r in ct.classical().index_set()]
            [1, 3, 2, 3, 4, 2]
        """
        if self._cartan_type.type() == 'E':
            if self._cartan_type.classical().rank() == 6:
                #       6   2 - 5
                #      /   /
                # 0 - 1 - 3 - 4
                if self._r == 1:
                    return 1
                if self._r == [3, 6]:
                    return 2
                if self._r in [2, 4]:
                    return 3
                if self._r == 5:
                    return 4
            if self._cartan_type.classical().rank() == 7:
                #     1-2-3
                #    /
                # 0-7-6-5-4
                if self._r <= 3:
                    return self._r + 1
                return 8 - self._r
            if self._cartan_type.classical().rank() == 8:
                #     1-2-3
                #    /
                # 0-8-7-6-5-4
                if self._r <= 3:
                    return self._r + 1
                return 9 - self._r
        if not self._cartan_type.is_untwisted_affine():
            if self._cartan_type.dual().type() == 'G':
                return self._r
        return len(self.module_generators[0]) // self._s

    Element = KRTableauxTypeFromRCElement
