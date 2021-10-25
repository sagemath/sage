r"""
Littelmann paths

AUTHORS:

- Mark Shimozono, Anne Schilling (2012): Initial version
- Anne Schilling (2013): Implemented
  :class:`~sage.combinat.crystals.littelmann_path.CrystalOfProjectedLevelZeroLSPaths`
- Travis Scrimshaw (2016): Implemented
  :class:`~sage.combinat.crystals.littelmann_path.InfinityCrystalOfLSPaths`
"""
#****************************************************************************
#       Copyright (C) 2012 Mark Shimozono
#                          Anne Schilling
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
#****************************************************************************

from sage.misc.cachefunc import cached_in_parent_method, cached_method
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.parent import Parent
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.categories.regular_crystals import RegularCrystals
from sage.categories.classical_crystals import ClassicalCrystals
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.loop_crystals import (RegularLoopCrystals,
                                           KirillovReshetikhinCrystals)
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.weyl_group import WeylGroup
from sage.rings.integer import Integer
from sage.rings.rational_field import QQ
from sage.combinat.root_system.root_system import RootSystem
from sage.arith.misc import integer_floor as floor
from sage.misc.latex import latex


class CrystalOfLSPaths(UniqueRepresentation, Parent):
    r"""
    Crystal graph of LS paths generated from the straight-line path to a given weight.

    INPUT:

    - ``cartan_type`` -- (optional) the Cartan type of a finite or affine root system
    - ``starting_weight`` -- a weight; if ``cartan_type`` is given, then the weight should
      be given as a list of coefficients of the fundamental weights, otherwise it should
      be given in the ``weight_space`` basis; for affine highest weight crystals, one needs
      to use the extended weight space.

    The crystal class of piecewise linear paths in the weight space,
    generated from a straight-line path from the origin to a given
    element of the weight lattice.

    OUTPUT:

    - a tuple of weights defining the directions of the piecewise linear segments

    EXAMPLES::

        sage: R = RootSystem(['A',2,1])
        sage: La = R.weight_space(extended = True).basis()
        sage: B = crystals.LSPaths(La[2]-La[0]); B
        The crystal of LS paths of type ['A', 2, 1] and weight -Lambda[0] + Lambda[2]

        sage: C = crystals.LSPaths(['A',2,1],[-1,0,1]); C
        The crystal of LS paths of type ['A', 2, 1] and weight -Lambda[0] + Lambda[2]
        sage: B == C
        True
        sage: c = C.module_generators[0]; c
        (-Lambda[0] + Lambda[2],)
        sage: [c.f(i) for i in C.index_set()]
        [None, None, (Lambda[1] - Lambda[2],)]

        sage: R = C.R; R
        Root system of type ['A', 2, 1]
        sage: Lambda = R.weight_space().basis(); Lambda
        Finite family {0: Lambda[0], 1: Lambda[1], 2: Lambda[2]}
        sage: b=C(tuple([-Lambda[0]+Lambda[2]]))
        sage: b==c
        True
        sage: b.f(2)
        (Lambda[1] - Lambda[2],)

    For classical highest weight crystals we can also compare the results with the tableaux implementation::

        sage: C = crystals.LSPaths(['A',2],[1,1])
        sage: sorted(C, key=str)
        [(-2*Lambda[1] + Lambda[2],), (-Lambda[1] + 1/2*Lambda[2], Lambda[1] - 1/2*Lambda[2]),
         (-Lambda[1] + 2*Lambda[2],), (-Lambda[1] - Lambda[2],),
         (1/2*Lambda[1] - Lambda[2], -1/2*Lambda[1] + Lambda[2]), (2*Lambda[1] - Lambda[2],),
         (Lambda[1] + Lambda[2],), (Lambda[1] - 2*Lambda[2],)]
        sage: C.cardinality()
        8
        sage: B = crystals.Tableaux(['A',2],shape=[2,1])
        sage: B.cardinality()
        8
        sage: B.digraph().is_isomorphic(C.digraph())
        True

    Make sure you use the weight space and not the weight lattice for your weights::

        sage: R = RootSystem(['A',2,1])
        sage: La = R.weight_lattice(extended = True).basis()
        sage: B = crystals.LSPaths(La[2]); B
        Traceback (most recent call last):
        ...
        ValueError: Please use the weight space, rather than weight lattice for your weights

    REFERENCES:

    .. [Li1995b]_
    """

    @staticmethod
    def __classcall_private__(cls, starting_weight, cartan_type = None, starting_weight_parent = None):
        """
        Classcall to mend the input.

        Internally, the
        :class:`~sage.combinat.crystals.littelmann_path.CrystalOfLSPaths` code
        works with a ``starting_weight`` that is in the weight space associated
        to the crystal. The user can, however, also input a ``cartan_type``
        and the coefficients of the fundamental weights as
        ``starting_weight``. This code transforms the input into the right
        format (also necessary for UniqueRepresentation).

        TESTS::

            sage: crystals.LSPaths(['A',2,1],[-1,0,1])
            The crystal of LS paths of type ['A', 2, 1] and weight -Lambda[0] + Lambda[2]

            sage: R = RootSystem(['B',2,1])
            sage: La = R.weight_space(extended=True).basis()
            sage: C = crystals.LSPaths(['B',2,1],[0,0,1])
            sage: B = crystals.LSPaths(La[2])
            sage: B is C
            True
        """
        if cartan_type is not None:
            cartan_type, starting_weight = CartanType(starting_weight), cartan_type
            if cartan_type.is_affine():
                extended = True
            else:
                extended = False

            R = RootSystem(cartan_type)
            P = R.weight_space(extended = extended)
            Lambda = P.basis()
            offset = R.index_set()[Integer(0)]
            starting_weight = P.sum(starting_weight[j-offset]*Lambda[j] for j in R.index_set())
        if starting_weight_parent is None:
            starting_weight_parent = starting_weight.parent()
        else:
            # Both the weight and the parent of the weight are passed as arguments of init to be able
            # to distinguish between crystals with the extended and non-extended weight lattice!
            if starting_weight.parent() != starting_weight_parent:
                raise ValueError("The passed parent is not equal to parent of the inputted weight!")

        return super(CrystalOfLSPaths, cls).__classcall__(cls, starting_weight, starting_weight_parent = starting_weight_parent)

    def __init__(self, starting_weight, starting_weight_parent):
        """
        EXAMPLES::

            sage: C = crystals.LSPaths(['A',2,1],[-1,0,1]); C
            The crystal of LS paths of type ['A', 2, 1] and weight -Lambda[0] + Lambda[2]
            sage: C.R
            Root system of type ['A', 2, 1]
            sage: C.weight
            -Lambda[0] + Lambda[2]
            sage: C.weight.parent()
            Extended weight space over the Rational Field of the Root system of type ['A', 2, 1]
            sage: C.module_generators
            ((-Lambda[0] + Lambda[2],),)

        TESTS::

            sage: C = crystals.LSPaths(['A',2,1], [-1,0,1])
            sage: TestSuite(C).run() # long time
            sage: C = crystals.LSPaths(['E',6], [1,0,0,0,0,0])
            sage: TestSuite(C).run()

            sage: R = RootSystem(['C',3,1])
            sage: La = R.weight_space().basis()
            sage: LaE = R.weight_space(extended=True).basis()
            sage: B = crystals.LSPaths(La[0])
            sage: BE = crystals.LSPaths(LaE[0])
            sage: B is BE
            False
            sage: B.weight_lattice_realization()
            Weight space over the Rational Field of the Root system of type ['C', 3, 1]
            sage: BE.weight_lattice_realization()
            Extended weight space over the Rational Field of the Root system of type ['C', 3, 1]
        """
        cartan_type = starting_weight.parent().cartan_type()
        self.R = RootSystem(cartan_type)
        self.weight = starting_weight
        if not self.weight.parent().base_ring().has_coerce_map_from(QQ):
             raise ValueError("Please use the weight space, rather than weight lattice for your weights")
        self._cartan_type = cartan_type
        self._name = "The crystal of LS paths of type %s and weight %s"%(cartan_type,starting_weight)
        if cartan_type.is_affine():
            if all(i>=0 for i in starting_weight.coefficients()):
                Parent.__init__( self, category=(RegularCrystals(),
                                                 HighestWeightCrystals(),
                                                 InfiniteEnumeratedSets()) )
            elif starting_weight.parent().is_extended():
                Parent.__init__(self, category=(RegularCrystals(), InfiniteEnumeratedSets()))
            else:
                cl = self._cartan_type.classical().index_set()
                if sum(self.weight[i] for i in cl) == 1:
                    cat = KirillovReshetikhinCrystals()
                else:
                    cat = RegularLoopCrystals().Finite()
                Parent.__init__(self, category=cat)
        else:
            Parent.__init__(self, category=ClassicalCrystals())

        if starting_weight == starting_weight.parent().zero():
            initial_element = self(())
        else:
            initial_element = self((starting_weight,))
        self.module_generators = (initial_element,)

    def _repr_(self):
        """
        EXAMPLES::

            sage: crystals.LSPaths(['B',3],[1,1,0]) # indirect doctest
            The crystal of LS paths of type ['B', 3] and weight Lambda[1] + Lambda[2]
        """
        return self._name

    def weight_lattice_realization(self):
        r"""
        Return weight lattice realization of ``self``.

        EXAMPLES::

            sage: B = crystals.LSPaths(['B',3],[1,1,0])
            sage: B.weight_lattice_realization()
            Weight space over the Rational Field of the Root system of type ['B', 3]
            sage: B = crystals.LSPaths(['B',3,1],[1,1,1,0])
            sage: B.weight_lattice_realization()
            Extended weight space over the Rational Field of the Root system of type ['B', 3, 1]
        """
        return self.weight.parent()

    class Element(ElementWrapper):
        """
        TESTS::

            sage: C = crystals.LSPaths(['E',6],[1,0,0,0,0,0])
            sage: c = C.an_element()
            sage: TestSuite(c).run()
        """

        def endpoint(self):
            r"""
            Computes the endpoint of the path.

            EXAMPLES::

                sage: C = crystals.LSPaths(['A',2],[1,1])
                sage: b = C.module_generators[0]
                sage: b.endpoint()
                Lambda[1] + Lambda[2]
                sage: b.f_string([1,2,2,1])
                (-Lambda[1] - Lambda[2],)
                sage: b.f_string([1,2,2,1]).endpoint()
                -Lambda[1] - Lambda[2]
                sage: b.f_string([1,2])
                (1/2*Lambda[1] - Lambda[2], -1/2*Lambda[1] + Lambda[2])
                sage: b.f_string([1,2]).endpoint()
                0
                sage: b = C([])
                sage: b.endpoint()
                0
            """
            if len(self.value) > 0:
                return sum(self.value)
            return self.parent().weight.parent().zero()
            #return self.parent().R.weight_space(extended = self.parent().extended).zero()

        def compress(self):
            r"""
            Merges consecutive positively parallel steps present in the path.

            EXAMPLES::

                sage: C = crystals.LSPaths(['A',2],[1,1])
                sage: Lambda = C.R.weight_space().fundamental_weights(); Lambda
                Finite family {1: Lambda[1], 2: Lambda[2]}
                sage: c = C(tuple([1/2*Lambda[1]+1/2*Lambda[2], 1/2*Lambda[1]+1/2*Lambda[2]]))
                sage: c.compress()
                (Lambda[1] + Lambda[2],)
            """
            if not self.value:
                return self
            q = []
            curr = self.value[0]
            for v in self.value[1:]:
                if positively_parallel_weights(curr, v):
                    curr = curr + v
                else:
                    q.append(curr)
                    curr = v
            q.append(curr)
            return self.parent()(tuple(q))

        def split_step(self, which_step, r):
            r"""
            Splits indicated step into two parallel steps of relative lengths `r` and `1-r`.

            INPUT:

            - ``which_step`` -- a position in the tuple ``self``
            - ``r`` -- a rational number between 0 and 1

            EXAMPLES::

                sage: C = crystals.LSPaths(['A',2],[1,1])
                sage: b = C.module_generators[0]
                sage: b.split_step(0,1/3)
                (1/3*Lambda[1] + 1/3*Lambda[2], 2/3*Lambda[1] + 2/3*Lambda[2])
            """
            assert 0 <= which_step and which_step <= len(self.value)
            v = self.value[which_step]
            return self.parent()(self.value[:which_step] + (r*v,(1-r)*v) + self.value[which_step+1:])

        def reflect_step(self, which_step, i):
            r"""
            Apply the `i`-th simple reflection to the indicated step in ``self``.

            EXAMPLES::

                sage: C = crystals.LSPaths(['A',2],[1,1])
                sage: b = C.module_generators[0]
                sage: b.reflect_step(0,1)
                (-Lambda[1] + 2*Lambda[2],)
                sage: b.reflect_step(0,2)
                (2*Lambda[1] - Lambda[2],)
            """
            assert i in self.index_set()
            assert 0 <= which_step and which_step <= len(self.value)
            return self.parent()(self.value[:which_step]+tuple([self.value[which_step].simple_reflection(i)])+self.value[which_step+1:])

        def _string_data(self, i):
            r"""
            Computes the `i`-string data of ``self``.

            TESTS::

                sage: C = crystals.LSPaths(['A',2],[1,1])
                sage: b = C.module_generators[0]
                sage: b._string_data(1)
                ()
                sage: b._string_data(2)
                ()
                sage: b.f(1)._string_data(1)
                ((0, -1, -1),)
                sage: b.f(1).f(2)._string_data(2)
                ((0, -1, -1),)
            """
            if not self.value:
                return ()
            # get the i-th simple coroot
            alv = self.value[0].parent().alphacheck()[i]
            # Compute the i-heights of the steps of vs
            steps = [v.scalar(alv) for v in self.value]
            # Get the wet step data
            minima_pos = []
            ps = 0
            psmin = 0
            for ix, step in enumerate(steps):
                ps = ps + step
                if ps < psmin:
                    minima_pos.append((ix,ps,step))
                    psmin = ps
            return tuple(minima_pos)

        def epsilon(self, i):
            r"""
            Returns the distance to the beginning of the `i`-string.

            This method overrides the generic implementation in the category of crystals
            since this computation is more efficient.

            EXAMPLES::

                sage: C = crystals.LSPaths(['A',2],[1,1])
                sage: [c.epsilon(1) for c in C]
                [0, 1, 0, 0, 1, 0, 1, 2]
                sage: [c.epsilon(2) for c in C]
                [0, 0, 1, 2, 1, 1, 0, 0]
            """
            return self.e(i,length_only=True)

        def phi(self, i):
            r"""
            Returns the distance to the end of the `i`-string.

            This method overrides the generic implementation in the category of crystals
            since this computation is more efficient.

            EXAMPLES::

                sage: C = crystals.LSPaths(['A',2],[1,1])
                sage: [c.phi(1) for c in C]
                [1, 0, 0, 1, 0, 2, 1, 0]
                sage: [c.phi(2) for c in C]
                [1, 2, 1, 0, 0, 0, 0, 1]
            """
            return self.f(i,length_only=True)

        def e(self, i, power=1, to_string_end=False, length_only=False):
            r"""
            Returns the `i`-th crystal raising operator on ``self``.

            INPUT:

            - ``i`` -- element of the index set of the underlying root system
            - ``power`` -- positive integer; specifies the power of the raising operator
              to be applied (default: 1)
            - ``to_string_end`` -- boolean; if set to True, returns the dominant end of the
              `i`-string of ``self``. (default: False)
            - ``length_only`` -- boolean; if set to True, returns the distance to the dominant
              end of the `i`-string of ``self``.

            EXAMPLES::

                sage: C = crystals.LSPaths(['A',2],[1,1])
                sage: c = C[2]; c
                (1/2*Lambda[1] - Lambda[2], -1/2*Lambda[1] + Lambda[2])
                sage: c.e(1)
                sage: c.e(2)
                (-Lambda[1] + 2*Lambda[2],)
                sage: c.e(2,to_string_end=True)
                (-Lambda[1] + 2*Lambda[2],)
                sage: c.e(1,to_string_end=True)
                (1/2*Lambda[1] - Lambda[2], -1/2*Lambda[1] + Lambda[2])
                sage: c.e(1,length_only=True)
                0
            """
            assert i in self.index_set()
            data = self._string_data(i)
            # compute the minimum i-height M on the path
            if not data:
                M = 0
            else:
                M = data[-1][1]
            max_raisings = floor(-M)
            if length_only:
                return max_raisings
            # set the power of e_i to apply
            if to_string_end:
                p = max_raisings
            else:
                p = power
            if p > max_raisings:
                return None

            # copy the vector sequence into a working vector sequence ws
            #!!! ws only needs to be the actual vector sequence, not some
            #!!! fancy crystal graph element
            ws = self.parent()(self.value)

            ix = len(data)-1
            while ix >= 0 and data[ix][1] < M + p:
            # get the index of the current step to be processed
                j = data[ix][0]
                # find the i-height where the current step might need to be split
                if ix == 0:
                    prev_ht = M + p
                else:
                    prev_ht = min(data[ix-1][1],M+p)
                # if necessary split the step. Then reflect the wet part.
                if data[ix][1] - data[ix][2] > prev_ht:
                    ws = ws.split_step(j,1-(prev_ht-data[ix][1])/(-data[ix][2]))
                    ws = ws.reflect_step(j+1,i)
                else:
                    ws = ws.reflect_step(j,i)
                ix = ix - 1
            #!!! at this point we should return the fancy crystal graph element
            #!!! corresponding to the humble vector sequence ws
            return self.parent()(ws.compress())

        def dualize(self):
            r"""
            Returns dualized path.

            EXAMPLES::

                sage: C = crystals.LSPaths(['A',2],[1,1])
                sage: for c in C:
                ....:     print("{} {}".format(c, c.dualize()))
                (Lambda[1] + Lambda[2],) (-Lambda[1] - Lambda[2],)
                (-Lambda[1] + 2*Lambda[2],) (Lambda[1] - 2*Lambda[2],)
                (1/2*Lambda[1] - Lambda[2], -1/2*Lambda[1] + Lambda[2]) (1/2*Lambda[1] - Lambda[2], -1/2*Lambda[1] + Lambda[2])
                (Lambda[1] - 2*Lambda[2],) (-Lambda[1] + 2*Lambda[2],)
                (-Lambda[1] - Lambda[2],) (Lambda[1] + Lambda[2],)
                (2*Lambda[1] - Lambda[2],) (-2*Lambda[1] + Lambda[2],)
                (-Lambda[1] + 1/2*Lambda[2], Lambda[1] - 1/2*Lambda[2]) (-Lambda[1] + 1/2*Lambda[2], Lambda[1] - 1/2*Lambda[2])
                (-2*Lambda[1] + Lambda[2],) (2*Lambda[1] - Lambda[2],)
            """
            if not self.value:
                return self
            dual_path = [-v for v in reversed(self.value)]
            return self.parent()(tuple(dual_path))

        def f(self, i, power=1, to_string_end=False, length_only=False):
            r"""
            Returns the `i`-th crystal lowering operator on ``self``.

            INPUT:

            - ``i`` -- element of the index set of the underlying root system
            - ``power`` -- positive integer; specifies the power of the lowering operator
              to be applied (default: 1)
            - ``to_string_end`` -- boolean; if set to True, returns the anti-dominant end of the
              `i`-string of ``self``. (default: False)
            - ``length_only`` -- boolean; if set to True, returns the distance to the anti-dominant
              end of the `i`-string of ``self``.

            EXAMPLES::

                sage: C = crystals.LSPaths(['A',2],[1,1])
                sage: c = C.module_generators[0]
                sage: c.f(1)
                (-Lambda[1] + 2*Lambda[2],)
                sage: c.f(1,power=2)
                sage: c.f(2)
                (2*Lambda[1] - Lambda[2],)
                sage: c.f(2,to_string_end=True)
                (2*Lambda[1] - Lambda[2],)
                sage: c.f(2,length_only=True)
                1

                sage: C = crystals.LSPaths(['A',2,1],[-1,-1,2])
                sage: c = C.module_generators[0]
                sage: c.f(2,power=2)
                (Lambda[0] + Lambda[1] - 2*Lambda[2],)
            """
            dual_path = self.dualize()
            dual_path = dual_path.e(i, power, to_string_end, length_only)
            if length_only:
                return dual_path
            if dual_path is None:
                return None
            return dual_path.dualize()

        def s(self, i):
            r"""
            Computes the reflection of ``self`` along the `i`-string.

            This method is more efficient than the generic implementation since it uses
            powers of `e` and `f` in the Littelmann model directly.

            EXAMPLES::

                sage: C = crystals.LSPaths(['A',2],[1,1])
                sage: c = C.module_generators[0]
                sage: c.s(1)
                (-Lambda[1] + 2*Lambda[2],)
                sage: c.s(2)
                (2*Lambda[1] - Lambda[2],)

                sage: C = crystals.LSPaths(['A',2,1],[-1,0,1])
                sage: c = C.module_generators[0]; c
                (-Lambda[0] + Lambda[2],)
                sage: c.s(2)
                (Lambda[1] - Lambda[2],)
                sage: c.s(1)
                (-Lambda[0] + Lambda[2],)
                sage: c.f(2).s(1)
                (Lambda[0] - Lambda[1],)
            """
            ph = self.phi(i)
            ep = self.epsilon(i)
            diff = ph - ep
            if diff >= 0:
                return self.f(i, power=diff)
            else:
                return self.e(i, power=-diff)

        def weight(self):
            """
            Return the weight of ``self``.

            EXAMPLES::

                sage: B = crystals.LSPaths(['A',1,1],[1,0])
                sage: b = B.highest_weight_vector()
                sage: b.f(0).weight()
                -Lambda[0] + 2*Lambda[1] - delta
            """
            P = self.parent().weight_lattice_realization()
            return sum([p for p in self.value], P.zero())

        def _latex_(self):
            r"""
            Latex method for ``self``.

            EXAMPLES::

                sage: C = crystals.LSPaths(['A',2],[1,1])
                sage: c = C.module_generators[0]
                sage: c._latex_()
                [\Lambda_{1} + \Lambda_{2}]
            """
            return [latex(p) for p in self.value]


#####################################################################
## Projected level-zero


class CrystalOfProjectedLevelZeroLSPaths(CrystalOfLSPaths):
    r"""
    Crystal of projected level zero LS paths.

    INPUT:

    - ``weight`` -- a dominant weight of the weight space of an affine
      Kac-Moody root system

    When ``weight`` is just a single fundamental weight `\Lambda_r`, this
    crystal is isomorphic to a Kirillov-Reshetikhin (KR) crystal, see also
    :meth:`sage.combinat.crystals.kirillov_reshetikhin.KirillovReshetikhinFromLSPaths`.
    For general weights, it is isomorphic to a tensor product of
    single-column KR crystals.

    EXAMPLES::

        sage: R = RootSystem(['C',3,1])
        sage: La = R.weight_space().basis()
        sage: LS = crystals.ProjectedLevelZeroLSPaths(La[1]+La[3])
        sage: LS.cardinality()
        84
        sage: GLS = LS.digraph()

        sage: K1 = crystals.KirillovReshetikhin(['C',3,1],1,1)
        sage: K3 = crystals.KirillovReshetikhin(['C',3,1],3,1)
        sage: T = crystals.TensorProduct(K3,K1)
        sage: T.cardinality()
        84
        sage: GT = T.digraph() # long time
        sage: GLS.is_isomorphic(GT, edge_labels = True) # long time
        True

    TESTS::

        sage: ct = CartanType(['A',4,2]).dual()
        sage: P = RootSystem(ct).weight_space()
        sage: La = P.fundamental_weights()
        sage: C = crystals.ProjectedLevelZeroLSPaths(La[1])
        sage: sorted(C, key=str)
        [(-Lambda[0] + Lambda[1],),
         (-Lambda[1] + 2*Lambda[2],),
         (1/2*Lambda[1] - Lambda[2], -1/2*Lambda[1] + Lambda[2]),
         (Lambda[0] - Lambda[1],),
         (Lambda[1] - 2*Lambda[2],)]
    """

    @staticmethod
    def __classcall_private__(cls, weight):
        """
        Classcall to mend the input.

        Internally, the
        :class:`~sage.combinat.crystals.littelmann_path.CrystalOfProjectedLevelZeroLSPaths`
        uses a level zero weight, which is passed on to
        :class:`~sage.combinat.crystals.littelmann_path.CrystalOfLSPaths`.
        ``weight`` is first coerced to a level zero weight.

        TESTS::

            sage: R = RootSystem(['C',3,1])
            sage: La = R.weight_space().basis()
            sage: C = crystals.ProjectedLevelZeroLSPaths(La[1] + La[2])
            sage: C2 = crystals.ProjectedLevelZeroLSPaths(La[1] + La[2])
            sage: C is C2
            True

            sage: R = RootSystem(['C',3,1])
            sage: La = R.weight_space(extended = True).basis()
            sage: crystals.ProjectedLevelZeroLSPaths(La[1] + La[2])
            Traceback (most recent call last):
            ...
            ValueError: The weight should be in the non-extended weight lattice!
        """
        if weight.parent().is_extended():
            raise ValueError("The weight should be in the non-extended weight lattice!")
        La = weight.parent().basis()
        weight = weight - weight.level() * La[0] / La[0].level()
        return super(CrystalOfLSPaths, cls).__classcall__(cls, weight, starting_weight_parent = weight.parent())

    @cached_method
    def maximal_vector(self):
        """
        Return the maximal vector of ``self``.

        EXAMPLES::

            sage: R = RootSystem(['A',2,1])
            sage: La = R.weight_space().basis()
            sage: LS = crystals.ProjectedLevelZeroLSPaths(2*La[1]+La[2])
            sage: LS.maximal_vector()
            (-3*Lambda[0] + 2*Lambda[1] + Lambda[2],)
        """
        return self.module_generators[0]

    @cached_method
    def classically_highest_weight_vectors(self):
        r"""
        Return the classically highest weight vectors of ``self``.

        EXAMPLES::

            sage: R = RootSystem(['A',2,1])
            sage: La = R.weight_space().basis()
            sage: LS = crystals.ProjectedLevelZeroLSPaths(2*La[1])
            sage: LS.classically_highest_weight_vectors()
            ((-2*Lambda[0] + 2*Lambda[1],),
             (-Lambda[0] + Lambda[1], -Lambda[1] + Lambda[2]))
        """
        I0 = self.cartan_type().classical().index_set()
        return tuple([x for x in self.list() if x.is_highest_weight(I0)])

    def one_dimensional_configuration_sum(self, q=None, group_components=True):
        r"""
        Compute the one-dimensional configuration sum.

        INPUT:

        - ``q`` -- (default: ``None``) a variable or ``None``; if ``None``,
          a variable ``q`` is set in the code
        - ``group_components`` -- (default: ``True``) boolean; if ``True``,
          then the terms are grouped by classical component

        The one-dimensional configuration sum is the sum of the weights
        of all elements in the crystal weighted by the energy function.
        For untwisted types it uses the parabolic quantum Bruhat graph,
        see [LNSSS2013]_. In the dual-of-untwisted case, the parabolic
        quantum Bruhat graph is defined by exchanging the roles of roots
        and coroots (which is still conjectural at this point).

        EXAMPLES::

            sage: R = RootSystem(['A',2,1])
            sage: La = R.weight_space().basis()
            sage: LS = crystals.ProjectedLevelZeroLSPaths(2*La[1])
            sage: LS.one_dimensional_configuration_sum() # long time
            B[-2*Lambda[1] + 2*Lambda[2]] + (q+1)*B[-Lambda[1]]
             + (q+1)*B[Lambda[1] - Lambda[2]] + B[2*Lambda[1]]
             + B[-2*Lambda[2]] + (q+1)*B[Lambda[2]]
            sage: R.<t> = ZZ[]
            sage: LS.one_dimensional_configuration_sum(t, False) # long time
            B[-2*Lambda[1] + 2*Lambda[2]] + (t+1)*B[-Lambda[1]]
             + (t+1)*B[Lambda[1] - Lambda[2]] + B[2*Lambda[1]]
             + B[-2*Lambda[2]] + (t+1)*B[Lambda[2]]

        TESTS::

            sage: R = RootSystem(['B',3,1])
            sage: La = R.weight_space().basis()
            sage: LS = crystals.ProjectedLevelZeroLSPaths(La[1]+La[2])
            sage: LS.one_dimensional_configuration_sum() == LS.one_dimensional_configuration_sum(group_components=False) # long time
            True
            sage: K1 = crystals.KirillovReshetikhin(['B',3,1],1,1)
            sage: K2 = crystals.KirillovReshetikhin(['B',3,1],2,1)
            sage: T = crystals.TensorProduct(K2,K1)
            sage: T.one_dimensional_configuration_sum() == LS.one_dimensional_configuration_sum() # long time
            True

            sage: R = RootSystem(['D',4,2])
            sage: La = R.weight_space().basis()
            sage: LS = crystals.ProjectedLevelZeroLSPaths(La[1]+La[2])
            sage: K1 = crystals.KirillovReshetikhin(['D',4,2],1,1)
            sage: K2 = crystals.KirillovReshetikhin(['D',4,2],2,1)
            sage: T = crystals.TensorProduct(K2,K1)
            sage: T.one_dimensional_configuration_sum() == LS.one_dimensional_configuration_sum() # long time
            True

            sage: R = RootSystem(['A',5,2])
            sage: La = R.weight_space().basis()
            sage: LS = crystals.ProjectedLevelZeroLSPaths(3*La[1])
            sage: K1 = crystals.KirillovReshetikhin(['A',5,2],1,1)
            sage: T = crystals.TensorProduct(K1,K1,K1)
            sage: T.one_dimensional_configuration_sum() == LS.one_dimensional_configuration_sum() # long time
            True
        """
        if q is None:
            from sage.rings.rational_field import QQ
            q = QQ['q'].gens()[0]
        #P0 = self.weight_lattice_realization().classical()
        P0 = RootSystem(self.cartan_type().classical()).weight_lattice()
        B = P0.algebra(q.parent())
        def weight(x):
            w = x.weight()
            return P0.sum(int(c)*P0.basis()[i] for i,c in w if i in P0.index_set())
        if group_components:
            G = self.digraph(index_set = self.cartan_type().classical().index_set())
            C = G.connected_components()
            return sum(q**(c[0].energy_function())*B.sum(B(weight(b)) for b in c) for c in C)
        return B.sum(q**(b.energy_function())*B(weight(b)) for b in self)

    def is_perfect(self, level=1):
        r"""
        Check whether the crystal ``self`` is perfect (of level ``level``).

        INPUT:

        - ``level`` -- (default: 1) positive integer

        A crystal `\mathcal{B}` is perfect of level `\ell` if:

        #. `\mathcal{B}` is isomorphic to the crystal graph of a
           finite-dimensional `U_q^{'}(\mathfrak{g})`-module.
        #. `\mathcal{B}\otimes \mathcal{B}` is connected.
        #. There exists a `\lambda\in X`, such that
           `\mathrm{wt}(\mathcal{B}) \subset \lambda + \sum_{i\in I} \ZZ_{\le 0} \alpha_i`
           and there is a unique element in
           `\mathcal{B}` of classical weight `\lambda`.
        #. For all `b \in \mathcal{B}`,
           `\mathrm{level}(\varepsilon (b)) \geq \ell`.
        #. For all `\Lambda` dominant weights of level `\ell`, there exist
           unique elements `b_{\Lambda}, b^{\Lambda} \in \mathcal{B}`, such
           that `\varepsilon (b_{\Lambda}) = \Lambda = \varphi(b^{\Lambda})`.

        Points (1)-(3) are known to hold. This method checks points (4) and (5).

        EXAMPLES::

            sage: C = CartanType(['C',2,1])
            sage: R = RootSystem(C)
            sage: La = R.weight_space().basis()
            sage: LS = crystals.ProjectedLevelZeroLSPaths(La[1])
            sage: LS.is_perfect()
            False
            sage: LS = crystals.ProjectedLevelZeroLSPaths(La[2])
            sage: LS.is_perfect()
            True

            sage: C = CartanType(['E',6,1])
            sage: R = RootSystem(C)
            sage: La = R.weight_space().basis()
            sage: LS = crystals.ProjectedLevelZeroLSPaths(La[1])
            sage: LS.is_perfect()
            True
            sage: LS.is_perfect(2)
            False

            sage: C = CartanType(['D',4,1])
            sage: R = RootSystem(C)
            sage: La = R.weight_space().basis()
            sage: all(crystals.ProjectedLevelZeroLSPaths(La[i]).is_perfect() for i in [1,2,3,4])
            True

            sage: C = CartanType(['A',6,2])
            sage: R = RootSystem(C)
            sage: La = R.weight_space().basis()
            sage: LS = crystals.ProjectedLevelZeroLSPaths(La[1]+La[2])
            sage: LS.is_perfect()
            True
            sage: LS.is_perfect(2)
            False
        """
        MPhi = []
        for b in self:
            p = b.Phi().level()
            assert p == b.Epsilon().level()
            if p < level:
                return False
            if p == level:
                MPhi += [b]
        weights = []
        I = self.index_set()
        rank = len(I)
        La = self.weight_lattice_realization().basis()
        from sage.combinat.integer_vector import IntegerVectors
        for n in range(1,level+1):
            for c in IntegerVectors(n, rank):
                w = sum(c[i]*La[i] for i in I)
                if w.level() == level:
                    weights.append(w)
        return sorted([b.Phi() for b in MPhi]) == sorted(weights)

    class Element(CrystalOfLSPaths.Element):
        """
        Element of a crystal of projected level zero LS paths.
        """

        @cached_in_parent_method
        def scalar_factors(self):
            r"""
            Obtain the scalar factors for ``self``.

            Each LS path (or ``self``) can be written as a piecewise linear map

            .. MATH::

                \pi(t) = \sum_{u'=1}^{u-1} (\sigma_{u'} - \sigma_{u'-1}) \nu_{u'}
                + (t-\sigma_{u-1}) \nu_{u}

            for `0 < \sigma_1 < \sigma_2 < \cdots < \sigma_s=1` and
            `\sigma_{u-1} \le t \le \sigma_{u}` and `1 \le u \le s`.
            This method returns the tuple of `(\sigma_1,\ldots,\sigma_s)`.

            EXAMPLES::

                sage: R = RootSystem(['C',3,1])
                sage: La = R.weight_space().basis()
                sage: LS = crystals.ProjectedLevelZeroLSPaths(La[1]+La[3])
                sage: b = LS.module_generators[0]
                sage: b.scalar_factors()
                [1]
                sage: c = b.f(1).f(3).f(2)
                sage: c.scalar_factors()
                [1/3, 1]
            """
            weight = self.parent().weight
            l = []
            s = 0
            for c in self.value:
                supp = c.support()
                if supp:
                    i = supp[0]
                    for w in weight._orbit_iter():
                        # Check whether the vectors c and w are positive scalar multiples of each other
                        # If i is not in the support of w, then the first
                        #   product is 0
                        if c[i] * w[i] > 0 and c[i] * w == w[i] * c:
                            s += c[i] / w[i]
                            l += [s]
                            break
            return l

        @cached_in_parent_method
        def weyl_group_representation(self):
            r"""
            Transform the weights in the LS path ``self`` to elements
            in the Weyl group.

            Each LS path can be written as the piecewise linear map:

            .. MATH::

                \pi(t) = \sum_{u'=1}^{u-1} (\sigma_{u'} - \sigma_{u'-1}) \nu_{u'}
                + (t-\sigma_{u-1}) \nu_{u}

            for `0 < \sigma_1 < \sigma_2 < \cdots < \sigma_s = 1` and
            `\sigma_{u-1} \le t \le \sigma_{u}` and `1 \le u \le s`.
            Each weight `\nu_u` is also associated to a Weyl group element.
            This method returns the list of Weyl group elements associated
            to the `\nu_u` for `1\le u\le s`.

            EXAMPLES::

                sage: R = RootSystem(['C',3,1])
                sage: La = R.weight_space().basis()
                sage: LS = crystals.ProjectedLevelZeroLSPaths(La[1]+La[3])
                sage: b = LS.module_generators[0]
                sage: c = b.f(1).f(3).f(2)
                sage: c.weyl_group_representation()
                [s2*s1*s3, s1*s3]
            """
            cartan = self.parent().weight.parent().cartan_type().classical()
            I = cartan.index_set()
            W = WeylGroup(cartan, prefix='s', implementation="permutation")
            return [W.from_reduced_word(x.to_dominant_chamber(index_set=I, reduced_word=True)[1]) for x in self.value]

        @cached_in_parent_method
        def energy_function(self):
            r"""
            Return the energy function of ``self``.

            The energy function `D(\pi)` of the level zero LS path
            `\pi \in \mathbb{B}_\mathrm{cl}(\lambda)` requires a series
            of definitions; for simplicity the root system is assumed to
            be untwisted affine.

            The LS path `\pi` is a piecewise linear map from the unit
            interval `[0,1]` to the weight lattice. It is specified by
            "times" `0 = \sigma_0 < \sigma_1 < \dotsm < \sigma_s = 1` and
            "direction vectors" `x_u \lambda` where `x_u \in W / W_J` for
            `1 \le u \le s`, and `W_J` is the stabilizer of `\lambda` in
            the finite Weyl group `W`. Precisely,

            .. MATH::

                \pi(t) = \sum_{u'=1}^{u-1} (\sigma_{u'}-\sigma_{u'-1})
                x_{u'} \lambda + (t-\sigma_{u-1}) x_{u} \lambda

            for `1 \le u \le s` and `\sigma_{u-1} \le t \le \sigma_{u}`.

            For any `x,y \in W / W_J`, let

            .. MATH::

                d: x = w_{0} \stackrel{\beta_{1}}{\leftarrow}
                w_{1} \stackrel{\beta_{2}}{\leftarrow} \cdots
                \stackrel{\beta_{n}}{\leftarrow} w_{n}=y

            be a shortest directed path in the parabolic quantum
            Bruhat graph. Define

            .. MATH::

                \mathrm{wt}(d) := \sum_{\substack{1 \le k \le n
                \\ \ell(w_{k-1}) < \ell(w_k)}}
                \beta_{k}^{\vee}.

            It can be shown that `\mathrm{wt}(d)` depends only on `x,y`;
            call its value `\mathrm{wt}(x,y)`. The energy function `D(\pi)`
            is defined by

            .. MATH::

                D(\pi) = -\sum_{u=1}^{s-1} (1-\sigma_{u}) \langle \lambda,
                \mathrm{wt}(x_u,x_{u+1}) \rangle.

            For more information, see [LNSSS2013]_.

            .. NOTE::

                In the dual-of-untwisted case the parabolic quantum
                Bruhat graph that is used is obtained by exchanging the
                roles of roots and coroots. Moreover, in the computation
                of the pairing the short roots must be doubled (or tripled
                for type `G`). This factor is determined by the translation
                factor of the corresponding root. Type `BC` is viewed as
                untwisted type, whereas the dual of `BC` is viewed as twisted.
                Except for the untwisted cases, these formulas are
                currently still conjectural.

            EXAMPLES::

                sage: R = RootSystem(['C',3,1])
                sage: La = R.weight_space().basis()
                sage: LS = crystals.ProjectedLevelZeroLSPaths(La[1]+La[3])
                sage: b = LS.module_generators[0]
                sage: c = b.f(1).f(3).f(2)
                sage: c.energy_function()
                0
                sage: c=b.e(0)
                sage: c.energy_function()
                1

                sage: R = RootSystem(['A',2,1])
                sage: La = R.weight_space().basis()
                sage: LS = crystals.ProjectedLevelZeroLSPaths(2*La[1])
                sage: b = LS.module_generators[0]
                sage: c = b.e(0)
                sage: c.energy_function()
                1
                sage: for c in sorted(LS, key=str):
                ....:     print("{} {}".format(c,c.energy_function()))
                (-2*Lambda[0] + 2*Lambda[1],)                    0
                (-2*Lambda[1] + 2*Lambda[2],)                    0
                (-Lambda[0] + Lambda[1], -Lambda[1] + Lambda[2]) 1
                (-Lambda[0] + Lambda[1], Lambda[0] - Lambda[2])  1
                (-Lambda[1] + Lambda[2], -Lambda[0] + Lambda[1]) 0
                (-Lambda[1] + Lambda[2], Lambda[0] - Lambda[2])  1
                (2*Lambda[0] - 2*Lambda[2],)                     0
                (Lambda[0] - Lambda[2], -Lambda[0] + Lambda[1])  0
                (Lambda[0] - Lambda[2], -Lambda[1] + Lambda[2])  0

            The next test checks that the energy function is constant
            on classically connected components::

                sage: R = RootSystem(['A',2,1])
                sage: La = R.weight_space().basis()
                sage: LS = crystals.ProjectedLevelZeroLSPaths(2*La[1]+La[2])
                sage: G = LS.digraph(index_set=[1,2])
                sage: C = G.connected_components()
                sage: [all(c[0].energy_function()==a.energy_function() for a in c) for c in C]
                [True, True, True, True]

                sage: R = RootSystem(['D',4,2])
                sage: La = R.weight_space().basis()
                sage: LS = crystals.ProjectedLevelZeroLSPaths(La[2])
                sage: J = R.cartan_type().classical().index_set()
                sage: hw = [x for x in LS if x.is_highest_weight(J)]
                sage: [(x.weight(), x.energy_function()) for x in hw]
                [(-2*Lambda[0] + Lambda[2], 0), (-2*Lambda[0] + Lambda[1], 1), (0, 2)]
                sage: G = LS.digraph(index_set=J)
                sage: C = G.connected_components()
                sage: [all(c[0].energy_function()==a.energy_function() for a in c) for c in C]
                [True, True, True]

                sage: R = RootSystem(CartanType(['G',2,1]).dual())
                sage: La = R.weight_space().basis()
                sage: LS = crystals.ProjectedLevelZeroLSPaths(La[1]+La[2])
                sage: G = LS.digraph(index_set=[1,2])
                sage: C = G.connected_components()
                sage: [all(c[0].energy_function()==a.energy_function() for a in c) for c in C] # long time
                [True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True]

                sage: ct = CartanType(['BC',2,2]).dual()
                sage: R = RootSystem(ct)
                sage: La = R.weight_space().basis()
                sage: LS = crystals.ProjectedLevelZeroLSPaths(2*La[1]+La[2])
                sage: G = LS.digraph(index_set=R.cartan_type().classical().index_set())
                sage: C = G.connected_components()
                sage: [all(c[0].energy_function()==a.energy_function() for a in c) for c in C] # long time
                [True, True, True, True, True, True, True, True, True, True, True]

                sage: R = RootSystem(['BC',2,2])
                sage: La = R.weight_space().basis()
                sage: LS = crystals.ProjectedLevelZeroLSPaths(2*La[1]+La[2])
                sage: G = LS.digraph(index_set=R.cartan_type().classical().index_set())
                sage: C = G.connected_components()
                sage: [all(c[0].energy_function()==a.energy_function() for a in c) for c in C] # long time
                [True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True]
            """
            weight = self.parent().weight
            P = weight.parent()
            c_weight = P.classical()(weight)
            ct = P.cartan_type()
            cartan = ct.classical()
            Qv = RootSystem(cartan).coroot_lattice()
            W = WeylGroup(cartan, prefix='s', implementation="permutation")
            J = tuple(weight.weyl_stabilizer())
            L = self.weyl_group_representation()
            if ct.is_untwisted_affine() or ct.type() == 'BC':
                untwisted = True
                G = W.quantum_bruhat_graph(J)
            else:
                untwisted = False
                cartan_dual = cartan.dual()
                Wd = WeylGroup(cartan_dual, prefix='s', implementation="permutation")
                G = Wd.quantum_bruhat_graph(J)
                Qd = RootSystem(cartan_dual).root_lattice()
                dualize = lambda x: Qv.from_vector(x.to_vector())
                L = [Wd.from_reduced_word(x.reduced_word()) for x in L]
                def stretch_short_root(a):
                    # stretches roots by translation factor
                    if ct.dual().type() == 'BC':
                        return ct.c()[a.to_simple_root()]*a
                    return ct.dual().c()[a.to_simple_root()]*a
                    #if a.is_short_root():
                    #    if cartan_dual.type() == 'G':
                    #        return 3*a
                    #    else:
                    #        return 2*a
                    #return a
            paths = [G.shortest_path(L[i+1],L[i]) for i in range(len(L)-1)]
            paths_labels = [[G.edge_label(p[i],p[i+1]) for i in range(len(p)-1) if p[i].length()+1 != p[i+1].length()] for p in paths]
            scalars = self.scalar_factors()
            if untwisted:
                s = sum((1-scalars[i])*c_weight.scalar( Qv.sum(root.associated_coroot()
                       for root in paths_labels[i]) ) for i in range(len(paths_labels)))
                if ct.type() == 'BC':
                    return 2*s
                else:
                    return s
            else:
                s = sum((1-scalars[i])*c_weight.scalar( dualize (Qd.sum(stretch_short_root(root) for root in paths_labels[i])) ) for i in range(len(paths_labels)))
                if ct.dual().type() == 'BC':
                    return s/2
                else:
                    return s


#####################################################################
## B(\infty)


class InfinityCrystalOfLSPaths(UniqueRepresentation, Parent):
    r"""
    LS path model for `\mathcal{B}(\infty)`.

    Elements of `\mathcal{B}(\infty)` are equivalence classes of paths `[\pi]`
    in `\mathcal{B}(k\rho)` for `k\gg 0`, where `\rho` is the Weyl vector.  A
    canonical representative for an element of `\mathcal{B}(\infty)` is chosen
    by taking `k` to be minimal such that the endpoint of `\pi` is strictly
    dominant but its representative in `\mathcal{B}((k-1)\rho)` is on the wall
    of the dominant chamber.

    REFERENCES:

    - [LZ2011]_
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: B1 = crystals.infinity.LSPaths(['A',4])
            sage: B2 = crystals.infinity.LSPaths('A4')
            sage: B3 = crystals.infinity.LSPaths(CartanType(['A',4]))
            sage: B1 is B2 and B2 is B3
            True
        """
        cartan_type = CartanType(cartan_type)
        return super(InfinityCrystalOfLSPaths, cls).__classcall__(cls, cartan_type)

    def __init__(self, cartan_type):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.LSPaths(['D',4,3])
            sage: TestSuite(B).run(max_runs=500)
            sage: B = crystals.infinity.LSPaths(['B',3])
            sage: TestSuite(B).run() # long time
        """
        Parent.__init__(self, category=(HighestWeightCrystals(),
                                        InfiniteEnumeratedSets()))
        self._cartan_type = cartan_type
        self.module_generators = (self.module_generator(),)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: crystals.infinity.LSPaths(['A',4])
            The infinity crystal of LS paths of type ['A', 4]
        """
        return "The infinity crystal of LS paths of type %s" % self._cartan_type

    @cached_method
    def module_generator(self):
        r"""
        Return the module generator (or highest weight element) of ``self``.

        The module generator is the unique path
        `\pi_\infty\colon t \mapsto t\rho`, for `t \in [0,\infty)`.

        EXAMPLES::

            sage: B = crystals.infinity.LSPaths(['A',6,2])
            sage: mg = B.module_generator(); mg
            (Lambda[0] + Lambda[1] + Lambda[2] + Lambda[3],)
            sage: mg.weight()
            0
        """
        rho = self.weight_lattice_realization().rho()
        return self((rho,))

    def weight_lattice_realization(self):
        """
        Return the weight lattice realization of ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.LSPaths(['C',4])
            sage: B.weight_lattice_realization()
            Weight space over the Rational Field of the Root system of type ['C', 4]
        """
        if self._cartan_type.is_affine():
            return self._cartan_type.root_system().weight_space(extended=True)
        return self._cartan_type.root_system().weight_space()

    class Element(CrystalOfLSPaths.Element):

        def e(self, i, power=1, length_only=False):
            r"""
            Return the `i`-th crystal raising operator on ``self``.

            INPUT:

            - ``i`` -- element of the index set
            - ``power`` -- (default: 1) positive integer; specifies the
              power of the lowering operator to be applied
            - ``length_only`` -- (default: ``False``) boolean; if ``True``,
              then return the distance to the anti-dominant end of the
              `i`-string of ``self``

            EXAMPLES::

                sage: B = crystals.infinity.LSPaths(['B',3,1])
                sage: mg = B.module_generator()
                sage: mg.e(0)
                sage: mg.e(1)
                sage: mg.e(2)
                sage: x = mg.f_string([1,0,2,1,0,2,1,1,0])
                sage: all(x.f(i).e(i) == x for i in B.index_set())
                True
                sage: all(x.e(i).f(i) == x for i in B.index_set() if x.epsilon(i) > 0)
                True

            TESTS:

            Check that this works in affine types::

                sage: B = crystals.infinity.LSPaths(['A',3,1])
                sage: mg = B.highest_weight_vector()
                sage: x = mg.f_string([0,1,2,3])
                sage: x.e_string([3,2,1,0]) == mg
                True

            We check that :meth:`epsilon` works::

                sage: B = crystals.infinity.LSPaths(['D',4])
                sage: mg = B.highest_weight_vector()
                sage: x = mg.f_string([1,3,4,2,4,3,2,1,4])
                sage: [x.epsilon(i) for i in B.index_set()]
                [1, 1, 0, 1]

            Check that :trac:`21671` is fixed::

                sage: B = crystals.infinity.LSPaths(['G',2])
                sage: len(B.subcrystal(max_depth=7))
                116
            """
            ret = super(InfinityCrystalOfLSPaths.Element, self).e(i, power=power,
                                                                  length_only=length_only)
            if ret is None:
                return None
            if length_only:
                return ret
            WLR = self.parent().weight_lattice_realization()
            value = list(ret.value)
            endpoint = sum(p for p in value)
            rho = WLR.rho()
            h = WLR.simple_coroots()

            if not positively_parallel_weights(value[-1], rho):
                value.append(rho)
                endpoint += rho

            while any(endpoint.scalar(alc) < 1 for alc in h):
                value[-1] += rho
                endpoint += rho
            while all(endpoint.scalar(alc) > 1 for alc in h) and value[-1] != WLR.zero():
                value[-1] -= rho
                endpoint -= rho
            while value[-1] == WLR.zero():
                value.pop()
            ret.value = tuple(value)
            return ret

        def f(self, i, power=1, length_only=False):
            r"""
            Return the `i`-th crystal lowering operator on ``self``.

            INPUT:

            - ``i`` -- element of the index set
            - ``power`` -- (default: 1) positive integer; specifies the
              power of the lowering operator to be applied
            - ``length_only`` -- (default: ``False``) boolean; if ``True``,
              then return the distance to the anti-dominant end of the
              `i`-string of ``self``

            EXAMPLES::

                sage: B = crystals.infinity.LSPaths(['D',3,2])
                sage: mg = B.highest_weight_vector()
                sage: mg.f(1)
                (3*Lambda[0] - Lambda[1] + 3*Lambda[2],
                 2*Lambda[0] + 2*Lambda[1] + 2*Lambda[2])
                sage: mg.f(2)
                (Lambda[0] + 2*Lambda[1] - Lambda[2],
                 2*Lambda[0] + 2*Lambda[1] + 2*Lambda[2])
                sage: mg.f(0)
                (-Lambda[0] + 2*Lambda[1] + Lambda[2] - delta,
                 2*Lambda[0] + 2*Lambda[1] + 2*Lambda[2])
            """
            dual_path = self.dualize()
            dual_path = super(InfinityCrystalOfLSPaths.Element, dual_path).e(i, power, length_only=length_only)
            if length_only:
                return dual_path
            if dual_path is None:
                return None
            ret = dual_path.dualize()
            WLR = self.parent().weight_lattice_realization()
            value = list(ret.value)
            endpoint = sum(p for p in value)
            rho = WLR.rho()
            h = WLR.simple_coroots()

            if not positively_parallel_weights(value[-1], rho):
                value.append(rho)
                endpoint += rho

            while any(endpoint.scalar(alc) < 1 for alc in h):
                value[-1] += rho
                endpoint += rho
            while all(endpoint.scalar(alc) > 1 for alc in h) and value[-1] != WLR.zero():
                value[-1] -= rho
                endpoint -= rho
            while value[-1] == WLR.zero():
                value.pop()
            ret.value = tuple(value)
            return ret

        @cached_method
        def weight(self):
            """
            Return the weight of ``self``.

            .. TODO::

                This is a generic algorithm. We should find a better
                description and implement it.

            EXAMPLES::

                sage: B = crystals.infinity.LSPaths(['E',6])
                sage: mg = B.highest_weight_vector()
                sage: f_seq = [1,4,2,6,4,2,3,1,5,5]
                sage: x = mg.f_string(f_seq)
                sage: x.weight()
                -3*Lambda[1] - 2*Lambda[2] + 2*Lambda[3] + Lambda[4] - Lambda[5]

                sage: al = B.cartan_type().root_system().weight_space().simple_roots()
                sage: x.weight() == -sum(al[i] for i in f_seq)
                True
            """
            WLR = self.parent().weight_lattice_realization()
            alpha = WLR.simple_roots()
            return -WLR.sum(alpha[i] for i in self.to_highest_weight()[1])

        def phi(self,i):
            r"""
            Return `\varphi_i` of ``self``.

            Let `\pi \in \mathcal{B}(\infty)`. Define

            .. MATH::

                \varphi_i(\pi) := \varepsilon_i(\pi) + \langle h_i,
                \mathrm{wt}(\pi) \rangle,

            where `h_i` is the `i`-th simple coroot and `\mathrm{wt}(\pi)`
            is the :meth:`weight` of `\pi`.

            INPUT:

            - ``i`` -- element of the index set

            EXAMPLES::

                sage: B = crystals.infinity.LSPaths(['D',4])
                sage: mg = B.highest_weight_vector()
                sage: x = mg.f_string([1,3,4,2,4,3,2,1,4])
                sage: [x.phi(i) for i in B.index_set()]
                [-1, 4, -2, -3]
            """
            WLR = self.parent().weight_lattice_realization()
            h = WLR.simple_coroots()
            return self.epsilon(i) + WLR(self.weight()).scalar(h[i])


#####################################################################
## Helper functions


def positively_parallel_weights(v, w):
    """
    Check whether the vectors ``v`` and ``w`` are positive scalar
    multiples of each other.

    EXAMPLES::

        sage: from sage.combinat.crystals.littelmann_path import positively_parallel_weights
        sage: La = RootSystem(['A',5,2]).weight_space(extended=True).fundamental_weights()
        sage: rho = sum(La)
        sage: positively_parallel_weights(rho, 4*rho)
        True
        sage: positively_parallel_weights(4*rho, rho)
        True
        sage: positively_parallel_weights(rho, -rho)
        False
        sage: positively_parallel_weights(rho, La[1] + La[2])
        False
    """
    supp = v.support()
    if len(supp) > 0:
        i = supp[0]
        if v[i]*w[i] > 0 and v[i]*w == w[i]*v:
            return True
    return False

