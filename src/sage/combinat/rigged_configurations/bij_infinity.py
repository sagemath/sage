# sage.doctest: needs sage.combinat sage.modules
r"""
Bijection between rigged configurations for `B(\infty)` and marginally large tableaux

AUTHORS:

- Travis Scrimshaw (2015-07-01): Initial version

REFERENCES:

.. [RC-MLT] Ben Salisbury and Travis Scrimshaw. *Connecting marginally
   large tableaux and rigged configurations via crystals*.
   Preprint. :arxiv:`1505.07040`.
"""

# ****************************************************************************
#       Copyright (C) 2015 Travis Scrimshaw <tscrim@ucdavis.edu>
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

from sage.combinat.rigged_configurations.rigged_configurations import RiggedConfigurations

from sage.combinat.rigged_configurations.bij_type_B import (KRTToRCBijectionTypeB,
                                                            RCToKRTBijectionTypeB)
from sage.combinat.rigged_configurations.bij_type_D import (KRTToRCBijectionTypeD,
                                                            RCToKRTBijectionTypeD)
from sage.combinat.rigged_configurations.bij_type_A import (KRTToRCBijectionTypeA,
                                                            RCToKRTBijectionTypeA)
from sage.combinat.rigged_configurations.bij_type_C import (KRTToRCBijectionTypeC,
                                                            RCToKRTBijectionTypeC)
from sage.combinat.rigged_configurations.tensor_product_kr_tableaux import TensorProductOfKirillovReshetikhinTableaux
from sage.combinat.crystals.letters import CrystalOfLetters
from sage.categories.morphism import Morphism
from sage.categories.homset import Hom
from sage.misc.flatten import flatten


class FromTableauIsomorphism(Morphism):
    r"""
    Crystal isomorphism of `B(\infty)` in the tableau model to the
    rigged configuration model.
    """

    def _repr_type(self):
        r"""
        Return the type of morphism of ``self``.

        EXAMPLES::

            sage: RC = crystals.infinity.RiggedConfigurations(['A',3])
            sage: T = crystals.infinity.Tableaux(['A',3])
            sage: phi = RC.coerce_map_from(T)
            sage: phi._repr_type()
            'Crystal Isomorphism'
        """
        return "Crystal Isomorphism"

    def __invert__(self):
        r"""
        Return the inverse of ``self``.

        EXAMPLES::

            sage: RC = crystals.infinity.RiggedConfigurations(['A',3])
            sage: T = crystals.infinity.Tableaux(['A',3])
            sage: phi = RC.coerce_map_from(T)
            sage: ~phi
            Crystal Isomorphism morphism:
              From: The infinity crystal of rigged configurations of type ['A', 3]
              To:   The infinity crystal of tableaux of type ['A', 3]
        """
        return FromRCIsomorphism(Hom(self.codomain(), self.domain()))

    def _call_(self, x):
        r"""
        Return the image of ``x`` in the rigged configuration model
        of `B(\infty)`.

        EXAMPLES::

            sage: RC = crystals.infinity.RiggedConfigurations(['A',3])
            sage: T = crystals.infinity.Tableaux(['A',3])
            sage: phi = RC.coerce_map_from(T)
            sage: x = T.an_element().f_string([2,2,1,1,3,2,1,2,1,3])
            sage: y = phi(x); ascii_art(y)
            -4[ ][ ][ ][ ]-2  -3[ ][ ][ ]-1  -1[ ][ ]-1
                              -2[ ]-1
            sage: (~phi)(y) == x
            True
        """
        conj = x.to_tableau().conjugate()
        ct = self.domain().cartan_type()
        act = ct.affine()
        TP = TensorProductOfKirillovReshetikhinTableaux(act, [[r,1] for r in conj.shape()])
        elt = TP(pathlist=[reversed(row) for row in conj])

        if ct.type() == 'A':
            bij = KRTToRCBijectionTypeA(elt)
        elif ct.type() == 'B':
            bij = MLTToRCBijectionTypeB(elt)
        elif ct.type() == 'C':
            bij = KRTToRCBijectionTypeC(elt)
        elif ct.type() == 'D':
            bij = MLTToRCBijectionTypeD(elt)
        else:
            raise NotImplementedError("bijection of type {} not yet implemented".format(ct))
        return self.codomain()(bij.run())


class FromRCIsomorphism(Morphism):
    r"""
    Crystal isomorphism of `B(\infty)` in the rigged configuration model
    to the tableau model.
    """

    def _repr_type(self):
        r"""
        Return the type of morphism of ``self``.

        EXAMPLES::

            sage: T = crystals.infinity.Tableaux(['A',3])
            sage: RC = crystals.infinity.RiggedConfigurations(['A',3])
            sage: phi = T.coerce_map_from(RC)
            sage: phi._repr_type()
            'Crystal Isomorphism'
        """
        return "Crystal Isomorphism"

    def __invert__(self):
        r"""
        Return the inverse of ``self``.

        EXAMPLES::

            sage: T = crystals.infinity.Tableaux(['A',3])
            sage: RC = crystals.infinity.RiggedConfigurations(['A',3])
            sage: phi = T.coerce_map_from(RC)
            sage: ~phi
            Crystal Isomorphism morphism:
              From: The infinity crystal of tableaux of type ['A', 3]
              To:   The infinity crystal of rigged configurations of type ['A', 3]
        """
        return FromTableauIsomorphism(Hom(self.codomain(), self.domain()))

    def _call_(self, x):
        r"""
        Return the image of ``x`` in the tableau model of `B(\infty)`.

        EXAMPLES::

            sage: T = crystals.infinity.Tableaux(['A',3])
            sage: RC = crystals.infinity.RiggedConfigurations(['A',3])
            sage: phi = T.coerce_map_from(RC)
            sage: x = RC.an_element().f_string([2,2,1,1,3,2,1,2,1,3])
            sage: y = phi(x); y.pp()
              1  1  1  1  1  2  2  3  4
              2  2  3  4
              3
            sage: (~phi)(y) == x
            True
        """
        lam = [sum(nu)+1 for nu in x]
        ct = self.domain().cartan_type()
        I = ct.index_set()
        if ct.type() == 'D':
            lam[-2] = max(lam[-2], lam[-1])
            lam.pop()
            l = sum([[[r+1, 1]]*v for r, v in enumerate(lam[:-1])], [])
            n = len(I)
            l = l + sum([[[n,1], [n-1,1]] for k in range(lam[-1])], [])
        else:
            if ct.type() == 'B':
                lam[-1] *= 2
            l = sum([[[r, 1]]*lam[i] for i, r in enumerate(I)], [])

        RC = RiggedConfigurations(ct.affine(), reversed(l))
        elt = RC(x)
        if ct.type() == 'A':
            bij = RCToKRTBijectionTypeA(elt)
        elif ct.type() == 'B':
            bij = RCToMLTBijectionTypeB(elt)
        elif ct.type() == 'C':
            bij = RCToKRTBijectionTypeC(elt)
        elif ct.type() == 'D':
            bij = RCToMLTBijectionTypeD(elt)
        else:
            raise NotImplementedError("bijection of type {} not yet implemented".format(ct))
        y = bij.run()

        # Now make the result marginally large
        y = [list(c) for c in y]
        cur = []
        L = CrystalOfLetters(ct)
        for i in I:
            cur.insert(0, L(i))
            c = y.count(cur)
            while c > 1:
                y.remove(cur)
                c -= 1
        return self.codomain()(*flatten(y))


class MLTToRCBijectionTypeB(KRTToRCBijectionTypeB):
    def run(self):
        r"""
        Run the bijection from a marginally large tableaux to a rigged
        configuration.

        EXAMPLES::

            sage: vct = CartanType(['B',4]).as_folding()
            sage: RC = crystals.infinity.RiggedConfigurations(vct)
            sage: T = crystals.infinity.Tableaux(['B',4])
            sage: Psi = T.crystal_morphism({T.module_generators[0]: RC.module_generators[0]})
            sage: TS = [x.value for x in T.subcrystal(max_depth=4)]
            sage: all(Psi(b) == RC(b) for b in TS) # long time # indirect doctest
            True
        """
        for cur_crystal in reversed(self.tp_krt):
            cur_column = list(cur_crystal)
            self.cur_path.insert(0, []) # Prepend an empty list
            self.cur_dims.insert(0, [0, 1])

            for letter in reversed(cur_column):
                self.cur_dims[0][0] += 1

                val = letter.value # Convert from a CrystalOfLetter to an Integer

                # Build the next state
                self.cur_path[0].insert(0, [letter]) # Prepend the value
                if self.cur_dims[0][0] == self.n:
                    # Spinor case, we go from \Lambda_{n-1} -> 2\Lambda_n
                    self.cur_dims.insert(1, [self.n,1])
                    self.cur_path.insert(1, self.cur_path[0])

                self.next_state(val)

        self.ret_rig_con.set_immutable() # Return it to immutable
        return self.ret_rig_con


class RCToMLTBijectionTypeB(RCToKRTBijectionTypeB):
    def run(self):
        r"""
        Run the bijection from rigged configurations to a marginally large
        tableau.

        EXAMPLES::

            sage: vct = CartanType(['B',4]).as_folding()
            sage: RC = crystals.infinity.RiggedConfigurations(vct)
            sage: T = crystals.infinity.Tableaux(['B',4])
            sage: Psi = RC.crystal_morphism({RC.module_generators[0]: T.module_generators[0]})
            sage: RCS = [x.value for x in RC.subcrystal(max_depth=4)]
            sage: all(Psi(nu) == T(nu) for nu in RCS) # long time # indirect doctest
            True
        """
        letters = CrystalOfLetters(self.rigged_con.parent()._cartan_type.classical())
        ret_crystal_path = []

        while self.cur_dims:
            dim = self.cur_dims[0]
            ret_crystal_path.append([])

            # Assumption: all factors are single columns
            if dim[0] == self.n:
                # Spinor case, since we've done 2\Lambda_n -> \Lambda_{n-1}
                self.cur_dims.pop(1)

            while dim[0] > 0:
                dim[0] -= 1 # This takes care of the indexing
                b = self.next_state(dim[0])

                # Make sure we have a crystal letter
                ret_crystal_path[-1].append(letters(b)) # Append the rank

            self.cur_dims.pop(0) # Pop off the leading column

        return ret_crystal_path


class MLTToRCBijectionTypeD(KRTToRCBijectionTypeD):
    def run(self):
        r"""
        Run the bijection from a marginally large tableaux to a rigged
        configuration.

        EXAMPLES::

            sage: RC = crystals.infinity.RiggedConfigurations(['D',4])
            sage: T = crystals.infinity.Tableaux(['D',4])
            sage: Psi = T.crystal_morphism({T.module_generators[0]: RC.module_generators[0]})
            sage: TS = [x.value for x in T.subcrystal(max_depth=4)]
            sage: all(Psi(b) == RC(b) for b in TS) # long time # indirect doctest
            True
        """
        for cur_crystal in reversed(self.tp_krt):
            # Iterate through the columns
            cur_column = list(cur_crystal)
            self.cur_path.insert(0, []) # Prepend an empty list

            self.cur_dims.insert(0, [0, 1])

            for letter in reversed(cur_column):
                self.cur_dims[0][0] += 1

                val = letter.value # Convert from a CrystalOfLetter to an Integer

                # Build the next state
                self.cur_path[0].insert(0, [letter]) # Prepend the value
                self.next_state(val)

                if self.cur_dims[0][0] == self.n - 1:
                    # Spinor case, we go from \Lambda_{n-2} -> \Lambda_{n-1} + \Lambda_n
                    self.cur_dims.insert(1, [self.n,1])
                    self.cur_path.insert(1, self.cur_path[0] + [None])

        self.ret_rig_con.set_immutable() # Return it to immutable
        return self.ret_rig_con


class RCToMLTBijectionTypeD(RCToKRTBijectionTypeD):
    def run(self):
        r"""
        Run the bijection from rigged configurations to a marginally large
        tableau.

        EXAMPLES::

            sage: RC = crystals.infinity.RiggedConfigurations(['D',4])
            sage: T = crystals.infinity.Tableaux(['D',4])
            sage: Psi = RC.crystal_morphism({RC.module_generators[0]: T.module_generators[0]})
            sage: RCS = [x.value for x in RC.subcrystal(max_depth=4)]
            sage: all(Psi(nu) == T(nu) for nu in RCS) # long time # indirect doctest
            True
        """
        letters = CrystalOfLetters(self.rigged_con.parent()._cartan_type.classical())
        ret_crystal_path = []

        while self.cur_dims:
            dim = self.cur_dims[0]
            ret_crystal_path.append([])

            # Assumption: all factors are single columns
            if dim[0] == self.n - 1:
                # Spinor case, since we've done \Lambda_n + \Lambda_{n-1} -> \Lambda_{n-2}
                self.cur_dims.pop(1)

            while dim[0] > 0:
                dim[0] -= 1 # This takes care of the indexing
                b = self.next_state(dim[0])

                # Make sure we have a crystal letter
                ret_crystal_path[-1].append(letters(b)) # Append the rank

            self.cur_dims.pop(0) # Pop off the leading column

        return ret_crystal_path
