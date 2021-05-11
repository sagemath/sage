"""
Empty Species
"""
#*****************************************************************************
#       Copyright (C) 2008 Florent Hivert <Florent.Hivert@univ-rouen,fr>,
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
from .species import GenericCombinatorialSpecies
from .series_order import inf
from sage.structure.unique_representation import UniqueRepresentation

class EmptySpecies(GenericCombinatorialSpecies, UniqueRepresentation):
    """
    Returns the empty species. This species has no structure at all.
    It is the zero of the semi-ring of species.

    EXAMPLES::

        sage: X = species.EmptySpecies(); X
        Empty species
        sage: X.structures([]).list()
        []
        sage: X.structures([1]).list()
        []
        sage: X.structures([1,2]).list()
        []
        sage: X.generating_series().coefficients(4)
        [0, 0, 0, 0]
        sage: X.isotype_generating_series().coefficients(4)
        [0, 0, 0, 0]
        sage: X.cycle_index_series().coefficients(4)
        [0, 0, 0, 0]

    The empty species is the zero of the semi-ring of species.
    The following tests that it is neutral with respect to addition::

        sage: Empt  = species.EmptySpecies()
        sage: S = species.CharacteristicSpecies(2)
        sage: X = S + Empt
        sage: X == S    # TODO: Not Implemented
        True
        sage: (X.generating_series().coefficients(4) ==
        ....:  S.generating_series().coefficients(4))
        True
        sage: (X.isotype_generating_series().coefficients(4) ==
        ....:  S.isotype_generating_series().coefficients(4))
        True
        sage: (X.cycle_index_series().coefficients(4) ==
        ....:  S.cycle_index_series().coefficients(4))
        True

    The following tests that it is the zero element with respect to
    multiplication::

        sage: Y = Empt*S
        sage: Y == Empt   # TODO: Not Implemented
        True
        sage: Y.generating_series().coefficients(4)
        [0, 0, 0, 0]
        sage: Y.isotype_generating_series().coefficients(4)
        [0, 0, 0, 0]
        sage: Y.cycle_index_series().coefficients(4)
        [0, 0, 0, 0]

    TESTS::

        sage: Empt  = species.EmptySpecies()
        sage: Empt2 = species.EmptySpecies()
        sage: Empt is Empt2
        True
    """
    def __init__(self, min=None, max=None, weight=None):
        """
        Initializer for the empty species.

        EXAMPLES::

            sage: F = species.EmptySpecies()
            sage: F._check()
            True
            sage: F == loads(dumps(F))
            True
        """
        # There is no structure at all, so we set min and max accordingly.
        GenericCombinatorialSpecies.__init__(self, weight=weight)
        self._name = "Empty species"

    def _gs(self, series_ring, base_ring):
        """
        Return the generating series for self.

        EXAMPLES::

            sage: F = species.EmptySpecies()
            sage: F.generating_series().coefficients(5) # indirect doctest
            [0, 0, 0, 0, 0]
            sage: F.generating_series().count(3)
            0
            sage: F.generating_series().count(4)
            0
        """
        return series_ring.zero()

    _itgs = _gs
    _cis  = _gs

    def _order(self):
        """
        Returns the order of the generating series.

        EXAMPLES::

            sage: F = species.EmptySpecies()
            sage: F._order()
            Infinite series order
        """
        return inf

    def _structures(self, structure_class, labels):
        """
        Thanks to the counting optimisation, this is never called... Otherwise
        this should return an empty iterator.

        EXAMPLES::

            sage: F = species.EmptySpecies()
            sage: F.structures([]).list()      # indirect doctest
            []
            sage: F.structures([1,2,3]).list() # indirect doctest
            []
        """
        assert False, "This should never be called"

    _default_structure_class = 0
    _isotypes = _structures

    def _equation(self, var_mapping):
        """
        Returns the right hand side of an algebraic equation satisfied by
        this species. This is a utility function called by the
        algebraic_equation_system method.

        EXAMPLES::

            sage: C = species.EmptySpecies()
            sage: Qz = QQ['z']
            sage: R.<node0> = Qz[]
            sage: var_mapping = {'z':Qz.gen(), 'node0':R.gen()}
            sage: C._equation(var_mapping)
            0
        """
        return 0

EmptySpecies_class = EmptySpecies
