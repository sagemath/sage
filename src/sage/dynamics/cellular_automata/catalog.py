r"""
Catalog of Cellular Automata

The ``cellular_automata`` object may be used to access examples of various
cellular automata currently implemented in Sage. Using tab-completion on
this object is an easy way to discover and quickly create the cellular
automata that are available (as listed here).

Let ``<tab>`` indicate pressing the :kbd:`Tab` key.  So begin by typing
``cellular_automata.<tab>`` to the see the currently implemented
named cellular automata.

- :class:`cellular_automata.Elementary
  <sage.dynamics.cellular_automata.elementary.ElementaryCellularAutomata>`
- :class:`cellular_automata.GraftalLace
  <sage.dynamics.cellular_automata.glca.GraftalLaceCellularAutomata>`
- :class:`cellular_automata.PeriodicSoliton
  <sage.dynamics.cellular_automata.solitons.PeriodicSolitonCellularAutomata>`
- :class:`cellular_automata.Soliton
  <sage.dynamics.cellular_automata.solitons.SolitonCellularAutomata>`
"""

from sage.misc.lazy_import import lazy_import
lazy_import('sage.dynamics.cellular_automata.elementary',
            'ElementaryCellularAutomata', 'Elementary',)
lazy_import('sage.dynamics.cellular_automata.glca',
            'GraftalLaceCellularAutomata', 'GraftalLace',)
lazy_import('sage.dynamics.cellular_automata.solitons',
            'SolitonCellularAutomata', 'Soliton')
lazy_import('sage.dynamics.cellular_automata.solitons',
            'PeriodicSolitonCellularAutomata', 'PeriodicSoliton')

del lazy_import # We remove the object from here so it doesn't appear under tab completion
