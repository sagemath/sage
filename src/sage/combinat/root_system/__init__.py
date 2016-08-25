"""
Root Systems
============

Quickref
--------

- ``T = CartanType(["A", 3]), T.is_finite()``     -- Cartan types
- ``T.dynkin_diagram(), DynkinDiagram(["G",2])``  -- Dynkin diagrams
- ``T.cartan_matrix(),  CartanMatrix(["F",4])``   -- Cartan matrices
- ``RootSystem(T).weight_lattice()``              -- Root systems
- ``WeylGroup(["B", 6, 1]).simple_reflections()`` -- Affine Weyl groups
- ``WeylCharacterRing(["D", 4])``                 -- Weyl character rings

Introductory material
---------------------

- :ref:`sage.combinat.root_system`                -- This overview
- :class:`CartanType`                             -- An introduction to Cartan types
- :class:`RootSystem`                             -- An introduction to root systems
- :ref:`sage.combinat.root_system.plot`           -- A root system visualization tutorial

- The `Lie Methods and Related Combinatorics <../../../../../thematic_tutorials/lie.html>`_ thematic tutorial


Related material
----------------

- :ref:`sage.combinat.crystals`                   -- Crystals

Cartan datum
------------

- :ref:`sage.combinat.root_system.cartan_type`
- :ref:`sage.combinat.root_system.dynkin_diagram`
- :ref:`sage.combinat.root_system.cartan_matrix`
- :ref:`sage.combinat.root_system.coxeter_matrix`
- :ref:`sage.combinat.root_system.coxeter_type`

Root systems
------------

- :ref:`sage.combinat.root_system.root_system`
- :ref:`sage.combinat.root_system.plot`
- :ref:`sage.combinat.root_system.root_lattice_realizations`
- :ref:`sage.combinat.root_system.root_lattice_realization_algebras`
- :ref:`sage.combinat.root_system.weight_lattice_realizations`
- :ref:`sage.combinat.root_system.root_space`
- :ref:`sage.combinat.root_system.weight_space`
- :ref:`sage.combinat.root_system.ambient_space`

Coxeter groups
--------------

- :ref:`sage.combinat.root_system.coxeter_group`
- :ref:`sage.combinat.root_system.weyl_group`
- :ref:`sage.combinat.root_system.extended_affine_weyl_group`
- :ref:`sage.combinat.root_system.fundamental_group`

.. SEEALSO::

    The categories :class:`CoxeterGroups` and :class:`WeylGroups`

Finite reflection groups
------------------------

- :ref:`sage.combinat.root_system.reflection_group_complex`
- :ref:`sage.combinat.root_system.reflection_group_real`

.. SEEALSO::

    The category :class:`~sage.categories.complex_reflection_groups.ComplexReflectionGroups`

Representation theory
---------------------

- :ref:`sage.combinat.root_system.weyl_characters`
- :ref:`sage.combinat.root_system.integrable_representations`
- :ref:`sage.combinat.root_system.branching_rules`
- :ref:`sage.combinat.root_system.hecke_algebra_representation`
- :ref:`sage.combinat.root_system.non_symmetric_macdonald_polynomials`

Root system data and code for specific families of Cartan types
---------------------------------------------------------------

- :ref:`sage.combinat.root_system.type_affine`
- :ref:`sage.combinat.root_system.type_dual`
- :ref:`sage.combinat.root_system.type_folded`
- :ref:`sage.combinat.root_system.type_reducible`
- :ref:`sage.combinat.root_system.type_relabel`

Root system data and code for specific Cartan types
---------------------------------------------------

- :ref:`sage.combinat.root_system.type_A`
- :ref:`sage.combinat.root_system.type_B`
- :ref:`sage.combinat.root_system.type_C`
- :ref:`sage.combinat.root_system.type_D`
- :ref:`sage.combinat.root_system.type_E`
- :ref:`sage.combinat.root_system.type_F`
- :ref:`sage.combinat.root_system.type_G`
- :ref:`sage.combinat.root_system.type_H`
- :ref:`sage.combinat.root_system.type_I`
- :ref:`sage.combinat.root_system.type_A_affine`
- :ref:`sage.combinat.root_system.type_B_affine`
- :ref:`sage.combinat.root_system.type_C_affine`
- :ref:`sage.combinat.root_system.type_D_affine`
- :ref:`sage.combinat.root_system.type_E_affine`
- :ref:`sage.combinat.root_system.type_F_affine`
- :ref:`sage.combinat.root_system.type_G_affine`
- :ref:`sage.combinat.root_system.type_BC_affine`
- :ref:`sage.combinat.root_system.type_A_infinity`
"""
# currently needed to activate the backward compatibility
# register_unpickle_override
from . import type_A
from . import type_B
from . import type_C
from . import type_D
from . import type_E
from . import type_F
from . import type_G
