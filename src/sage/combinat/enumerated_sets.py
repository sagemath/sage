"""
Enumerated sets and combinatorial objects

.. TODO:: Proofread / point to the main classes rather than the modules

Categories
----------

- :class:`EnumeratedSets`, :class:`FiniteEnumeratedSets`

Basic enumerated sets
---------------------

- :class:`~sage.combinat.subset.Subsets`, :class:`~sage.combinat.combination.Combinations`
- :class:`~sage.combinat.permutation.Arrangements`, :class:`~sage.combinat.tuple.Tuples`
- :class:`~sage.sets.finite_enumerated_set.FiniteEnumeratedSet`
- :class:`~DisjointUnionEnumeratedSets`

Integer lists
-------------

- :ref:`sage.combinat.partition`
  (see also: :ref:`sage.combinat.catalog_partitions`)
- :ref:`sage.combinat.composition`
- :class:`~sage.combinat.composition_signed.SignedCompositions`
- :class:`IntegerListsLex`
- :ref:`sage.combinat.superpartition`

- :class:`~sage.combinat.integer_vector.IntegerVectors`
- :func:`~sage.combinat.integer_vector_weighted.WeightedIntegerVectors`
- :class:`~sage.combinat.integer_vectors_mod_permgroup.IntegerVectorsModPermutationGroup`

- :ref:`sage.combinat.parking_functions`
- :ref:`sage.combinat.non_decreasing_parking_function`

- :ref:`sage.combinat.sidon_sets`

Words
-----

- :class:`Words`
- :ref:`sage.combinat.subword`
- :ref:`sage.combinat.necklace`
- :ref:`sage.combinat.words.lyndon_word`
- :ref:`sage.combinat.dyck_word`
- :ref:`sage.combinat.debruijn_sequence`
- :ref:`sage.combinat.shuffle`

Permutations, ...
-----------------

- :ref:`sage.combinat.permutation`
- :ref:`sage.combinat.permutation_cython`
- :ref:`sage.combinat.affine_permutation`
- :class:`~sage.combinat.permutation.Arrangements`
- :ref:`sage.combinat.derangements`
- :ref:`sage.combinat.baxter_permutations`

.. SEEALSO::

    - :class:`SymmetricGroup`, :func:`PermutationGroup`, :ref:`sage.groups.perm_gps.permutation_groups_catalog`
    - :class:`FiniteSetMaps`
    - :ref:`sage.combinat.integer_vectors_mod_permgroup`
    - :ref:`sage.combinat.rsk`

Partitions, tableaux, ...
-------------------------

See: :ref:`sage.combinat.catalog_partitions`

Polyominoes
-----------

See: :ref:`sage.combinat.parallelogram_polyomino`

Integer matrices, ...
---------------------

- :ref:`sage.combinat.integer_matrices`
- :ref:`sage.combinat.matrices.hadamard_matrix`
- :ref:`sage.combinat.matrices.latin`
- :ref:`sage.combinat.alternating_sign_matrix`
- :ref:`sage.combinat.six_vertex_model`
- :ref:`sage.combinat.similarity_class_type`
- :ref:`sage.combinat.restricted_growth`
- :ref:`sage.combinat.vector_partition`

.. SEEALSO::

    - :class:`MatrixSpace`
    - :ref:`sage.groups.matrix_gps.catalog`

Subsets and set partitions
--------------------------

- :class:`~sage.combinat.subset.Subsets`, :class:`~sage.combinat.combination.Combinations`
- :class:`~sage.combinat.subsets_pairwise.PairwiseCompatibleSubsets`
- :ref:`sage.combinat.subsets_hereditary`
- :ref:`sage.combinat.set_partition_ordered`
- :ref:`sage.combinat.set_partition`
- :ref:`sage.combinat.diagram_algebras`
- :class:`~sage.combinat.multiset_partition_into_sets_ordered.OrderedMultisetPartitionsIntoSets`,
  :class:`~sage.combinat.multiset_partition_into_sets_ordered.OrderedMultisetPartitionIntoSets`

Trees
-----

- :ref:`sage.combinat.abstract_tree`
- :ref:`sage.combinat.ordered_tree`
- :ref:`sage.combinat.binary_tree`
- :ref:`sage.combinat.rooted_tree`

Enumerated sets related to graphs
---------------------------------

- :ref:`sage.combinat.degree_sequences`
- :ref:`sage.combinat.graph_path`
- :ref:`sage.combinat.perfect_matching`

Backtracking solvers and generic enumerated sets
------------------------------------------------

.. TODO::

    Do we want a separate section, possibly more proeminent, for
    backtracking solvers?

- :func:`~sage.sets.recursively_enumerated_set.RecursivelyEnumeratedSet`
- :class:`~sage.combinat.backtrack.GenericBacktracker`
- :mod:`sage.parallel.map_reduce`
- :ref:`sage.combinat.tiling`
- :ref:`sage.combinat.dlx`
- :ref:`sage.combinat.matrices.dlxcpp`
- :ref:`sage.combinat.species.all`
- :class:`~sage.combinat.integer_lists.IntegerListsLex`
- :class:`~sage.combinat.integer_vectors_mod_permgroup.IntegerVectorsModPermutationGroup`

Low level enumerated sets
-------------------------

- :ref:`sage.combinat.gray_codes`

Misc enumerated sets
--------------------

- :class:`~sage.combinat.gelfand_tsetlin_patterns.GelfandTsetlinPattern`, :class:`~sage.combinat.gelfand_tsetlin_patterns.GelfandTsetlinPatterns`
- :class:`~sage.combinat.knutson_tao_puzzles.KnutsonTaoPuzzleSolver`
- :func:`LatticePolytope`
"""
