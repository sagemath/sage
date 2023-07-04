database_mutation_class: Database of exceptional mutation classes of quivers
============================================================================

Description
-----------

Contains a database of all exceptional mutation classes of quivers.

Every file in the database is of the form
``mutation_classes_n.dig6`` for some ``n`` and

-  contains a ``cPickle.dump`` of a dictionary where
-  the keys are tuples representing irreducible exceptional quiver
   mutation types of rank ``n``, and
-  the values are all quivers in the given mutation class stored in
   canonical form as ``(dig6,edges)`` where
-  ``dig6`` is the dig6 data of the given ``DiGraph``, and
-  ``edges`` are the non-simply-laced edges thereof.
-  is obtained by running the function

   ``sage.combinat.cluster_algebra_quiver.quiver_mutation_type._save_data_dig6(n, types='Exceptional', verbose=False)``


SPKG Maintainers
----------------

-  C. Stump <christian.stump@gmail.com>
