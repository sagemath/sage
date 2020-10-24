Asymptotic Expansions
=====================


The Asymptotic Ring
-------------------

The asymptotic ring, as well as its main documentation is contained in
the module

- :doc:`sage/rings/asymptotic/asymptotic_ring`.


Asymptotic Expansion Generators
-------------------------------

Some common asymptotic expansions can be generated in

- :doc:`sage/rings/asymptotic/asymptotic_expansion_generators`.


Supplements
-----------

Behind the scenes of working with asymptotic expressions a couple of
additional classes and tools turn up. For instance the growth of each
summand is managed in growth groups, see below.


Growth Groups
^^^^^^^^^^^^^

The growth of a summand of an asymptotic expression is managed in

- :doc:`sage/rings/asymptotic/growth_group` and

- :doc:`sage/rings/asymptotic/growth_group_cartesian`.


Term Monoids
^^^^^^^^^^^^

A summand of an asymptotic expression is basically a term out of the following monoid:

- :doc:`sage/rings/asymptotic/term_monoid`.


Miscellaneous
^^^^^^^^^^^^^

Various useful functions and tools are collected in

- :doc:`sage/rings/asymptotic/misc`.


Asymptotic Expansions --- Table of Contents
-------------------------------------------

.. toctree::

   sage/rings/asymptotic/asymptotic_ring
   sage/rings/asymptotic/asymptotic_expansion_generators
   sage/rings/asymptotic/growth_group
   sage/rings/asymptotic/growth_group_cartesian
   sage/rings/asymptotic/term_monoid
   sage/rings/asymptotic/misc
   sage/rings/asymptotic/asymptotics_multivariate_generating_functions

.. include:: ../footer.txt
