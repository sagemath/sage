.. _sage.coding:

Coding Theory
=============

Coding theory is the mathematical theory for algebraic and combinatorial codes
used for forward error correction in communications theory. Sage provides an
extensive library of objects and algorithms in coding theory.

Basic objects in coding theory are codes, channels, encoders, and
decoders. The following modules provide the base classes defining them.

.. toctree::
   :maxdepth: 1

   sage/coding/abstract_code
   sage/coding/channel
   sage/coding/encoder
   sage/coding/decoder

Catalogs for available constructions of the basic objects and for bounds on
the parameters of linear codes are provided.

.. toctree::
   :maxdepth: 1

   sage/coding/channels_catalog
   sage/coding/codes_catalog
   sage/coding/decoders_catalog
   sage/coding/encoders_catalog
   sage/coding/bounds_catalog


Linear Codes
------------

The following module is a base class for linear code objects regardless their
metric.

.. toctree::
   :maxdepth: 1

   sage/coding/linear_code_no_metric

There is a number of representatives of linear codes over a specific metric.

.. toctree::
   :maxdepth: 1

   sage/coding/linear_code
   sage/coding/linear_rank_metric

Families of Linear Codes
------------------------

Famous families of codes, listed below, are represented in Sage by their own
classes. For some of them, implementations of special decoding algorithms or
computations for structural invariants are available.

.. toctree::
   :maxdepth: 1

   sage/coding/parity_check_code
   sage/coding/hamming_code
   sage/coding/cyclic_code
   sage/coding/bch_code
   sage/coding/golay_code
   sage/coding/reed_muller_code
   sage/coding/grs_code
   sage/coding/goppa_code
   sage/coding/kasami_codes
   sage/coding/ag_code

.. toctree::
   :hidden:

   sage/coding/ag_code_decoders

In contrast, for some code families Sage can only construct their generator
matrix and has no other a priori knowledge on them:

.. toctree::
   :maxdepth: 1

   sage/coding/code_constructions
   sage/coding/guava
   sage/coding/self_dual_codes
   sage/coding/binary_code

Derived Code Constructions
--------------------------

Sage supports the following derived code constructions. If the constituent code
is from a special code family, the derived codes inherit structural properties
like decoding radius or minimum distance:

.. toctree::
   :maxdepth: 1

   sage/coding/subfield_subcode
   sage/coding/punctured_code
   sage/coding/extended_code

Other derived constructions that simply produce the modified generator matrix
can be found among the methods of a constructed code.

Decoding
--------

Information-set decoding for linear codes:

.. toctree::
   :maxdepth: 1

   sage/coding/information_set_decoder

Guruswami-Sudan interpolation-based list decoding for Reed-Solomon codes:

.. toctree::
   :maxdepth: 1

   sage/coding/guruswami_sudan/gs_decoder
   sage/coding/guruswami_sudan/interpolation
   sage/coding/guruswami_sudan/utils

Automorphism Groups of Linear Codes
-----------------------------------

.. toctree::
   :maxdepth: 1

   sage/coding/codecan/codecan
   sage/coding/codecan/autgroup_can_label

Bounds for Parameters of Linear Codes
-------------------------------------

.. toctree::
   :maxdepth: 1

   sage/coding/code_bounds
   sage/coding/delsarte_bounds

Databases for Coding Theory
---------------------------

.. toctree::
   :maxdepth: 1

   sage/coding/databases
   sage/coding/two_weight_db

Miscellaneous Modules
---------------------

There is at least one module in Sage for source coding in communications theory:

.. toctree::
   :maxdepth: 1

   sage/coding/source_coding/huffman

.. include:: ../footer.txt
