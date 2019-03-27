.. _sage.coding:

Coding Theory
=============

Basic Coding Theory objects
---------------------------


.. toctree::
   :maxdepth: 1

   sage/coding/linear_code
   sage/coding/channel_constructions
   sage/coding/decoder
   sage/coding/encoder

Catalogs
--------

.. toctree::
   :maxdepth: 2

   sage/coding/channels_catalog
   sage/coding/codes_catalog
   sage/coding/decoders_catalog
   sage/coding/encoders_catalog
   sage/coding/bounds_catalog
   sage/coding/databases
   sage/coding/two_weight_db

Code constructions
------------------

.. toctree::
   :maxdepth: 1

   sage/coding/linear_code

The named code families below are represented in Sage by their own classes,
allowing specialised implementations of e.g. decoding or computation of properties:

.. toctree::
   :maxdepth: 2

   sage/coding/grs
   sage/coding/hamming_code
   sage/coding/golay_code
   sage/coding/parity_check_code
   sage/coding/reed_muller_code
   sage/coding/cyclic_code
   sage/coding/bch

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
is from a special code family, the derived codes inherit e.g. decoding or
minimum distance capabilities:

.. toctree::
   :maxdepth: 2

   sage/coding/subfield_subcode
   sage/coding/punctured_code
   sage/coding/extended_code

Other derived constructions that simply produce the modified generator matrix
can be found among the methods of a constructed code.

Methods and Operations related to Linear Codes
----------------------------------------------


.. toctree::
   :maxdepth: 2

   sage/coding/codecan/codecan
   sage/coding/codecan/autgroup_can_label

Source coding
-------------

.. toctree::
   :maxdepth: 1

   sage/coding/source_coding/huffman

Other modules
-------------

.. toctree::
   :maxdepth: 1

   sage/coding/relative_finite_field_extension
   sage/coding/guruswami_sudan/gs_decoder
   sage/coding/guruswami_sudan/interpolation
   sage/coding/guruswami_sudan/utils
   sage/coding/information_set_decoder
   sage/coding/code_bounds
   sage/coding/delsarte_bounds


.. include:: ../footer.txt
