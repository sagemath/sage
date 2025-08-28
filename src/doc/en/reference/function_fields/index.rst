Algebraic Function Fields
=========================

Sage allows basic computations with elements and ideals in orders of algebraic
function fields over arbitrary constant fields. Advanced computations, like
computing the genus or a basis of the Riemann-Roch space of a divisor, are
available for function fields over finite fields, number fields, and the
algebraic closure of `\QQ`.

.. toctree::
   :maxdepth: 1

   sage/rings/function_field/function_field
   sage/rings/function_field/function_field_rational
   sage/rings/function_field/function_field_polymod
   sage/rings/function_field/element
   sage/rings/function_field/element_rational
   sage/rings/function_field/element_polymod
   sage/rings/function_field/order
   sage/rings/function_field/order_rational
   sage/rings/function_field/order_basis
   sage/rings/function_field/order_polymod
   sage/rings/function_field/ideal
   sage/rings/function_field/ideal_rational
   sage/rings/function_field/ideal_polymod
   sage/rings/function_field/place
   sage/rings/function_field/place_rational
   sage/rings/function_field/place_polymod
   sage/rings/function_field/divisor
   sage/rings/function_field/differential
   sage/rings/function_field/valuation_ring
   sage/rings/function_field/derivations
   sage/rings/function_field/derivations_rational
   sage/rings/function_field/derivations_polymod
   sage/rings/function_field/maps
   sage/rings/function_field/extensions
   sage/rings/function_field/constructor

A basic reference for the theory of algebraic function fields is [Stich2009]_.

Jacobians of function fields
----------------------------

Arithmetic in Jacobians of function fields are available in two flavors.

.. toctree::
   :maxdepth: 1

   sage/rings/function_field/jacobian_base
   sage/rings/function_field/jacobian_hess
   sage/rings/function_field/jacobian_khuri_makdisi
   sage/rings/function_field/khuri_makdisi

A Support Module
----------------

.. toctree::
   :maxdepth: 1

   sage/rings/function_field/hermite_form_polynomial

.. include:: ../footer.txt
