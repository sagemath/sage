# sage.doctest: needs sage.libs.flint
r"""
This module is now called ``ring.py`` (see :issue:`31559`). Do not import from here as
it will generate a deprecation warning.

TESTS::

    sage: from sage.modular.modform.ring import find_generators
    sage: find_generators(ModularFormsRing(1))
    doctest:warning
    ...
    DeprecationWarning: find_generators is deprecated. Please use sage.modular.modform.ring.generators instead.
    See https://github.com/sagemath/sage/issues/31559 for details.
    [(4,
    1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + 60480*q^6 + 82560*q^7 + 140400*q^8 + 181680*q^9 + O(q^10)),
    (6,
    1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 - 4058208*q^6 - 8471232*q^7 - 17047800*q^8 - 29883672*q^9 + O(q^10))]

::

    sage: from sage.modular.modform.find_generators import find_generators
    sage: find_generators(ModularFormsRing(1))
    doctest:warning
    ...
    DeprecationWarning: Importing find_generators from here is deprecated; please use "from sage.modular.modform.ring import find_generators" instead.
    See https://github.com/sagemath/sage/issues/31559 for details.
    [(4,
    1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + 60480*q^6 + 82560*q^7 + 140400*q^8 + 181680*q^9 + O(q^10)),
    (6,
    1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 - 4058208*q^6 - 8471232*q^7 - 17047800*q^8 - 29883672*q^9 + O(q^10))]

::

    sage: from sage.modular.modform.ring import basis_for_modform_space
    sage: basis_for_modform_space(ModularFormsRing(1), 4)
    doctest:warning
    ...
    DeprecationWarning: basis_for_modform_space is deprecated. Please use sage.modular.modform.ring.q_expansion_basis instead.
    See https://github.com/sagemath/sage/issues/31559 for details.
    [1 + 240*q + O(q^2)]

::

    sage: from sage.modular.modform.find_generators import _span_of_forms_in_weight
    sage: forms = [(4, 240*eisenstein_series_qexp(4,5)), (6,504*eisenstein_series_qexp(6,5))]
    sage: _span_of_forms_in_weight(forms, 12, prec=5)
    doctest:warning
    ...
    DeprecationWarning: Importing _span_of_forms_in_weight from here is deprecated; please use "from sage.modular.modform.ring import _span_of_forms_in_weight" instead.
    See https://github.com/sagemath/sage/issues/31559 for details.
    Vector space of degree 5 and dimension 2 over Rational Field
    Basis matrix:
    [        1         0    196560  16773120 398034000]
    [        0         1       -24       252     -1472]
"""

from sage.misc.lazy_import import lazy_import
lazy_import('sage.modular.modform.ring', ('_span_of_forms_in_weight', 'find_generators', 'basis_for_modform_space', 'ModularFormsRing'), deprecation=31559)
