from cypari2.convert cimport new_gen_from_double

cpdef Gen new_gen_from_real_double_element(RealDoubleElement self):
    """
    Return a PARI representation of ``self``.

    EXAMPLES::

        sage: from sage.libs.pari.convert_sage_real_double import new_gen_from_real_double_element
        sage: new_gen_from_real_double_element(RDF(-2.5))
        -2.50000000000000
    """
    return new_gen_from_double(self._value)
