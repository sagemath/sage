"""
QQbar decorators

Python decorators for use with the algebraic field QQbar.

AUTHORS:

- Brent Baccala (7 Jun 2018) -- handle_AA_and_QQbar

Decorators
==========
"""

from sage.misc.decorators import decorator_keywords, sage_wraps

@decorator_keywords
def handle_AA_and_QQbar(func):
    r"""
    Decorator to call a function that only accepts arguments in number fields.

    The argument list is scanned for ideals and/or polynomials over algebraic
    fields (``QQbar`` or ``AA``).  If any exist, they are converted to a common
    number field before calling the function, and the results are converted back.
    Lists, dictionaries (values only), sets, and tuples are converted recursively.

    This decorator can not used with methods that depend on factoring, since
    factorization might require larger number fields than those required to
    express the polynomials.  No means is provided to check whether factoring
    is being attempted by a wrapped method, and if a method invoked a library
    or subprocess (like Singular), it's hard to imagine how such a check could
    be performed.

    See https://mathoverflow.net/questions/304525 for a discussion of why a
    simple attempt to overcome this limitation didn't work.
    """

    @sage_wraps(func)
    def wrapper(*args, **kwds):

        """
        TESTS::

            sage: from sage.rings.qqbar_decorators import handle_AA_and_QQbar
            sage: @handle_AA_and_QQbar
            ....: def return_base_ring(x):
            ....:     return x.base_ring()

            sage: P.<x> = QQbar[]
            sage: return_base_ring(x)
            Rational Field

            sage: P.<y,z> = QQbar[]
            sage: return_base_ring(y)
            Rational Field

            sage: return_base_ring(ideal(y,z))
            Rational Field

        Check that :trac:`29468` is fixed::

            sage: J = QQbar['x,y'].ideal('x^2 - y')
            sage: type(J.groebner_basis())
            <class 'sage.rings.polynomial.multi_polynomial_sequence.PolynomialSequence_generic'>
            sage: J.groebner_basis().is_immutable()
            True

        ::

            sage: @handle_AA_and_QQbar
            ....: def f(x):
            ....:     print(x.ring().base_ring())
            ....:     return x
            sage: R.<x,y> = QQbar[]
            sage: s = Sequence([x, R(sqrt(2)) * y], immutable=True)
            sage: t = f(s)
            Number Field in a with defining polynomial y^2 - 2
            sage: t.ring().base_ring()
            Algebraic Field
            sage: t.is_immutable()
            True
            sage: s == t
            True
        """

        from sage.misc.flatten import flatten
        from sage.rings.polynomial.polynomial_element import Polynomial
        from sage.rings.polynomial.multi_polynomial import MPolynomial
        from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence, is_PolynomialSequence
        from sage.rings.ideal import Ideal, Ideal_generic
        from sage.rings.qqbar import AlgebraicField_common, number_field_elements_from_algebraics

        if not any(isinstance(a, (Polynomial, MPolynomial, Ideal_generic))
                   and isinstance(a.base_ring(), AlgebraicField_common)
                   or is_PolynomialSequence(a)
                   and isinstance(a.ring().base_ring(), AlgebraicField_common) for a in args):
            return func(*args, **kwds)

        polynomials = []

        for a in flatten(args, ltypes=(list, tuple, set)):
            if isinstance(a, Ideal_generic):
                polynomials.extend(a.gens())
            elif isinstance(a, Polynomial):
                polynomials.append(a)
            elif isinstance(a, MPolynomial):
                polynomials.append(a)

        orig_elems = flatten([p.coefficients() for p in polynomials])

        # We need minimal=True if these elements are over AA, because
        # same_field=True might trigger an exception otherwise.

        numfield, new_elems, morphism = number_field_elements_from_algebraics(orig_elems, same_field=True, minimal=True)

        elem_dict = dict(zip(orig_elems, new_elems))

        def forward_map(item):
            if isinstance(item, Ideal_generic):
                return Ideal([forward_map(g) for g in item.gens()])
            elif isinstance(item, Polynomial):
                return item.map_coefficients(elem_dict.__getitem__, new_base_ring=numfield)
            elif isinstance(item, MPolynomial):
                return item.map_coefficients(elem_dict.__getitem__, new_base_ring=numfield)
            elif is_PolynomialSequence(item):
                return PolynomialSequence(map(forward_map, item),
                                          immutable=item.is_immutable())
            elif isinstance(item, list):
                return list(map(forward_map, item))
            elif isinstance(item, dict):
                return {k: forward_map(v) for k,v in item.items()}
            elif isinstance(item, tuple):
                return tuple(map(forward_map, item))
            elif isinstance(item, set):
                return set(map(forward_map, list(item)))
            else:
                return item

        def reverse_map(item):
            if isinstance(item, Ideal_generic):
                return Ideal([reverse_map(g) for g in item.gens()])
            elif isinstance(item, Polynomial):
                return item.map_coefficients(morphism)
            elif isinstance(item, MPolynomial):
                return item.map_coefficients(morphism)
            elif is_PolynomialSequence(item):
                return PolynomialSequence(map(reverse_map, item),
                                          immutable=item.is_immutable())
            elif isinstance(item, list):
                return list(map(reverse_map, item))
            elif isinstance(item, tuple):
                return tuple(map(reverse_map, item))
            elif isinstance(item, set):
                return set(map(reverse_map, list(item)))
            else:
                return item

        args = forward_map(args)
        kwds = forward_map(kwds)

        return reverse_map(func(*args, **kwds))

    return wrapper
