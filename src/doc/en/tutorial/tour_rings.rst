.. _section-rings:

Basic Rings
===========

When defining matrices, vectors, or polynomials, it is sometimes
useful and sometimes necessary to specify the "ring" over which it is
defined.  A *ring* is a mathematical construction in which there are
well-behaved notions of addition and multiplication; if you've never
heard of them before, you probably just need to know about these
four commonly used rings:

* the integers `\{..., -1, 0, 1, 2, ...\}`, called ``ZZ`` in Sage.
* the rational numbers -- i.e., fractions, or ratios, of integers --
  called ``QQ`` in Sage.
* the real numbers, called ``RR`` in Sage.
* the complex numbers, called ``CC`` in Sage.

You may need to know about these distinctions because the same
polynomial, for example, can be treated differently depending on the
ring over which it is defined.  For instance, the polynomial `x^2-2`
has two roots, `\pm \sqrt{2}`.  Those roots are not rational, so if
you are working with polynomials with rational coefficients, the
polynomial won't factor.  With real coefficients, it will.  Therefore
you may want to specify the ring to insure that you are getting the
information you expect.  The following two commands defines the sets
of polynomials with rational coefficients and real coefficients,
respectively.  The sets are named "ratpoly" and "realpoly", but these
aren't important here; however, note that the strings ".<t>" and
".<z>" name the *variables* used in the two cases. ::

    sage: ratpoly.<t> = PolynomialRing(QQ)
    sage: realpoly.<z> = PolynomialRing(RR)

Now we illustrate the point about factoring `x^2-2`:

.. link

::

    sage: factor(t^2-2)
    t^2 - 2
    sage: factor(z^2-2)
    (z - 1.41421356237310) * (z + 1.41421356237310)

Similar comments apply to matrices: the row-reduced form of a matrix
can depend on the ring over which it is defined, as can its
eigenvalues and eigenvectors.  For more about constructing
polynomials, see :ref:`section-poly`, and for more about matrices, see
:ref:`section-linalg`.

The symbol ``I`` represents the square root of :math:`-1`; ``i`` is a
synonym for ``I``. Of course, this is not a rational number::

    sage: i  # square root of -1
    I
    sage: i in QQ
    False

Note: The above code may not work as expected if the variable ``i``
has been assigned a different value, for example, if it was used
as a loop variable. If this is the case, type ::

    sage: reset('i')

to get the original complex value of ``i``.

There is one subtlety in defining complex numbers: as mentioned above,
the symbol ``i`` represents a square root of `-1`, but it is a
*formal* square root of `-1` as an algebraic number.  Calling ``CC(i)``
or ``CC.0`` or ``CC.gen(0)`` returns the *complex* square root of `-1`.
Arithmetic involving different kinds of numbers is possible by
so-called coercion, see :ref:`section-coercion`.

::

    sage: i = CC(i)       # floating point complex number
    sage: i == CC.0
    True
    sage: a, b = 4/3, 2/3
    sage: z = a + b*i
    sage: z
    1.33333333333333 + 0.666666666666667*I
    sage: z.imag()        # imaginary part
    0.666666666666667
    sage: z.real() == a   # automatic coercion before comparison
    True
    sage: a + b
    2
    sage: 2*b == a
    True
    sage: parent(2/3)
    Rational Field
    sage: parent(4/2)
    Rational Field
    sage: 2/3 + 0.1       # automatic coercion before addition
    0.766666666666667
    sage: 0.1 + 2/3       # coercion rules are symmetric in Sage
    0.766666666666667

Here are more examples of basic rings in Sage. As noted above, the
ring of rational numbers may be referred to using ``QQ``, or also
``RationalField()`` (a *field* is a ring in
which the multiplication is commutative and in which every nonzero
element has a reciprocal in that ring, so the rationals form a field,
but the integers don't)::

    sage: RationalField()
    Rational Field
    sage: QQ
    Rational Field
    sage: 1/2 in QQ
    True

The decimal number ``1.2`` is considered to be in ``QQ``: decimal numbers
which happen to also be rational can be "coerced" into the rational
numbers (see :ref:`section-coercion`).  The numbers `\pi` and `\sqrt{2}`
are not rational, though::

    sage: 1.2 in QQ
    True
    sage: pi in QQ
    False
    sage: pi in RR
    True
    sage: sqrt(2) in QQ
    False
    sage: sqrt(2) in CC
    True

For use in higher mathematics, Sage also knows about other rings, such
as finite fields, `p`-adic integers, the ring of algebraic numbers,
polynomial rings, and matrix rings.  Here are constructions of some of
these::

    sage: GF(3)
    Finite Field of size 3
    sage: GF(27, 'a')  # need to name the generator if not a prime field
    Finite Field in a of size 3^3
    sage: Zp(5)
    5-adic Ring with capped relative precision 20
    sage: sqrt(3) in QQbar # algebraic closure of QQ
    True
