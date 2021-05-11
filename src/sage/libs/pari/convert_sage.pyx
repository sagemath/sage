r"""
Convert PARI objects to Sage types
"""

#*****************************************************************************
#       Copyright (C) 2016 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#       Copyright (C) 2016 Luca De Feo <luca.defeo@polytechnique.edu>
#       Copyright (C) 2016 Vincent Delecroix <vincent.delecroix@u-bordeaux.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cypari2.types cimport (GEN, typ, t_INT, t_FRAC, t_REAL, t_COMPLEX,
                            t_INTMOD, t_PADIC, t_INFINITY, t_VEC, t_COL,
                            t_VECSMALL, t_MAT, t_STR,
                            lg, precp)
from cypari2.pari_instance cimport prec_words_to_bits
from cypari2.paridecl cimport gel, inf_get_sign

from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.all import RealField, ComplexField, QuadraticField
from sage.matrix.args cimport MatrixArgs
from sage.rings.padics.factory import Qp
from sage.rings.infinity import Infinity


cpdef gen_to_sage(Gen z, locals=None):
    """
    Convert a PARI gen to a Sage/Python object.

    INPUT:

    - ``z`` -- PARI ``gen``

    - ``locals`` -- optional dictionary used in fallback cases that
      involve :func:`sage_eval`

    OUTPUT:

    One of the following depending on the PARI type of ``z``

    - a :class:`~sage.rings.integer.Integer` if ``z`` is an integer (type ``t_INT``)

    - a :class:`~sage.rings.rational.Rational` if ``z`` is a rational (type ``t_FRAC``)

    - a :class:`~sage.rings.real_mpfr.RealNumber` if ``z`` is a real
      number (type ``t_REAL``). The precision will be equivalent.

    - a :class:`~sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic`
      or a :class:`~sage.rings.complex_mpfr.ComplexNumber` if ``z`` is a complex
      number (type ``t_COMPLEX``). The former is used when the real and imaginary parts are
      integers or rationals and the latter when they are floating point numbers. In that
      case The precision will be the maximal precision of the real and imaginary parts.

    - a Python list if ``z`` is a vector or a list (type ``t_VEC``, ``t_COL``)

    - a Python string if ``z`` is a string (type ``t_STR``)

    - a Python list of Python integers if ``z`` is a small vector (type ``t_VECSMALL``)

    - a matrix if ``z`` is a matrix (type ``t_MAT``)

    - a padic element (type ``t_PADIC``)

    - a :class:`~sage.rings.infinity.Infinity` if ``z`` is an infinity
      (type ``t_INF``)

    EXAMPLES::

        sage: from sage.libs.pari.convert_sage import gen_to_sage

    Converting an integer::

        sage: z = pari('12'); z
        12
        sage: z.type()
        't_INT'
        sage: a = gen_to_sage(z); a
        12
        sage: a.parent()
        Integer Ring

        sage: gen_to_sage(pari('7^42'))
        311973482284542371301330321821976049

    Converting a rational number::

        sage: z = pari('389/17'); z
        389/17
        sage: z.type()
        't_FRAC'
        sage: a = gen_to_sage(z); a
        389/17
        sage: a.parent()
        Rational Field

        sage: gen_to_sage(pari('5^30 / 3^50'))
        931322574615478515625/717897987691852588770249

    Converting a real number::

        sage: pari.set_real_precision(70)
        15
        sage: z = pari('1.234'); z
        1.234000000000000000000000000000000000000000000000000000000000000000000
        sage: a = gen_to_sage(z); a
        1.234000000000000000000000000000000000000000000000000000000000000000000000000
        sage: a.parent()
        Real Field with 256 bits of precision
        sage: pari.set_real_precision(15)
        70
        sage: a = gen_to_sage(pari('1.234')); a
        1.23400000000000000
        sage: a.parent()
        Real Field with 64 bits of precision

    For complex numbers, the parent depends on the PARI type::

        sage: z = pari('(3+I)'); z
        3 + I
        sage: z.type()
        't_COMPLEX'
        sage: a = gen_to_sage(z); a
        i + 3
        sage: a.parent()
        Number Field in i with defining polynomial x^2 + 1 with i = 1*I

        sage: z = pari('(3+I)/2'); z
        3/2 + 1/2*I
        sage: a = gen_to_sage(z); a
        1/2*i + 3/2
        sage: a.parent()
        Number Field in i with defining polynomial x^2 + 1 with i = 1*I

        sage: z = pari('1.0 + 2.0*I'); z
        1.00000000000000 + 2.00000000000000*I
        sage: a = gen_to_sage(z); a
        1.00000000000000000 + 2.00000000000000000*I
        sage: a.parent()
        Complex Field with 64 bits of precision

        sage: z = pari('1 + 1.0*I'); z
        1 + 1.00000000000000*I
        sage: a = gen_to_sage(z); a
        1.00000000000000000 + 1.00000000000000000*I
        sage: a.parent()
        Complex Field with 64 bits of precision

        sage: z = pari('1.0 + 1*I'); z
        1.00000000000000 + I
        sage: a = gen_to_sage(z); a
        1.00000000000000000 + 1.00000000000000000*I
        sage: a.parent()
        Complex Field with 64 bits of precision

    Converting polynomials::

        sage: f = pari('(2/3)*x^3 + x - 5/7 + y')
        sage: f.type()
        't_POL'

        sage: R.<x,y> = QQ[]
        sage: gen_to_sage(f, {'x': x, 'y': y})
        2/3*x^3 + x + y - 5/7
        sage: parent(gen_to_sage(f, {'x': x, 'y': y}))
        Multivariate Polynomial Ring in x, y over Rational Field

        sage: x,y = SR.var('x,y')
        sage: gen_to_sage(f, {'x': x, 'y': y})
        2/3*x^3 + x + y - 5/7
        sage: parent(gen_to_sage(f, {'x': x, 'y': y}))
        Symbolic Ring

        sage: gen_to_sage(f)
        Traceback (most recent call last):
        ...
        NameError: name 'x' is not defined

    Converting vectors::

        sage: z1 = pari('[-3, 2.1, 1+I]'); z1
        [-3, 2.10000000000000, 1 + I]
        sage: z2 = pari('[1.0*I, [1,2]]~'); z2
        [1.00000000000000*I, [1, 2]]~
        sage: z1.type(), z2.type()
        ('t_VEC', 't_COL')
        sage: a1 = gen_to_sage(z1)
        sage: a2 = gen_to_sage(z2)
        sage: type(a1), type(a2)
        (<... 'list'>, <... 'list'>)
        sage: [parent(b) for b in a1]
        [Integer Ring,
         Real Field with 64 bits of precision,
         Number Field in i with defining polynomial x^2 + 1 with i = 1*I]
        sage: [parent(b) for b in a2]
        [Complex Field with 64 bits of precision, <... 'list'>]

        sage: z = pari('Vecsmall([1,2,3,4])')
        sage: z.type()
        't_VECSMALL'
        sage: a = gen_to_sage(z); a
        [1, 2, 3, 4]
        sage: type(a)
        <... 'list'>
        sage: [parent(b) for b in a]
        [<... 'int'>, <... 'int'>, <... 'int'>, <... 'int'>]

    Matrices::

        sage: z = pari('[1,2;3,4]')
        sage: z.type()
        't_MAT'
        sage: a = gen_to_sage(z); a
        [1 2]
        [3 4]
        sage: a.parent()
        Full MatrixSpace of 2 by 2 dense matrices over Integer Ring

    Conversion of p-adics::

        sage: z = pari('569 + O(7^8)'); z
        2 + 4*7 + 4*7^2 + 7^3 + O(7^8)
        sage: a = gen_to_sage(z); a
        2 + 4*7 + 4*7^2 + 7^3 + O(7^8)
        sage: a.parent()
        7-adic Field with capped relative precision 8

    Conversion of infinities::

        sage: gen_to_sage(pari('oo'))
        +Infinity
        sage: gen_to_sage(pari('-oo'))
        -Infinity

    Conversion of strings::

        sage: s = pari('"foo"').sage(); s
        'foo'
        sage: type(s)
        <type 'str'>
    """
    cdef GEN g = z.g
    cdef long t = typ(g)
    cdef long tx, ty
    cdef Gen real, imag, prec, xprec, yprec
    cdef Py_ssize_t i, j, nr, nc

    if t == t_INT:
        return Integer(z)
    elif t == t_FRAC:
        return Rational(z)
    elif t == t_REAL:
        prec = z.bitprecision()
        if typ(prec.g) == t_INFINITY:
            sage_prec = 53
        else:
            sage_prec = prec
        return RealField(sage_prec)(z)
    elif t == t_COMPLEX:
        real = z.real()
        imag = z.imag()
        tx = typ(real.g)
        ty = typ(imag.g)
        if tx in [t_INTMOD, t_PADIC] or ty in [t_INTMOD, t_PADIC]:
            raise NotImplementedError("No conversion to python available for t_COMPLEX with t_INTMOD or t_PADIC components")
        if tx == t_REAL or ty == t_REAL:
            xprec = real.bitprecision()  # will be infinite if exact
            yprec = imag.bitprecision()  # will be infinite if exact
            if typ(xprec.g) == t_INFINITY:
                if typ(yprec.g) == t_INFINITY:
                    sage_prec = 53
                else:
                    sage_prec = yprec
            elif typ(yprec.g) == t_INFINITY:
                sage_prec = xprec
            else:
                sage_prec = max(xprec, yprec)

            R = RealField(sage_prec)
            C = ComplexField(sage_prec)
            return C(R(real), R(imag))
        else:
            K = QuadraticField(-1, 'i')
            return K([gen_to_sage(real), gen_to_sage(imag)])
    elif t == t_VEC or t == t_COL:
        return [gen_to_sage(x, locals) for x in z.python_list()]
    elif t == t_VECSMALL:
        return z.python_list_small()
    elif t == t_MAT:
        nc = lg(g) - 1
        nr = 0 if nc == 0 else lg(gel(g,1)) - 1
        ma = MatrixArgs.__new__(MatrixArgs)
        ma.nrows = nr
        ma.ncols = nc
        ma.entries = [gen_to_sage(z[i,j], locals) for i in range(nr) for j in range(nc)]
        return ma.matrix()
    elif t == t_PADIC:
        p = z.padicprime()
        K = Qp(Integer(p), precp(g))
        return K(z.lift())
    elif t == t_STR:
        return str(z)
    elif t == t_INFINITY:
        if inf_get_sign(g) >= 0:
            return Infinity
        else:
            return -Infinity

    # Fallback (e.g. polynomials): use string representation
    from sage.misc.sage_eval import sage_eval
    locals = {} if locals is None else locals
    return sage_eval(str(z), locals=locals)
