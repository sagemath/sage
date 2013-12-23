from sage.rings.all import *

def pari(x):
    """
    Return the pari object constructed from a Sage object.

    The work is done by the __call__ method of the class PariInstance,
    which in turn passes the work to any class which has its own
    method _pari_().

    EXAMPLES::

        sage: pari([2,3,5])
        [2, 3, 5]
        sage: pari(Matrix(2,2,range(4)))
        [0, 1; 2, 3]
        sage: pari(x^2-3)
        x^2 - 3

    ::

        sage: a = pari(1); a, a.type()
        (1, 't_INT')
        sage: a = pari(1/2); a, a.type()
        (1/2, 't_FRAC')
        sage: a = pari(1/2); a, a.type()
        (1/2, 't_FRAC')

    Conversion from reals uses the real's own precision::

        sage: a = pari(1.2); a, a.type(), a.precision()
        (1.20000000000000, 't_REAL', 4) # 32-bit
        (1.20000000000000, 't_REAL', 3) # 64-bit

    Conversion from strings uses the current PARI real precision.
    By default, this is 96 bits on 32-bit systems and 128 bits on
    64-bit systems::

        sage: a = pari('1.2'); a, a.type(), a.precision()
        (1.20000000000000, 't_REAL', 5) # 32-bit
        (1.20000000000000, 't_REAL', 4) # 64-bit

    But we can change this precision::

        sage: pari.set_real_precision(35)  # precision in decimal digits
        15
        sage: a = pari('1.2'); a, a.type(), a.precision()
        (1.2000000000000000000000000000000000, 't_REAL', 6)  # 32-bit
        (1.2000000000000000000000000000000000, 't_REAL', 4)  # 64-bit

    Set the precision to 15 digits for the remaining tests::

        sage: pari.set_real_precision(15)
        35

    Conversion from matrices is supported, but not from vectors;
    use lists or tuples instead::

        sage: a = pari(matrix(2,3,[1,2,3,4,5,6])); a, a.type()
        ([1, 2, 3; 4, 5, 6], 't_MAT')

        sage: v = vector([1.2,3.4,5.6])
        sage: pari(v)
        Traceback (most recent call last):
        ...
        PariError: syntax error, unexpected ')', expecting )-> or ','
        sage: b = pari(list(v)); b,b.type()
        ([1.20000000000000, 3.40000000000000, 5.60000000000000], 't_VEC')
        sage: b = pari(tuple(v)); b, b.type()
        ([1.20000000000000, 3.40000000000000, 5.60000000000000], 't_VEC')

    Some more exotic examples::

        sage: K.<a> = NumberField(x^3 - 2)
        sage: pari(K)
        [y^3 - 2, [1, 1], -108, 1, [[1, 1.25992104989487, 1.58740105196820; 1, -0.629960524947437 - 1.09112363597172*I, -0.793700525984100 + 1.37472963699860*I], [1, 1.25992104989487, 1.58740105196820; 1, -1.72108416091916, 0.581029111014503; 1, 0.461163111024285, -2.16843016298270], [1, 1, 2; 1, -2, 1; 1, 0, -2], [3, 0, 0; 0, 0, 6; 0, 6, 0], [6, 0, 0; 0, 6, 0; 0, 0, 3], [2, 0, 0; 0, 0, 1; 0, 1, 0], [2, [0, 0, 2; 1, 0, 0; 0, 1, 0]]], [1.25992104989487, -0.629960524947437 - 1.09112363597172*I], [1, y, y^2], [1, 0, 0; 0, 1, 0; 0, 0, 1], [1, 0, 0, 0, 0, 2, 0, 2, 0; 0, 1, 0, 1, 0, 0, 0, 0, 2; 0, 0, 1, 0, 1, 0, 1, 0, 0]]

        sage: E = EllipticCurve('37a1')
        sage: pari(E)
        [0, 0, 1, -1, 0, 0, -2, 1, -1, 48, -216, 37, 110592/37, [0.837565435283323, 0.269594436405445, -1.10715987168877]~, 2.99345864623196, -2.45138938198679*I, 0.942638555913623, 1.32703057887968*I, 7.33813274078958]

    Conversion from basic Python types::

        sage: pari(int(-5))
        -5
        sage: pari(long(2**150))
        1427247692705959881058285969449495136382746624
        sage: pari(float(pi))
        3.14159265358979
        sage: pari(complex(exp(pi*I/4)))
        0.707106781186548 + 0.707106781186547*I
        sage: pari(False)
        0
        sage: pari(True)
        1

    Some commands are just executed without returning a value::

        sage: pari("dummy = 0; kill(dummy)")
        sage: type(pari("dummy = 0; kill(dummy)"))
        <type 'NoneType'>
    """
    from sage.libs.pari.all import pari
    return pari(x)

def python(z, locals=None):
    """
    Return the closest Python/Sage equivalent of the given pari object.

    INPUT:

        - `z` -- pari object

        - `locals` -- optional dictionary used in fallback cases that
          involve sage_eval

    The component parts of a t_COMPLEX may be t_INT, t_REAL, t_INTMOD,
    t_FRAC, t_PADIC.  The components need not have the same type
    (e.g. if z=2+1.2*I then z.real() is t_INT while z.imag() is
    t_REAL().  They are converted as follows:

    t_INT:    ZZ[i]
    t_FRAC:   QQ(i)
    t_REAL:   ComplexField(prec) for equivalent precision
    t_INTMOD, t_PADIC: raise NotImplementedError

    EXAMPLES::

        sage: a = pari('(3+I)').python(); a
        i + 3
        sage: a.parent()
        Maximal Order in Number Field in i with defining polynomial x^2 + 1

        sage: a = pari('2^31-1').python(); a
        2147483647
        sage: a.parent()
        Integer Ring

        sage: a = pari('12/34').python(); a
        6/17
        sage: a.parent()
        Rational Field

        sage: a = pari('1.234').python(); a
        1.23400000000000000
        sage: a.parent()
        Real Field with 64 bits of precision

        sage: a = pari('(3+I)/2').python(); a
        1/2*i + 3/2
        sage: a.parent()
        Number Field in i with defining polynomial x^2 + 1

    Conversion of complex numbers: the next example is converting from
    an element of the Symbolic Ring, which goes via the string
    representation::

        sage: I = SR(I)
        sage: a = pari(1.0+2.0*I).python(); a
        1.00000000000000000 + 2.00000000000000000*I
        sage: type(a)
        <type 'sage.rings.complex_number.ComplexNumber'>
        sage: a.parent()
        Complex Field with 64 bits of precision

    For architecture-independent complex numbers, start from a
    suitable ComplexField::

        sage: z = pari(CC(1.0+2.0*I)); z
        1.00000000000000 + 2.00000000000000*I
        sage: a=z.python(); a
        1.00000000000000000 + 2.00000000000000000*I
        sage: a.parent()
        Complex Field with 64 bits of precision

    Vectors and matrices::

        sage: a = pari('[1,2,3,4]')
        sage: a
        [1, 2, 3, 4]
        sage: a.type()
        't_VEC'
        sage: b = a.python(); b
        [1, 2, 3, 4]
        sage: type(b)
        <type 'list'>

        sage: a = pari('[1,2;3,4]')
        sage: a.type()
        't_MAT'
        sage: b = a.python(); b
        [1 2]
        [3 4]
        sage: b.parent()
        Full MatrixSpace of 2 by 2 dense matrices over Integer Ring

        sage: a = pari('Vecsmall([1,2,3,4])')
        sage: a.type()
        't_VECSMALL'
        sage: a.python()
        [1, 2, 3, 4]

    We use the locals dictionary::

        sage: f = pari('(2/3)*x^3 + x - 5/7 + y')
        sage: x,y=var('x,y')
        sage: import sage.libs.pari.gen_py
        sage: sage.libs.pari.gen_py.python(f, {'x':x, 'y':y})
        2/3*x^3 + x + y - 5/7
        sage: sage.libs.pari.gen_py.python(f)
        Traceback (most recent call last):
        ...
        NameError: name 'x' is not defined

    Conversion of p-adics::

        sage: K = Qp(11,5)
        sage: x = K(11^-10 + 5*11^-7 + 11^-6); x
        11^-10 + 5*11^-7 + 11^-6 + O(11^-5)
        sage: y = pari(x); y
        11^-10 + 5*11^-7 + 11^-6 + O(11^-5)
        sage: y.sage()
        11^-10 + 5*11^-7 + 11^-6 + O(11^-5)
        sage: pari(K(11^-5)).sage()
        11^-5 + O(11^0)
    """
    from sage.libs.pari.pari_instance import prec_words_to_bits

    t = z.type()
    if t == "t_REAL":
        return RealField(prec_words_to_bits(z.precision()))(z)
    elif t == "t_FRAC":
         Q = RationalField()
         return Q(z)
    elif t == "t_INT":
         Z = IntegerRing()
         return Z(z)
    elif t == "t_COMPLEX":
        tx = z.real().type()
        ty = z.imag().type()
        if tx in ["t_INTMOD", "t_PADIC"] or ty in ["t_INTMOD", "t_PADIC"]:
            raise NotImplementedError, "No conversion to python available for t_COMPLEX with t_INTMOD or t_PADIC components"
        if tx == "t_REAL" or ty == "t_REAL":
            xprec = z.real().precision() # will be 0 if exact
            yprec = z.imag().precision() # will be 0 if exact
            if xprec == 0:
                prec = prec_words_to_bits(yprec)
            elif yprec == 0:
                prec = prec_words_to_bits(xprec)
            else:
                prec = max(prec_words_to_bits(xprec), prec_words_to_bits(yprec))
            R = RealField(prec)
            C = ComplexField(prec)
            return C(R(z.real()), R(z.imag()))
        if tx == "t_FRAC" or ty == "t_FRAC":
            return QuadraticField(-1,'i')([python(c) for c in list(z)])
        if tx == "t_INT" or ty == "t_INT":
            return QuadraticField(-1,'i').ring_of_integers()([python(c) for c in list(z)])
        raise NotImplementedError, "No conversion to python available for t_COMPLEX with components %s"%(tx,ty)
    elif t == "t_VEC":
        return [python(x) for x in z.python_list()]
    elif t == "t_VECSMALL":
        return z.python_list_small()
    elif t == "t_MAT":
        from sage.matrix.constructor import matrix
        return matrix(z.nrows(),z.ncols(),[python(z[i,j]) for i in range(z.nrows()) for j in range(z.ncols())])
    elif t == "t_PADIC":
        from sage.rings.padics.factory import Qp
        Z = IntegerRing()
        p = z.padicprime()
        rprec = Z(z.padicprec(p)) - Z(z._valp())
        K = Qp(Z(p), rprec)
        return K(z.lift())
    else:
        from sage.misc.sage_eval import sage_eval
        return sage_eval(str(z), locals=locals)
