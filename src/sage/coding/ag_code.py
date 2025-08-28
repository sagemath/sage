# sage.doctest: needs sage.rings.finite_rings sage.schemes
"""
AG codes

Algebraic geometry codes or shortly AG codes are linear codes defined using
functions or differentials on algebraic curves over finite fields. Sage
implements evaluation AG codes and differential AG codes as Goppa defined in
[Gop1981]_ and provides decoding algorithms for them in full generality.

EXAMPLES::

    sage: k.<a> = GF(4)
    sage: A.<x,y> = AffineSpace(k, 2)
    sage: C = Curve(y^2 + y - x^3)
    sage: F = C.function_field()
    sage: pls = F.places()
    sage: p = C([0,0])
    sage: Q, = p.places()
    sage: pls.remove(Q)
    sage: G = 5*Q
    sage: codes.EvaluationAGCode(pls, G)
    [8, 5] evaluation AG code over GF(4)
    sage: codes.DifferentialAGCode(pls, G)
    [8, 3] differential AG code over GF(4)

As is well known, the two kinds of AG codes are dual to each other. ::

    sage: E = codes.EvaluationAGCode(pls, G)
    sage: D = codes.DifferentialAGCode(pls, G)
    sage: E.dual_code() == D
    True
    sage: D.dual_code() == E
    True

Decoders for both evaluation and differential AG codes are available.

.. toctree::

  ag_code_decoders

A natural generalization of classical Goppa codes is Cartier codes [Cou2014]_. Cartier codes are
subfield subcodes of differential AG codes.

EXAMPLES::

    sage: F.<a> = GF(9)
    sage: P.<x,y,z> = ProjectiveSpace(F, 2);
    sage: C = Curve(x^3*y + y^3*z + x*z^3)
    sage: F = C.function_field()
    sage: pls = F.places()
    sage: Z, = C([0,0,1]).places()
    sage: pls.remove(Z)
    sage: G = 3*Z
    sage: codes.CartierCode(pls, G)  # long time
    [9, 4] Cartier code over GF(3)

AUTHORS:

- Kwankyu Lee (2019-03): initial version
"""

# ****************************************************************************
#       Copyright (C) 2019 Kwankyu Lee <kwankyu@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.function_field.place import FunctionFieldPlace

from .linear_code import (AbstractLinearCode,
                          LinearCodeGeneratorMatrixEncoder,
                          LinearCodeSyndromeDecoder)

from .ag_code_decoders import (EvaluationAGCodeUniqueDecoder,
                               EvaluationAGCodeEncoder,
                               DifferentialAGCodeUniqueDecoder,
                               DifferentialAGCodeEncoder)


class AGCode(AbstractLinearCode):
    """
    Base class of algebraic geometry codes.

    A subclass of this class is required to define ``_function_field``
    attribute that refers to an abstract functiom field or the function field
    of the underlying curve used to construct a code of the class.
    """
    def base_function_field(self):
        """
        Return the function field used to construct the code.

        EXAMPLES::

            sage: k.<a> = GF(4)
            sage: A.<x,y> = AffineSpace(k, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: p = C([0,0])
            sage: Q, = p.places()
            sage: pls.remove(Q)
            sage: G = 5*Q
            sage: code = codes.EvaluationAGCode(pls, G)
            sage: code.base_function_field()
            Function field in y defined by y^2 + y + x^3
        """
        return self._function_field


class EvaluationAGCode(AGCode):
    """
    Evaluation AG code defined by rational places ``pls`` and a divisor ``G``.

    INPUT:

    - ``pls`` -- list of rational places of a function field

    - ``G`` -- a divisor whose support is disjoint from ``pls``

    If ``G`` is a place, then it is regarded as a prime divisor.

    EXAMPLES::

        sage: k.<a> = GF(4)
        sage: A.<x,y> = AffineSpace(k, 2)
        sage: C = Curve(y^2 + y - x^3)
        sage: F = C.function_field()
        sage: pls = F.places()
        sage: Q, = C.places_at_infinity()
        sage: pls.remove(Q)
        sage: G = 5*Q
        sage: codes.EvaluationAGCode(pls, G)
        [8, 5] evaluation AG code over GF(4)

        sage: G = F.get_place(5)
        sage: codes.EvaluationAGCode(pls, G)
        [8, 5] evaluation AG code over GF(4)
    """
    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, pls, G):
        """
        Initialize.

        TESTS::

            sage: k.<a> = GF(4)
            sage: A.<x,y> = AffineSpace(k, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: Q, = C.places_at_infinity()
            sage: pls = [pl for pl in pls if pl != Q]
            sage: G = 5*Q
            sage: code = codes.EvaluationAGCode(pls, G)
            sage: TestSuite(code).run()
        """
        if issubclass(type(G), FunctionFieldPlace):
            G = G.divisor()  # place is converted to a prime divisor

        F = G.parent().function_field()
        K = F.constant_base_field()
        n = len(pls)

        if any(p.degree() > 1 for p in pls):
            raise ValueError("there is a nonrational place among the places")

        if any(p in pls for p in G.support()):
            raise ValueError("the support of the divisor is not disjoint from the places")

        self._registered_encoders['evaluation'] = EvaluationAGCodeEncoder
        self._registered_decoders['K'] = EvaluationAGCodeUniqueDecoder

        super().__init__(K, n, default_encoder_name='evaluation',
                               default_decoder_name='K')

        # compute basis functions associated with a generator matrix
        basis_functions = G.basis_function_space()
        m = matrix([vector(K, [b.evaluate(p) for p in pls]) for b in basis_functions])
        I = MatrixSpace(K, m.nrows()).identity_matrix()
        mI = m.augment(I)
        mI.echelonize()
        M = mI.submatrix(0, 0, m.nrows(), m.ncols())
        T = mI.submatrix(0, m.ncols())
        r = M.rank()

        self._generator_matrix = M.submatrix(0, 0, r)
        self._basis_functions = [sum(c * b for c, b in zip(T[i], basis_functions))
                                     for i in range(r)]

        self._pls = tuple(pls)
        self._G = G

        self._function_field = F

    def __eq__(self, other):
        """
        Test equality of ``self`` with ``other``.

        EXAMPLES::

            sage: k.<a> = GF(4)
            sage: A.<x,y> = AffineSpace(k, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: Q, = C.places_at_infinity()
            sage: pls.remove(Q)
            sage: codes.EvaluationAGCode(pls, 5*Q) == codes.EvaluationAGCode(pls, 6*Q)
            False
        """
        if self is other:
            return True

        if not isinstance(other, EvaluationAGCode):
            return False

        return self._pls == other._pls and self._G == other._G

    def __hash__(self):
        """
        Return the hash value of ``self``.

        EXAMPLES::

            sage: k.<a> = GF(4)
            sage: A.<x,y> = AffineSpace(k, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: Q, = C.places_at_infinity()
            sage: pls.remove(Q)
            sage: code = codes.EvaluationAGCode(pls, 5*Q)
            sage: {code: 1}
            {[8, 5] evaluation AG code over GF(4): 1}
        """
        return hash((self._pls, self._G))

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: k.<a> = GF(4)
            sage: A.<x,y> = AffineSpace(k, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: Q, = C.places_at_infinity()
            sage: pls.remove(Q)
            sage: codes.EvaluationAGCode(pls, 7*Q)
            [8, 7] evaluation AG code over GF(4)
        """
        return "[{}, {}] evaluation AG code over GF({})".format(
                self.length(), self.dimension(), self.base_field().cardinality())

    def _latex_(self):
        r"""
        Return the latex representation of ``self``.

        EXAMPLES::

            sage: k.<a> = GF(4)
            sage: A.<x,y> = AffineSpace(k, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: Q, = C.places_at_infinity()
            sage: pls.remove(Q)
            sage: code = codes.EvaluationAGCode(pls, 3*Q)
            sage: latex(code)
            [8, 3]\text{ evaluation AG code over }\Bold{F}_{2^{2}}
        """
        return r"[{}, {}]\text{{ evaluation AG code over }}{}".format(
                self.length(), self.dimension(), self.base_field()._latex_())

    def basis_functions(self):
        r"""
        Return the basis functions associated with the generator matrix.

        EXAMPLES::

            sage: k.<a> = GF(4)
            sage: A.<x,y> = AffineSpace(k, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: Q, = C.places_at_infinity()
            sage: pls.remove(Q)
            sage: code = codes.EvaluationAGCode(pls, 3*Q)
            sage: code.basis_functions()
            (y + a*x + 1, y + x, (a + 1)*x)
            sage: matrix([[f.evaluate(p) for p in pls] for f in code.basis_functions()])
            [    1     0     0     1     a a + 1     1     0]
            [    0     1     0     1     1     0 a + 1     a]
            [    0     0     1     1     a     a a + 1 a + 1]
        """
        return tuple(self._basis_functions)

    def generator_matrix(self):
        r"""
        Return a generator matrix of the code.

        EXAMPLES::

            sage: k.<a> = GF(4)
            sage: A.<x,y> = AffineSpace(k, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: Q, = C.places_at_infinity()
            sage: pls.remove(Q)
            sage: code = codes.EvaluationAGCode(pls, 3*Q)
            sage: code.generator_matrix()
            [    1     0     0     1     a a + 1     1     0]
            [    0     1     0     1     1     0 a + 1     a]
            [    0     0     1     1     a     a a + 1 a + 1]
        """
        return self._generator_matrix

    def designed_distance(self):
        """
        Return the designed distance of the AG code.

        If the code is of dimension zero, then a :exc:`ValueError` is raised.

        EXAMPLES::

            sage: k.<a> = GF(4)
            sage: A.<x,y> = AffineSpace(k, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: Q, = C.places_at_infinity()
            sage: pls.remove(Q)
            sage: code = codes.EvaluationAGCode(pls, 3*Q)
            sage: code.designed_distance()
            5
        """
        if self.dimension() == 0:
            raise ValueError("not defined for zero code")

        d = self.length() - self._G.degree()
        return d if d > 0 else 1


class DifferentialAGCode(AGCode):
    """
    Differential AG code defined by rational places ``pls`` and a divisor ``G``.

    INPUT:

    - ``pls`` -- list of rational places of a function field

    - ``G`` -- a divisor whose support is disjoint from ``pls``

    If ``G`` is a place, then it is regarded as a prime divisor.

    EXAMPLES::

        sage: k.<a> = GF(4)
        sage: A.<x,y> = AffineSpace(k, 2)
        sage: C = A.curve(y^3 + y - x^4)
        sage: Q = C.places_at_infinity()[0]
        sage: O = C([0,0]).place()
        sage: pls = [p for p in C.places() if p not in [O, Q]]
        sage: G = -O + 3*Q
        sage: codes.DifferentialAGCode(pls, -O + Q)
        [3, 2] differential AG code over GF(4)

        sage: F = C.function_field()
        sage: G = F.get_place(1)
        sage: codes.DifferentialAGCode(pls, G)
        [3, 1] differential AG code over GF(4)
    """
    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, pls, G):
        """
        Initialize.

        TESTS::

            sage: k.<a> = GF(4)
            sage: A.<x,y> = AffineSpace(k, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: Q, = C.places_at_infinity()
            sage: pls.remove(Q)
            sage: code = codes.DifferentialAGCode(pls, 3*Q)
            sage: TestSuite(code).run()
        """
        if issubclass(type(G), FunctionFieldPlace):
            G = G.divisor()  # place is converted to a prime divisor

        F = G.parent().function_field()
        K = F.constant_base_field()
        n = len(pls)

        if any(p.degree() > 1 for p in pls):
            raise ValueError("there is a nonrational place among the places")

        if any(p in pls for p in G.support()):
            raise ValueError("the support of the divisor is not disjoint from the places")

        self._registered_encoders['residue'] = DifferentialAGCodeEncoder
        self._registered_decoders['K'] = DifferentialAGCodeUniqueDecoder

        super().__init__(K, n, default_encoder_name='residue',
                               default_decoder_name='K')

        # compute basis differentials associated with a generator matrix
        basis_differentials = (-sum(pls) + G).basis_differential_space()
        m = matrix([vector(K, [w.residue(p) for p in pls]) for w in basis_differentials])
        I = MatrixSpace(K, m.nrows()).identity_matrix()
        mI = m.augment(I)
        mI.echelonize()
        M = mI.submatrix(0, 0, m.nrows(), m.ncols())
        T = mI.submatrix(0, m.ncols())
        r = M.rank()

        self._generator_matrix = M.submatrix(0, 0, r)
        self._basis_differentials = [sum(c * w for c, w in zip(T[i], basis_differentials))
                                     for i in range(r)]

        self._pls = tuple(pls)
        self._G = G
        self._function_field = F

    def __eq__(self, other):
        """
        Test equality of ``self`` with ``other``.

        EXAMPLES::

            sage: k.<a> = GF(4)
            sage: A.<x,y> = AffineSpace(k, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: Q, = C.places_at_infinity()
            sage: pls.remove(Q)
            sage: c1 = codes.DifferentialAGCode(pls, 3*Q)
            sage: c2 = codes.DifferentialAGCode(pls, 3*Q)
            sage: c1 is c2
            False
            sage: c1 == c2
            True
        """
        if self is other:
            return True

        if not isinstance(other, DifferentialAGCode):
            return False

        return self._pls == other._pls and self._G == other._G

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: k.<a> = GF(4)
            sage: A.<x,y> = AffineSpace(k, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: Q, = C.places_at_infinity()
            sage: pls.remove(Q)
            sage: code = codes.DifferentialAGCode(pls, 3*Q)
            sage: {code: 1}
            {[8, 5] differential AG code over GF(4): 1}
        """
        return hash((self._pls, self._G))

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: k.<a> = GF(4)
            sage: A.<x,y> = AffineSpace(k, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: Q, = C.places_at_infinity()
            sage: pls.remove(Q)
            sage: codes.DifferentialAGCode(pls, 3*Q)
            [8, 5] differential AG code over GF(4)
        """
        return "[{}, {}] differential AG code over GF({})".format(
                self.length(), self.dimension(), self.base_field().cardinality())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: k.<a> = GF(4)
            sage: A.<x,y> = AffineSpace(k, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: Q, = C.places_at_infinity()
            sage: pls.remove(Q)
            sage: code = codes.DifferentialAGCode(pls, 3*Q)
            sage: latex(code)
            [8, 5]\text{ differential AG code over }\Bold{F}_{2^{2}}
        """
        return r"[{}, {}]\text{{ differential AG code over }}{}".format(
                self.length(), self.dimension(), self.base_field()._latex_())

    def basis_differentials(self):
        r"""
        Return the basis differentials associated with the generator matrix.

        EXAMPLES::

            sage: k.<a> = GF(4)
            sage: A.<x,y> = AffineSpace(k, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: Q, = C.places_at_infinity()
            sage: pls.remove(Q)
            sage: code = codes.DifferentialAGCode(pls, 3*Q)
            sage: matrix([[w.residue(p) for p in pls]
            ....:         for w in code.basis_differentials()])
            [    1     0     0     0     0 a + 1 a + 1     1]
            [    0     1     0     0     0 a + 1     a     0]
            [    0     0     1     0     0     a     1     a]
            [    0     0     0     1     0     a     0 a + 1]
            [    0     0     0     0     1     1     1     1]
        """
        return tuple(self._basis_differentials)

    def generator_matrix(self):
        """
        Return a generator matrix of the code.

        EXAMPLES::

            sage: k.<a> = GF(4)
            sage: A.<x,y> = AffineSpace(k, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: Q, = C.places_at_infinity()
            sage: pls.remove(Q)
            sage: code = codes.DifferentialAGCode(pls, 3*Q)
            sage: code.generator_matrix()
            [    1     0     0     0     0 a + 1 a + 1     1]
            [    0     1     0     0     0 a + 1     a     0]
            [    0     0     1     0     0     a     1     a]
            [    0     0     0     1     0     a     0 a + 1]
            [    0     0     0     0     1     1     1     1]
        """
        return self._generator_matrix

    def designed_distance(self):
        """
        Return the designed distance of the differential AG code.

        If the code is of dimension zero, then a :exc:`ValueError` is raised.

        EXAMPLES::

            sage: k.<a> = GF(4)
            sage: A.<x,y> = AffineSpace(k, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: Q, = C.places_at_infinity()
            sage: pls.remove(Q)
            sage: code = codes.DifferentialAGCode(pls, 3*Q)
            sage: code.designed_distance()
            3
        """
        if self.dimension() == 0:
            raise ValueError("not defined for zero code")

        d = self._G.degree() - 2 * self._function_field.genus() + 2
        return d if d > 0 else 1


class CartierCode(AGCode):
    r"""
    Cartier code defined by rational places ``pls`` and a divisor ``G`` of a function field.

    INPUT:

    - ``pls`` -- list of rational places

    - ``G`` -- a divisor whose support is disjoint from ``pls``

    - ``r`` -- integer (default: 1)

    - ``name`` -- string; name of the generator of the subfield `\GF{p^r}`

    OUTPUT: cartier code over `\GF{p^r}` where `p` is the characteristic of the
    base constant field of the function field

    Note that if ``r`` is 1 the default, then ``name`` can be omitted.

    EXAMPLES::

        sage: F.<a> = GF(9)
        sage: P.<x,y,z> = ProjectiveSpace(F, 2);
        sage: C = Curve(x^3*y + y^3*z + x*z^3)
        sage: F = C.function_field()
        sage: pls = F.places()
        sage: Z, = C(0,0,1).places()
        sage: pls.remove(Z)
        sage: G = 3*Z
        sage: code = codes.CartierCode(pls, G)  # long time
        sage: code.minimum_distance()           # long time
        2
    """
    def __init__(self, pls, G, r=1, name=None):
        """
        Initialize.

        TESTS::

            sage: F.<a> = GF(9)
            sage: P.<x,y,z> = ProjectiveSpace(F, 2);
            sage: C = Curve(x^3*y + y^3*z + x*z^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: Z, = C(0,0,1).places()
            sage: pls.remove(Z)
            sage: G = 3*Z
            sage: code = codes.CartierCode(pls, G)  # long time
            sage: TestSuite(code).run()             # long time
        """
        F = G.parent().function_field()
        K = F.constant_base_field()

        if any(p.degree() > 1 for p in pls):
            raise ValueError("there is a nonrational place among the places")

        if any(p in pls for p in G.support()):
            raise ValueError("the support of the divisor is not disjoint from the places")

        if K.degree() % r != 0:
            raise ValueError("{} does not divide the degree of the constant base field".format(r))

        n = len(pls)
        D = sum(pls)
        p = K.characteristic()

        subfield = K.subfield(r, name=name)

        # compute a basis R of the space of differentials in Omega(G - D)
        # fixed by the Cartier operator
        E = G - D

        Grp = E.parent()  # group of divisors
        V, fr_V, to_V = E.differential_space()

        EE = Grp(0)
        dic = E.dict()
        for place in dic:
            mul = dic[place]
            if mul > 0:
                mul = mul // p**r
            EE += mul * place

        W, fr_W, to_W = EE.differential_space()

        a = K.gen()
        field_basis = [a**i for i in range(K.degree())] # over prime subfield
        basis = E.basis_differential_space()

        m = []
        for w in basis:
            for c in field_basis:
                cw = F(c) * w # c does not coerce...
                carcw = cw
                for i in range(r): # apply cartier r times
                    carcw = carcw.cartier()
                m.append([f for e in to_W(carcw - cw) for f in vector(e)])

        ker = matrix(m).kernel()

        R = []
        s = len(field_basis)
        ncols = s * len(basis)
        for row in ker.basis():
            v = vector([K(row[d:d+s]) for d in range(0,ncols,s)])
            R.append(fr_V(v))

        # construct a generator matrix
        m = []
        col_index = D.support()
        for w in R:
            row = []
            for p in col_index:
                res = w.residue(p).trace() # lies in constant base field
                c = subfield(res) # as w is Cartier fixed
                row.append(c)
            m.append(row)

        self._generator_matrix = matrix(m).row_space().basis_matrix()

        self._pls = tuple(pls)
        self._G = G
        self._r = r
        self._function_field = F

        self._registered_encoders['GeneratorMatrix'] = LinearCodeGeneratorMatrixEncoder
        self._registered_decoders['Syndrome'] = LinearCodeSyndromeDecoder

        super().__init__(subfield, n,
                         default_encoder_name='GeneratorMatrix',
                         default_decoder_name='Syndrome')

    def __eq__(self, other):
        """
        Test equality of ``self`` with ``other``.

        EXAMPLES::

            sage: F.<a> = GF(9)
            sage: P.<x,y,z> = ProjectiveSpace(F, 2);
            sage: C = Curve(x^3*y + y^3*z + x*z^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: Z, = C(0,0,1).places()
            sage: pls.remove(Z)
            sage: c1 = codes.CartierCode(pls, 3*Z)  # long time
            sage: c2 = codes.CartierCode(pls, 1*Z)  # long time
            sage: c1 == c2                          # long time
            False
        """
        if self is other:
            return True

        if not isinstance(other, CartierCode):
            return False

        return self._pls == other._pls and self._G == other._G and self._r == other._r

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: F.<a> = GF(9)
            sage: P.<x,y,z> = ProjectiveSpace(F, 2);
            sage: C = Curve(x^3*y + y^3*z + x*z^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: Z, = C(0,0,1).places()
            sage: pls.remove(Z)
            sage: G = 3*Z
            sage: code = codes.CartierCode(pls, G)  # long time
            sage: {code: 1}                         # long time
            {[9, 4] Cartier code over GF(3): 1}
        """
        return hash((self._pls, self._G, self._r))

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: F.<a> = GF(9)
            sage: P.<x,y,z> = ProjectiveSpace(F, 2);
            sage: C = Curve(x^3*y + y^3*z + x*z^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: Z, = C(0,0,1).places()
            sage: pls.remove(Z)
            sage: G = 3*Z
            sage: codes.CartierCode(pls, G)  # long time
            [9, 4] Cartier code over GF(3)
        """
        return "[{}, {}] Cartier code over GF({})".format(
                self.length(), self.dimension(), self.base_field().cardinality())

    def _latex_(self):
        r"""
        Return the latex representation of ``self``.

        EXAMPLES::

            sage: F.<a> = GF(9)
            sage: P.<x,y,z> = ProjectiveSpace(F, 2);
            sage: C = Curve(x^3*y + y^3*z + x*z^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: Z, = C(0,0,1).places()
            sage: pls.remove(Z)
            sage: G = 3*Z
            sage: code = codes.CartierCode(pls, G)  # long time
            sage: latex(code)                       # long time
            [9, 4]\text{ Cartier code over }\Bold{F}_{3}
        """
        return r"[{}, {}]\text{{ Cartier code over }}{}".format(
                self.length(), self.dimension(), self.base_field()._latex_())

    def generator_matrix(self):
        r"""
        Return a generator matrix of the Cartier code.

        EXAMPLES::

            sage: F.<a> = GF(9)
            sage: P.<x,y,z> = ProjectiveSpace(F, 2);
            sage: C = Curve(x^3*y + y^3*z + x*z^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: Z, = C(0,0,1).places()
            sage: pls.remove(Z)
            sage: G = 3*Z
            sage: code = codes.CartierCode(pls, G)  # long time
            sage: code.generator_matrix()           # long time
            [1 0 0 2 2 0 2 2 0]
            [0 1 0 2 2 0 2 2 0]
            [0 0 1 0 0 0 0 0 2]
            [0 0 0 0 0 1 0 0 2]
        """
        return self._generator_matrix

    def designed_distance(self):
        """
        Return the designed distance of the Cartier code.

        The designed distance is that of the differential code of which the
        Cartier code is a subcode.

        EXAMPLES::

            sage: F.<a> = GF(9)
            sage: P.<x,y,z> = ProjectiveSpace(F, 2);
            sage: C = Curve(x^3*y + y^3*z + x*z^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: Z, = C(0,0,1).places()
            sage: pls.remove(Z)
            sage: G = 3*Z
            sage: code = codes.CartierCode(pls, G)  # long time
            sage: code.designed_distance()          # long time
            1
        """
        if self.dimension() == 0:
            raise ValueError("not defined for zero code")

        d = self._G.degree() - 2 * self._function_field.genus() + 2
        return d if d > 0 else 1
