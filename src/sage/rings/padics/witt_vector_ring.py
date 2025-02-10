"""
Witt vector rings: implementation
"""
from itertools import product

from sage.categories.commutative_rings import CommutativeRings
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.integer_ring import ZZ
from sage.rings.padics.factory import Zp, Zq
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.ring import CommutativeRing
from sage.sets.primes import Primes
from sage.structure.unique_representation import UniqueRepresentation

from .witt_vector import WittVector_base, WittVector_non_p_typical, WittVector_p_typical

_Primes = Primes()


def _fast_char_p_power(x, n, p=None):
    r"""
    Raise x^n power in characteristic p.

    If x is not an element of a ring of characteristic p, this throws an error.

    If x is an element of GF(p^k), this is already fast.
    However, is x is a polynomial, this seems to be slow?

    EXAMPLES::

        sage: from sage.rings.padics.witt_vector_ring import _fast_char_p_power as fpow
        sage: t = GF(1913)(33)
        sage: fpow(t, 77)
        1371
    """
    if n not in ZZ:
        raise ValueError(f'Exponent {n} is not an integer')
    if n == 0 or x == 1:
        return x.parent().one()
    if x.parent().characteristic() not in _Primes:
        raise ValueError(f'{x} is not in a ring of prime characteristic')

    x_is_Polynomial = isinstance(x, Polynomial)
    x_is_MPolynomial = isinstance(x, MPolynomial)

    if not (x_is_Polynomial or x_is_MPolynomial):
        return x**n
    if x.is_gen():
        return x**n
    if n < 0:
        x = x**-1  # This may throw an error.
        n = -n

    P = x.parent()
    if p is None:
        p = P.characteristic()
    base_p_digits = ZZ(n).digits(base=p)

    x_to_the_n = 1

    for p_exp, digit in enumerate(base_p_digits):
        if digit == 0:
            continue
        inner_term = x**digit
        term_dict = {}
        for e_int_or_tuple, c in inner_term.dict().items():
            power = p**p_exp
            new_c = _fast_char_p_power(c, power)
            new_e_tuple = None
            if x_is_Polynomial:  # Then the dict keys are ints
                new_e_tuple = e_int_or_tuple * power
            elif x_is_MPolynomial:  # Then the dict keys are ETuples
                new_e_tuple = e_int_or_tuple.emul(power)
            term_dict[new_e_tuple] = new_c
        term = P(term_dict)
        x_to_the_n *= term

    return x_to_the_n


_fcppow = _fast_char_p_power


class WittVectorRing_base(CommutativeRing, UniqueRepresentation):

    Element = WittVector_base

    def __init__(self, base_ring, prec, prime, algorithm='none', category=None):
        self.prec = prec
        self.prime = prime

        self._algorithm = algorithm
        self.sum_polynomials = None
        self.prod_polynomials = None

        if algorithm == 'standard':
            self._generate_sum_and_product_polynomials(base_ring)

        if category is None:
            category = CommutativeRings()
        CommutativeRing.__init__(self, base_ring, category=category)

    def __iter__(self):
        for t in product(self.base(), repeat=self.prec):
            yield self(t)

    def _repr_(self):
        return f"Ring of Witt Vectors of length {self.prec} over {self.base()}"

    def _coerce_map_from_(self, S):
        # Question: do we return True is S == self.base()?
        # We have the teichmuller lift, but I don't think that's
        # a "coercion" map, per se.
        return (S is ZZ)

    def _element_constructor_(self, x):
        if x in ZZ:
            return self.element_class(self, self._int_to_vector(x))
        elif isinstance(x, tuple) or isinstance(x, list):
            return self.element_class(self, x)
        else:
            return NotImplemented

    def _int_to_vector(self, k):
        p = self.prime

        should_negate = False
        if k < 0:
            k = -k
            should_negate = True

        vec_k = [k]
        for n in range(1, self.prec):
            total = k - k**(p**n) - sum(p**(n-i) * vec_k[n-i]**(p**i) for i in range(1, n))
            total //= p**n
            vec_k.append(total)

        if should_negate:
            if p == 2:
                return NotImplemented
            else:
                vec_k = [-x for x in vec_k]

        return vec_k

    def _generate_sum_and_product_polynomials(self, base):
        p = self.prime
        prec = self.prec
        x_var_names = [f'X{i}' for i in range(prec)]
        y_var_names = [f'Y{i}' for i in range(prec)]
        var_names = x_var_names + y_var_names

        # Okay, what's going on here? Sage, by default, relies on
        # Singular for Multivariate Polynomial Rings, but Singular uses
        # only SIXTEEN bits (unsigned) to store its exponents. So if we
        # want exponents larger than 2^16 - 1, we have to use the
        # generic implementation. However, after some experimentation,
        # it seems like the generic implementation is faster?
        #
        # After trying to compute S_4 for p=5, it looks like generic is
        # faster for  very small polys, and MUCH slower for large polys.
        # So we'll default to singular unless we can't use it.
        #
        # Remark: Since when is SIXTEEN bits sufficient for anyone???
        #
        if p**(prec-1) >= 2**16:
            implementation = 'generic'
        else:
            implementation = 'singular'

        # We first generate the "universal" polynomials and then project
        # to the base ring.
        R = PolynomialRing(ZZ, var_names, implementation=implementation)
        x_y_vars = R.gens()
        x_vars = x_y_vars[:prec]
        y_vars = x_y_vars[prec:]

        self.sum_polynomials = [0]*(self.prec)
        for n in range(prec):
            s_n = x_vars[n] + y_vars[n]
            for i in range(n):
                s_n += (x_vars[i]**(p**(n-i)) + y_vars[i]**(p**(n-i)) - self.sum_polynomials[i]**(p**(n-i))) / p**(n-i)
            self.sum_polynomials[n] = R(s_n)

        self.prod_polynomials = [x_vars[0] * y_vars[0]] + [0]*(self.prec)
        for n in range(1, prec):
            x_poly = sum([p**i * x_vars[i]**(p**(n-i)) for i in range(n+1)])
            y_poly = sum([p**i * y_vars[i]**(p**(n-i)) for i in range(n+1)])
            p_poly = sum([p**i * self.prod_polynomials[i]**(p**(n-i)) for i in range(n)])
            p_n = (x_poly*y_poly - p_poly) // p**n
            self.prod_polynomials[n] = p_n

        # We have to use generic here, because Singular doesn't support
        # Polynomial Rings over Polynomial Rings. For example,
        # ``PolynomialRing(GF(5)['x'], ['X', 'Y'], implementation='singular')``
        # will fail.
        S = PolynomialRing(base, x_y_vars, implementation='generic')
        for n in range(prec):
            self.sum_polynomials[n] = S(self.sum_polynomials[n])
            self.prod_polynomials[n] = S(self.prod_polynomials[n])

    def characteristic(self):
        p = self.prime
        if self.base()(p).is_unit():
            # If p is invertible, W_n(R) is isomorphic to R^n.
            return self.base().characteristic()

        # This is a conjecture. It's known for char(R) == p.
        return p**(self.prec-1) * self.base().characteristic()

    def precision(self):
        return self.prec

    def random_element(self):
        return self.element_class(self, tuple(self.base().random_element()
                                              for _ in range(self.prec)))

    def teichmuller_lift(self, x):
        if x not in self.base():
            raise Exception(f'{x} not in {self.base()}')
        return self.element_class(self, (x,) + tuple(0 for _ in range(self.prec-1)))

    def is_finite(self):
        return self.base().is_finite()

    def cardinality(self):
        return self.base().cardinality()**(self.prec)


class WittVectorRing_p_typical(WittVectorRing_base):

    Element = WittVector_p_typical

    def __init__(self, base_ring, prec, prime, algorithm=None, category=None):
        WittVectorRing_base.__init__(self, base_ring, prec, prime,
                               algorithm=algorithm, category=category)

        if algorithm == 'finotti':
            self.generate_binomial_table()

    def _int_to_vector(self, k):
        p = self.prime
        R = Zp(p, prec=self.prec+1, type='fixed-mod')
        F = GF(p)
        B = self.base()

        series = R(k)
        witt_vector = []
        for _ in range(self.prec):
            # Probably slightly faster to do "series % p," but this way, temp is in F_p
            temp = F(series)
            witt_vector.append(B(temp))  # make sure elements of vector are in base
            series = (series - R.teichmuller(temp)) // p
        return witt_vector

    def generate_binomial_table(self):
        import numpy as np
        p = self.prime
        R = Zp(p, prec=self.prec+1, type='fixed-mod')
        v_p = ZZ.valuation(p)
        table = [[0]]
        for k in range(1, self.prec+1):
            row = np.empty(p**k, dtype=int)
            row[0] = 0
            prev_bin = 1
            for i in range(1, p**k // 2 + 1):
                val = v_p(i)
                # Instead of calling binomial each time, we compute the coefficients
                # recursively. This is MUCH faster.
                next_bin = prev_bin * (p**k - (i-1)) // i
                prev_bin = next_bin
                series = R(-next_bin // p**(k-val))
                for _ in range(val):
                    temp = series % p
                    series = (series - R.teichmuller(temp)) // p
                row[i] = ZZ(series % p)
                row[p**k - i] = row[i]  # binomial coefficients are symmetric
            table.append(row)
        self.binomial_table = table

    def eta_bar(self, vec, eta_index):
        vec = tuple(x for x in vec if x != 0)  # strip zeroes

        # special cases
        if len(vec) <= 1:
            return 0
        if eta_index == 0:
            return sum(vec)

        # renaming to match notation in paper
        k = eta_index
        p = self.prime
        # if vec = (x,y), we know what to do: Theorem 8.6
        if len(vec) == 2:
            # Here we have to check if we've pre-computed already
            x, y = vec
            scriptN = [[None] for _ in range(k+1)]  # each list starts with None, so that indexing matches paper
            # calculate first N_t scriptN's
            for t in range(1, k+1):
                for i in range(1, p**t):
                    scriptN[t].append(self.binomial_table[t][i] * _fcppow(x, i) * _fcppow(y, p**t - i))
            indexN = [p**i - 1 for i in range(k+1)]
            for t in range(2, k+1):
                for l in range(1, t):
                    # append scriptN_{t, N_t+l}
                    next_scriptN = self.eta_bar(scriptN[t-l][1:indexN[t-l]+t-l], l)
                    scriptN[t].append(next_scriptN)
            return sum(scriptN[k][1:])

        # if vec is longer, we split and recurse: Proposition 5.4
        # This is where we need to using multiprocessing.
        else:
            m = len(vec) // 2
            v_1 = vec[:m]
            v_2 = vec[m:]
            s_1 = sum(v_1)
            s_2 = sum(v_2)
            scriptM = [[] for _ in range(k+1)]
            for t in range(1, k+1):
                scriptM[t].append(self.eta_bar(v_1, t))
                scriptM[t].append(self.eta_bar(v_2, t))
                scriptM[t].append(self.eta_bar((s_1, s_2), t))
            for t in range(2, k+1):
                for s in range(1, t):
                    result = self.eta_bar(scriptM[t-s], s)
                    scriptM[t].append(result)
            return sum(scriptM[k])


class WittVectorRing_finite_field(WittVectorRing_p_typical):
    def __init__(self, base_ring, prec, prime, category=None):
        WittVectorRing_p_typical.__init__(self, base_ring, prec, prime,
                                    algorithm='Zq_isomorphism',
                                    category=category)

    def _series_to_vector(self, series):
        F = self.base()  # known to be finite
        R = Zq(F.cardinality(), prec=self.prec, type='fixed-mod', modulus=F.polynomial(), names=['z'])
        K = R.residue_field()
        p = self.prime

        series = R(series)
        witt_vector = []
        for i in range(self.prec):
            temp = K(series)
            elem = temp.polynomial()(F.gen())  # hack to convert to F (K != F for some reason)
            witt_vector.append(elem**(p**i))
            series = (series - R.teichmuller(temp)) // p
        return witt_vector

    def _vector_to_series(self, vec):
        F = self.base()
        R = Zq(F.cardinality(), prec=self.prec, type='fixed-mod', modulus=F.polynomial(), names=['z'])
        K = R.residue_field()
        p = self.prime

        series = R.zero()
        for i in range(self.prec):
            temp = vec[i].nth_root(p**i)
            elem = temp.polynomial()(K.gen())
            # hack to convert to K (F != K for some reason)

            series += p**i * R.teichmuller(elem)
        return series


class WittVectorRing_non_p_typical(WittVectorRing_base):

    Element = WittVector_non_p_typical

    def __init__(self, base_ring, prec, prime, algorithm=None, category=None):
        WittVectorRing_base.__init__(self, base_ring, prec, prime,
                               algorithm=algorithm, category=category)

    def _repr_(self):
        return f"Ring of {self.prime}-Witt Vectors of length {self.prec} over {self.base()}"


class WittVectorRing_p_invertible(WittVectorRing_non_p_typical):
    def __init__(self, base_ring, prec, prime, category=None):
        WittVectorRing_non_p_typical.__init__(self, base_ring, prec, prime,
                                        algorithm='standard_otf',
                                        category=category)
