# sage.doctest:           optional - sage.modules           (because __init__ constructs a vector space)
# some tests are marked # optional - sage.libs.pari         (because they use finite fields)
r"""
Maximal orders of global function fields
"""

from sage.misc.cachefunc import cached_method

from .ideal import FunctionFieldIdeal_global
from .order_polymod import FunctionFieldMaximalOrder_polymod


class FunctionFieldMaximalOrder_global(FunctionFieldMaximalOrder_polymod):
    """
    Maximal orders of global function fields.

    INPUT:

    - ``field`` -- function field to which this maximal order belongs

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
        sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
        sage: L.maximal_order()
        Maximal order of Function field in y defined by y^4 + x*y + 4*x + 1
    """

    def __init__(self, field):
        """
        Initialize.

        TESTS::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.maximal_order()
            sage: TestSuite(O).run()
        """
        FunctionFieldMaximalOrder_polymod.__init__(self, field, ideal_class=FunctionFieldIdeal_global)

    @cached_method
    def p_radical(self, prime):
        """
        Return the ``prime``-radical of the maximal order.

        INPUT:

        - ``prime`` -- prime ideal of the maximal order of the base
          rational function field

        The algorithm is outlined in Section 6.1.3 of [Coh1993]_.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2 * (x^2 + x + 1)^2)
            sage: o = K.maximal_order()
            sage: O = F.maximal_order()
            sage: p = o.ideal(x + 1)
            sage: O.p_radical(p)
            Ideal (x + 1) of Maximal order of Function field in y
            defined by y^3 + x^6 + x^4 + x^2
        """
        from sage.matrix.constructor import matrix
        from sage.modules.free_module_element import vector

        g = prime.gens()[0]

        if not (g.denominator() == 1 and g.numerator().is_irreducible()):
            raise ValueError('not a prime ideal')

        F = self.function_field()
        n = F.degree()
        o = prime.ring()
        p = g.numerator()

        # Fp is isomorphic to the residue field o/p
        Fp, fr_Fp, to_Fp = o._residue_field_global(p)

        # exp = q^j should be at least extension degree where q is
        # the order of the residue field o/p
        q = F.constant_base_field().order()**p.degree()
        exp = q
        while exp <= F.degree():
            exp = exp**q

        # radical equals to the kernel of the map x |-> x^exp
        mat = []
        for g in self.basis():
            v = [to_Fp(c) for c in self._coordinate_vector(g**exp)]
            mat.append(v)
        mat = matrix(Fp, mat)
        ker = mat.kernel()

        # construct module generators of the p-radical
        vecs = []
        for i in range(n):
            v = vector([p if j == i else 0 for j in range(n)])
            vecs.append(v)
        for b in ker.basis():
            v = vector([fr_Fp(c) for c in b])
            vecs.append(v)

        return self._ideal_from_vectors(vecs)

    @cached_method
    def decomposition(self, ideal):
        """
        Return the decomposition of the prime ideal.

        INPUT:

        - ``ideal`` -- prime ideal of the base maximal order

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); R.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: o = K.maximal_order()
            sage: O = F.maximal_order()
            sage: p = o.ideal(x + 1)
            sage: O.decomposition(p)
            [(Ideal (x + 1, y + 1) of Maximal order
             of Function field in y defined by y^3 + x^6 + x^4 + x^2, 1, 1),
             (Ideal (x + 1, (1/(x^3 + x^2 + x))*y^2 + y + 1) of Maximal order
             of Function field in y defined by y^3 + x^6 + x^4 + x^2, 2, 1)]
        """
        from sage.matrix.constructor import matrix

        F = self.function_field()
        n = F.degree()

        p = ideal.gen().numerator()
        o = ideal.ring()

        # Fp is isomorphic to the residue field o/p
        Fp, fr, to = o._residue_field_global(p)
        P,X = Fp['X'].objgen()

        V = Fp**n # Ob = O/pO

        mtable = []
        for i in range(n):
            row = []
            for j in range(n):
                row.append( V([to(e) for e in self._mtable[i][j]]) )
            mtable.append(row)

        if p not in self._kummer_places:
            #####################################
            # Decomposition by Kummer's theorem #
            #####################################
            # gen is self._kummer_gen
            gen_vec_pow = self._kummer_gen_vec_pow
            mul_vecs = self._mul_vecs

            f = self._kummer_polynomial
            fp = P([to(c.numerator()) for c in f.list()])
            decomposition = []
            for q, exp in fp.factor():
                # construct O.ideal([p,q(gen)])
                gen_vecs = list(matrix.diagonal(n * [p]))
                c = q.list()

                # q(gen) in vector form
                qgen = sum(fr(c[i]) * gen_vec_pow[i] for i in range(len(c)))

                I = matrix.identity(o._ring, n)
                for i in range(n):
                    gen_vecs.append(mul_vecs(qgen,I[i]))
                prime = self._ideal_from_vectors_and_denominator(gen_vecs)

                # Compute an element beta in O but not in pO. How to find beta
                # is explained in Section 4.8.3 of [Coh1993]. We keep beta
                # as a vector over k[x] with respect to the basis of O.

                # p and qgen generates the prime; modulo pO, qgenb generates the prime
                qgenb = [to(qgen[i]) for i in range(n)]
                m =[]
                for i in range(n):
                    m.append(sum(qgenb[j] * mtable[i][j] for j in range(n)))
                beta  = [fr(coeff) for coeff in matrix(m).left_kernel().basis()[0]]

                prime.is_prime.set_cache(True)
                prime._prime_below = ideal
                prime._relative_degree = q.degree()
                prime._ramification_index = exp
                prime._beta = beta

                prime._kummer_form = (p, qgen)

                decomposition.append((prime, q.degree(), exp))
        else:
            #############################
            # Buchman-Lenstra algorithm #
            #############################
            from sage.matrix.special import block_matrix
            from sage.modules.free_module_element import vector

            pO = self.ideal(p)
            Ip = self.p_radical(ideal)
            Ob = matrix.identity(Fp, n)

            def bar(I): # transfer to O/pO
                m = []
                for v in I._hnf:
                    m.append([to(e) for e in v])
                h = matrix(m).echelon_form()
                return cut_last_zero_rows(h)

            def liftb(Ib):
                m = [vector([fr(e) for e in v]) for v in Ib]
                m += [v for v in pO._hnf]
                return self._ideal_from_vectors_and_denominator(m,1)

            def cut_last_zero_rows(h):
                i = h.nrows()
                while i > 0 and h.row(i-1).is_zero():
                    i -= 1
                return h[:i]

            def mul_vec(v1,v2):
                s = 0
                for i in range(n):
                    for j in range(n):
                        s += v1[i] * v2[j] * mtable[i][j]
                return s

            def pow(v, r): # r > 0
                m = v
                while r > 1:
                    m = mul_vec(m,v)
                    r -= 1
                return m

            # Algorithm 6.2.7 of [Coh1993]
            def div(Ib, Jb):
                # compute a basis of Jb/Ib
                sJb = Jb.row_space()
                sIb = Ib.row_space()
                sJbsIb,proj_sJbsIb,lift_sJbsIb = sJb.quotient_abstract(sIb)
                supplement_basis = [lift_sJbsIb(v) for v in sJbsIb.basis()]

                m = []
                for b in V.gens(): # basis of Ob = O/pO
                    b_row = [] # row vector representation of the map a -> a*b
                    for a in supplement_basis:
                        b_row += lift_sJbsIb(proj_sJbsIb( mul_vec(a,b) ))
                    m.append(b_row)
                return matrix(Fp,n,m).left_kernel().basis_matrix()

            # Algorithm 6.2.5 of [Coh1993]
            def mul(Ib, Jb):
                m = []
                for v1 in Ib:
                    for v2 in Jb:
                        m.append(mul_vec(v1,v2))
                h = matrix(m).echelon_form()
                return cut_last_zero_rows(h)

            def add(Ib,Jb):
                m = block_matrix([[Ib], [Jb]])
                h = m.echelon_form()
                return cut_last_zero_rows(h)

            # K_1, K_2, ...
            Lb = IpOb = bar(Ip+pO)
            Kb = [Lb]
            while not Lb.is_zero():
                Lb = mul(Lb,IpOb)
                Kb.append(Lb)

            # J_1, J_2, ...
            Jb =[Kb[0]] + [div(Kb[j],Kb[j-1]) for j in range(1,len(Kb))]

            # H_1, H_2, ...
            Hb = [div(Jb[j],Jb[j+1]) for j in range(len(Jb)-1)] + [Jb[-1]]

            q = Fp.order()

            def split(h):
                # VsW represents O/H as a vector space
                W = h.row_space() # H/pO
                VsW,to_VsW,lift_to_V = V.quotient_abstract(W)

                # compute the space K of elements in O/H that satisfy a^q-a=0
                l = [lift_to_V(b) for b in VsW.basis()]

                images = [to_VsW(pow(x, q) - x) for x in l]
                K = VsW.hom(images, VsW).kernel()

                if K.dimension() == 0:
                    return []
                if K.dimension() == 1: # h is prime
                    return [(liftb(h),VsW.dimension())] # relative degree

                # choose a such that a^q - a is 0 but a is not in Fp
                for a in K.basis():
                    # IMPORTANT: This criterion is based on the assumption
                    # that O.basis() starts with 1.
                    if a.support() != [0]:
                        break
                else:
                    raise AssertionError("no appropriate value found")

                a = lift_to_V(a)
                # compute the minimal polynomial of a
                m = [to_VsW(Ob[0])] # 1 in VsW
                apow = a
                while True:
                    v = to_VsW(apow)
                    try:
                        sol = matrix(m).solve_left(v)
                    except ValueError:
                        m.append(v)
                        apow = mul_vec(apow, a)
                        continue
                    break

                minpol = X**len(sol) - P(list(sol))

                # The minimal polynomial of a has only linear factors and at least two
                # of them. We set f to the first factor and g to the product of the rest.
                fac = minpol.factor()
                f = fac[0][0]
                g = (fac/f).expand()
                d,u,v = f.xgcd(g)

                assert d == 1, "Not relatively prime {} and {}".format(f,g)

                # finally, idempotent!
                e = lift_to_V(sum([c1*c2 for c1,c2 in zip(u*f,m)]))

                h1 = add(h, matrix([mul_vec(e,Ob[i]) for i in range(n)]))
                h2 = add(h, matrix([mul_vec(Ob[0]-e,Ob[i]) for i in range(n)]))

                return split(h1) + split(h2)

            decomposition = []
            for i in range(len(Hb)):
                index = i + 1 # Hb starts with H_1
                for prime, degree in split(Hb[i]):
                    # Compute an element beta in O but not in pO. How to find beta
                    # is explained in Section 4.8.3 of [Coh1993]. We keep beta
                    # as a vector over k[x] with respect to the basis of O.
                    m =[]
                    for i in range(n):
                        r = []
                        for g in prime._hnf:
                            r += sum(to(g[j]) * mtable[i][j] for j in range(n))
                        m.append(r)
                    beta = [fr(e) for e in matrix(m).left_kernel().basis()[0]]

                    prime.is_prime.set_cache(True)
                    prime._prime_below = ideal
                    prime._relative_degree = degree
                    prime._ramification_index = index
                    prime._beta = beta

                    decomposition.append((prime, degree, index))

        return decomposition
