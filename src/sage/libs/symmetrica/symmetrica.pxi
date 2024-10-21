# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules

# We put all definitions together, whether they appear in def.h or
# macro.h
cdef extern from 'symmetrica/def.h':
    pass

cdef extern from 'symmetrica/macro.h':
    ctypedef int INT
    ctypedef INT OBJECTKIND

    cdef struct vector:
        pass
    cdef struct bruch:
        pass
    cdef struct graph:
        pass
    cdef struct list:
        pass
    cdef struct matrix:
        pass
    cdef struct monom:
        pass
    cdef struct number:
        pass
    cdef struct partition:
        pass
    cdef struct permutation:
        pass
    cdef struct reihe:
        pass
    cdef struct skewpartition:
        pass
    cdef struct symchar:
        pass
    cdef struct tableaux:
        pass

    cdef enum:
        INFREELIST = -1
        EMPTY = 0
        INTEGER = 1
        VECTOR = 2
        PARTITION = 3
        FRACTION = 4
        BRUCH = 4
        PERMUTATION = 6
        SKEWPARTITION = 7
        TABLEAUX = 8
        POLYNOM = 9
        SCHUR = 10
        MATRIX = 11
        AUG_PART = 12
        HOM_SYM = 13
        HOMSYM = 13
        SCHUBERT = 14
        INTEGERVECTOR = 15
        INTEGER_VECTOR = 15
        INT_VECTOR = 15
        INTVECTOR = 15
        KOSTKA = 16
        INTINT = 17
        SYMCHAR = 18
        WORD =19
        LIST =20 # 210688 */
        MONOM =21 #230688*/
        LONGINT =22 # 170888 */
        GEN_CHAR =23 # 280888  nur fuer test-zwecke */
        BINTREE =24 # 291288 */
        GRAPH =25 # 210889 */
        COMP =26 # 300889 */
        COMPOSITION =26 # 300889 */
        KRANZTYPUS =27 # 280390 */
        POW_SYM =28
        POWSYM =28
        MONOMIAL =29  # 090992 */
        BTREE =30
        KRANZ =31
        GRAL =32 # 200691 */
        GROUPALGEBRA =32 # 170693 */
        ELM_SYM =33  # 090992 */
        ELMSYM =33  # 090992 */
        FINITEFIELD = 35 # 250193 */
        FF = 35 # 250193 */
        REIHE = 36 # 090393 */
        CHARPARTITION = 37 # 130593 */ # internal use */
        CHAR_AUG_PART = 38 # 170593 */ # internal use */
        INTEGERMATRIX =40 # AK 141293 */
        CYCLOTOMIC = 41 # MD */
        MONOPOLY = 42 # MD */
        SQ_RADICAL = 43 # MD */
        BITVECTOR = 44
        LAURENT =45
        SUBSET =47 # AK 220997 */
        FASTPOLYNOM =211093
        EXPONENTPARTITION =240298
        SKEWTABLEAUX =20398
        PARTTABLEAUX =10398
        BARPERM =230695
        PERMVECTOR =180998
        PERM_VECTOR =180998
        PERMUTATIONVECTOR =180998
        PERMUTATION_VECTOR =180998
        INTEGERBRUCH =220998
        INTEGER_BRUCH =220998
        INTEGERFRACTION =220998
        INTEGER_FRACTION =220998
        HASHTABLE  =120199


    cdef struct loc:
        INT w2, w1, w0
        loc *nloc

    cdef struct longint:
        loc *floc
        signed char signum #-1,0,+1
        INT laenge

    ctypedef union OBJECTSELF:
        INT ob_INT
        INT *ob_INTpointer
        char *ob_charpointer
        bruch *ob_bruch
        graph *ob_graph
        list *ob_list
        longint *ob_longint
        matrix *ob_matrix
        monom *ob_monom
        number *ob_number
        partition *ob_partition
        permutation *ob_permutation
        reihe *ob_reihe
        skewpartition *ob_skewpartition
        symchar *ob_symchar
        tableaux *ob_tableaux
        vector *ob_vector

    cdef struct obj:
        OBJECTKIND ob_kind
        OBJECTSELF ob_self

    cdef struct ganzdaten:
        INT basis, basislaenge, auspos, auslaenge, auszz

    cdef struct zahldaten:
        char ziffer[13]
        INT mehr
        INT ziffernzhal
        loc *fdez

    ctypedef obj *OP

    cdef struct vector:
        OP v_length
        OP v_self

    cdef struct REIHE_variablen:
        INT index
        INT potenz
        REIHE_variablen *weiter

    cdef struct REIHE_mon:
        OP coeff
        REIHE_variablen *zeiger
        REIHE_mon *ref

    cdef struct REIHE_poly:
        INT grad
        REIHE_mon *uten
        REIHE_poly *rechts

    cdef struct reihe:
        INT exist
        INT reihenart
        INT z
        reihe *x
        reihe *y
        reihe *p
        INT (*eingabefkt)()
        char ope
        REIHE_poly *infozeig

    ctypedef reihe* REIHE_ZEIGER

    cdef struct list:
        OP l_self
        OP l_next

    cdef struct partition:
        OBJECTKIND pa_kind
        OP pa_self

    cdef struct permutation:
        OBJECTKIND p_kind
        OP p_self

    cdef struct monom:
        OP mo_self
        OP mo_koeff

    cdef struct bruch:
        OP b_oben
        OP b_uten
        INT b_info

    cdef struct matrix:
        OP m_length
        OP m_height
        OP m_self

    cdef struct skewpartition:
        OP spa_gross
        OP spa_klein

    cdef struct tableaux:
        OP t_umriss
        OP t_self

    cdef struct symchar:
        OP sy_werte
        OP sy_parlist
        OP sy_dimension

    cdef struct graph:
        OBJECTKIND gr_kind
        OP gr_self

    cdef struct CYCLO_DATA:
        OP index, deg, poly, autos

    cdef struct FIELD_DATA:
        OP index, deg, poly

    ctypedef union data:
        CYCLO_DATA *c_data
        FIELD_DATA *f_data
        OP o_data

    cdef struct number:
        OP n_self
        data n_data


    #MACROS
    #S_PA_I(OP a, INT i)
    OBJECTKIND s_o_k(OP a)
    void* c_o_k(OP a, OBJECTKIND k)
    OBJECTSELF S_O_S(OP a)

    void* add(OP a, OP b, OP c)
    void* mult(OP a, OP b, OP c)
    void* sub(OP a, OP b, OP c)


    #Integers
    void* m_i_i(INT n, OP a)
    void* M_I_I(INT n, OP a)
    INT S_I_I(OP a)
    t_int_longint(OP a, OP b)
    void* m_i_longint(INT n, OP a)


    #Fractions
    OP S_B_O(OP a)
    OP S_B_U(OP a)
    OP m_ou_b(OP o, OP u, OP d)

    #Vectors
    void* M_IL_V(INT length, OP a)
    void* m_il_v(INT n, OP a )
    void* m_i_i(INT n, OP a)

    INT s_v_li(OP a)
    OP s_v_i(OP a, INT i)

    #Partitions
    OP s_pa_l(OP a)
    INT s_pa_li(OP a)
    INT s_pa_ii(OP a, INT i)
    OP s_pa_i(OP a, INT i)
    OP S_PA_S(OP a)
    OP S_PA_I(OP a, INT )
    void* b_ks_pa(INT kind, OP b, OP a)


    #Skew Partitions
    INT b_gk_spa(OP gross, OP klein, OP result)
    INT m_gk_spa(OP gross, OP klein, OP result)
    OP s_spa_g(OP spa)
    OP s_spa_k(OP spa)

    #Permutations
    OP s_p_i(OP a, INT i)
    INT s_p_li(OP a)
    INT s_p_ii(OP a, INT i)
    void* m_il_p(INT n, OP a)


    #Barred Permutations


    #Lists
    OP s_l_s(OP a)
    OP S_L_S(OP a)
    OP s_l_n(OP a)
    INT lastp_list(OP l)
    INT empty_listp(OP l)

    #Matrices
    INT S_M_HI(OP a)
    INT S_M_LI(OP a)
    OP S_M_IJ(OP a, INT i, INT j)
    void* m_ilih_m(INT l, INT h, OP a)

    #Schur polynomials
    OP s_s_s(OP a)
    OP s_s_k(OP a)
    OP s_s_n(OP a)
    void* m_skn_s(OP part, OP koeff, OP next, OP result)
    void* b_skn_s(OP part, OP koeff, OP next, OP result)

    #Schubert polynomials
    OP s_sch_s(OP a)
    OP s_sch_k(OP a)
    OP s_sch_n(OP a)
    void* m_skn_sch(OP perm, OP koeff, OP next, OP result)
    void* b_skn_sch(OP perm, OP koeff, OP next, OP result)

    #Polynomials
    OP s_po_n(OP a)
    OP s_po_sl(OP a)
    OP s_po_k(OP a)
    OP s_po_s(OP a)
    void* m_skn_po(OP s, OP k, OP next, OP polynom)

    #Tableaux
    OP S_T_S(OP t)

    #########

    INT insert(OP a, OP b, INT (*eq)(), INT (*comp)())


    INT nullp_sqrad(OP a)
    OP S_PO_K(OP a)
    OP S_PO_S(OP a)
    OP S_L_N(OP a)
    OP S_N_S(OP a)
    INT einsp(OP A)

    #Empty Object
    int EMPTYP(OP obj)

    #Functions
    INT anfang()
    INT ende()
    OP callocobject()
    INT sscan(INT, OP a)
    INT scan(INT, OP a)
    INT freeall(OP a)
    INT freeself(OP a)
    INT sprint(char* t, OP a)
    INT sprint_integer(char* t, OP A)
    INT println(OP a)

    #factorial
    INT fakul(OP a, OP b)

##########################################
cdef object matrix_constructor
cdef object Integer
cdef object Tableau, SkewTableau
cdef object SkewPartition
cdef object Partition
cdef object Permutation, Permutations
cdef object builtinlist
cdef object sqrt
cdef object Rational
cdef object QQ
cdef object ZZ
cdef object SymmetricFunctions
cdef object PolynomialRing
cdef object SchubertPolynomialRing, SchubertPolynomial_class
cdef object two, fifteen, thirty, zero, sage_maxint

cdef int maxint = 2147483647

cdef void late_import() noexcept:
    global matrix_constructor, \
           Integer, \
           Tableau, \
           SkewTableau, \
           SkewPartition, \
           Partition, \
           Permutation, Permutations,\
           prod, \
           PolynomialRing, \
           Rational, \
           QQ, \
           ZZ, \
           SymmetricFunctions, \
           sqrt, \
           builtinlist, \
           MPolynomialRing_base, MPolynomial,\
           SchubertPolynomialRing, SchubertPolynomial_class,\
           two, fifteen, thirty, zero, sage_maxint

    if Integer is not None:
        return

    import sage.rings.integer
    Integer = sage.rings.integer.Integer

    import sage.combinat.tableau
    Tableau = sage.combinat.tableau.Tableau

    import sage.combinat.skew_tableau
    SkewTableau = sage.combinat.skew_tableau.SkewTableau

    import sage.combinat.skew_partition
    SkewPartition = sage.combinat.skew_partition.SkewPartition

    import sage.combinat.partition
    Partition = sage.combinat.partition.Partition

    import sage.combinat.permutation
    Permutation = sage.combinat.permutation.Permutation
    Permutations = sage.combinat.permutation.Permutations

    from sage.misc.functional import sqrt
    from sage.misc.misc_c import prod
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    from sage.rings.rational_field import QQ
    from sage.rings.rational import Rational
    from sage.rings.integer_ring import ZZ

    import builtins
    builtinlist = builtins.list

    import sage.rings.polynomial.multi_polynomial_ring
    MPolynomialRing_base = sage.rings.polynomial.multi_polynomial_ring.MPolynomialRing_base
    import sage.rings.polynomial.multi_polynomial
    MPolynomial = sage.rings.polynomial.multi_polynomial.MPolynomial

    two = Integer(2)
    fifteen = Integer(15)
    thirty = Integer(30)
    zero = Integer(0)
    sage_maxint = Integer(maxint)


    #Symmetric Functions
    import sage.combinat.sf.sf
    SymmetricFunctions = sage.combinat.sf.sf.SymmetricFunctions

    import sage.matrix.constructor
    matrix_constructor = sage.matrix.constructor.matrix

    import sage.combinat.schubert_polynomial
    SchubertPolynomialRing = sage.combinat.schubert_polynomial.SchubertPolynomialRing
    SchubertPolynomial_class = sage.combinat.schubert_polynomial.SchubertPolynomial_class


##########################################
cdef object _py(OP a):
    cdef OBJECTKIND objk
    objk = s_o_k(a)
    if objk == INTEGER:
        return _py_int(a)
    elif objk == LONGINT:
        return _py_longint(a)
    elif objk == PARTITION:
        return _py_partition(a)
    elif objk == PERMUTATION:
        return _py_permutation(a)
    elif objk == SKEWPARTITION:
        return _py_skew_partition(a)
    elif objk == FRACTION:
        return _py_fraction(a)
    elif objk == SQ_RADICAL:
        return _py_sq_radical(a)
    elif objk == MATRIX or objk == KRANZTYPUS:
        return _py_matrix(a)
    elif objk == SCHUR:
        return _py_schur(a)
    elif objk == HOMSYM:
        return _py_homsym(a)
    elif objk == POWSYM:
        return _py_powsym(a)
    elif objk == ELMSYM:
        return _py_elmsym(a)
    elif objk == MONOMIAL:
        return _py_monomial(a)
    elif objk == LIST:
        return _py_list(a)
    elif objk == VECTOR:
        return _py_vector(a)
    elif objk == TABLEAUX:
        return _py_tableau(a)
    elif objk == EMPTY:
        return None
    elif objk == POLYNOM:
        return _py_polynom(a)
    elif objk == SCHUBERT:
        return _py_schubert(a)
    else:
        #println(a)
        raise NotImplementedError(str(objk))

cdef int _op(object a, OP result) except -1:
    late_import()
    if isinstance(a, Integer):
        _op_integer(a, result)
    elif isinstance(a, Partition):
        _op_partition(a, result)
    elif isinstance(a, Rational):
        _op_fraction(a, result)
    else:
        raise TypeError("cannot convert a (= %s) to OP" % a)


def test_integer(object x):
    """
    Test functionality for converting between Sage's integers
    and symmetrica's integers.

    EXAMPLES::

        sage: from sage.libs.symmetrica.symmetrica import test_integer
        sage: test_integer(1)
        1
        sage: test_integer(-1)
        -1
        sage: test_integer(2^33)
        8589934592
        sage: test_integer(-2^33)
        -8589934592
        sage: test_integer(2^100)
        1267650600228229401496703205376
        sage: test_integer(-2^100)
        -1267650600228229401496703205376
        sage: for i in range(100):
        ....:     if test_integer(2^i) != 2^i:
        ....:         print("Failure at {}".format(i))
    """
    cdef OP a = callocobject()
    _op_integer(x, a)
    res = _py(a)
    freeall(a)
    return res

##########
#Integers#
##########

cdef int _op_integer(object x, OP a) except -1:
    try:
        _op_int(x, a)
    except OverflowError:
        _op_longint(x, a)
    return 0


cdef int _op_int(object x, OP a) except -1:
    M_I_I(x, a)
    return 0

cdef object _py_int(OP a):
    late_import()
    return Integer(S_I_I(a))


cdef int _op_longint(object x, OP a) except -1:
    late_import()
    cdef OP op_maxint_long = callocobject(),
    cdef OP quo_long = callocobject()
    cdef OP rem_long = callocobject()

    qr = x.quo_rem(sage_maxint)

    m_i_longint(maxint, op_maxint_long)
    _op_integer(qr[0], a)
    _op_integer(qr[1], rem_long)

    #Multiply a by op_maxint_long
    mult(op_maxint_long, a, a)

    #Add rem to a
    add(a, rem_long, a)

    freeall(rem_long)
    freeall(quo_long)
    freeall(op_maxint_long)
    return 0

cdef object _py_longint(OP a):
    late_import()
    cdef longint *x = S_O_S(a).ob_longint
    cdef loc *l = x.floc
    cdef int sign = x.signum
    res = zero
    n = zero
    while l != NULL:
        res += Integer( l.w0 ) * two**n
        res += Integer( l.w1 ) * two**(n+fifteen)
        res += Integer( l.w2 ) * two**(n+thirty)
        n += thirty + fifteen
        l = l.nloc
    if sign < 0:
        res *= Integer(-1)

    return res


###########
#Fractions#
###########
cdef object _py_fraction(OP a):
    return _py(S_B_O(a))/_py(S_B_U(a))

cdef int _op_fraction(object f, OP a) except -1:
    cdef OP o = callocobject(), u = callocobject()
    _op_integer(f.numerator(), o)
    _op_integer(f.denominator(), u)
    m_ou_b(o, u, a)

#########
#Vectors#
#########
cdef object _py_vector(OP a):
    cdef INT i
    res = []
    for i from 0 <= i < s_v_li(a):
        res.append( _py(s_v_i(a, i)))
    return res

cdef void* _op_il_vector(object l, OP a) noexcept:
    cdef INT length, i
    length = len(l)

    m_il_v(length, a)
    for i from 0 <= i < length:
        m_i_i(l[i], s_v_i(a, i))

#########
#Numbers#
#########
cdef object _py_sq_radical(OP a):
    late_import()

    cdef OP ptr
    ptr = S_N_S(a)

    res = 0
    if nullp_sqrad(a):
        return res

    while ptr != NULL:

        if einsp(S_PO_S(ptr)):
            res += _py(S_PO_K(ptr))
        else:
            res += _py(S_PO_K(ptr))*sqrt(_py(S_PO_S(ptr)))

        ptr = S_L_N(ptr)

    return res.radical_simplify()

############
#Partitions#
############
cdef void* _op_partition(object p, OP a) noexcept:
    cdef int n, i, j

    if not EMPTYP(a):
        freeself(a)

    n = len(p)
    b_ks_pa(VECTOR, callocobject(), a)
    m_il_v(n, S_PA_S(a))

    j = 0
    for i from n > i >= 0:
        _op_integer(p[i], S_PA_I(a,j))
        j = j + 1

cdef object _py_partition(OP a):
    cdef INT n, i
    late_import()
    res = []
    n = s_pa_li(a)
    for i from n > i >=0:
        res.append(s_pa_ii(a, i))
    return Partition(res)

################
#Skew Partition#
################
cdef void* _op_skew_partition(object p, OP a) noexcept:
    cdef OP gross, klein
    gross = callocobject()
    klein = callocobject()
    _op_partition(p[0], gross)
    _op_partition(p[1], klein)
    b_gk_spa(gross, klein, a)

cdef object _py_skew_partition(OP a):
    late_import()
    return SkewPartition( [ _py_partition(s_spa_g(a)), _py_partition(s_spa_k(a)) ] )

##############
#Permutations#
##############
cdef void* _op_permutation(object p, OP a) noexcept:
    cdef int n, i, j

    if not EMPTYP(a):
        freeself(a)

    n = len(p)
    m_il_p(n, a)
    for i from 0 <= i < n:
        _op_integer(p[i], s_p_i(a,i))

cdef object _py_permutation(OP a):
    late_import()
    cdef INT n, i
    res = []
    n = s_p_li(a)
    for i from 0 <= i < n:
        res.append(s_p_ii(a, i))
    return Permutation(res)

#####################
#Barred Permutations#
#####################

#######
#Lists#
#######
cdef object _py_list(OP a):
    cdef OP x
    x = a
    res = []
    if S_L_S(a) == NULL:
        return []
    elif empty_listp(a):
        return []
    while x  != NULL:
        res.append(_py(s_l_s(x)))
        x = s_l_n(x)

    return res


#############
#Polynomials#
#############
cdef object _py_polynom(OP a):
    late_import()
    cdef int maxneeded = 0, i = 0
    cdef OP pointer = a

    if pointer == NULL:
        return 0


    #Find the maximum number of variables needed
    while pointer != NULL:
        l = _py(s_po_sl(pointer))
        if l > maxneeded:
            maxneeded = l
        pointer = s_po_n(pointer)

    pointer = a
    parent_ring = _py(s_po_k(pointer)).parent()
    if maxneeded == 1:
        P = PolynomialRing(parent_ring, 'x')
    else:
        P = PolynomialRing(parent_ring, maxneeded, 'x')
    d = {}
    while pointer != NULL:
        exps = tuple(_py_vector(s_po_s(pointer)))
        d[ exps ] = _py(s_po_k(pointer))
        pointer = s_po_n(pointer)

    return P(d)


cdef object _py_polynom_alphabet(OP a, object alphabet, object length):
    """
    Convert a symmetrica multivariate polynomial a to a Sage multivariate
    polynomials.  Alphabet specifies the names of the variables which are
    fed into PolynomialRing.  length specifies the number of variables; if
    it is set to 0, then the number of variables is autodetected based on
    the number of variables in alphabet or the result obtained from
    symmetrica.
    """
    late_import()
    cdef OP pointer = a

    if pointer == NULL:
        return 0

    parent_ring = _py(s_po_k(pointer)).parent()
    if length == 0:
        if isinstance(alphabet, (builtinlist, tuple)):
            l = len(alphabet)
        elif isinstance(alphabet, str) and ',' in alphabet:
            l = len(alphabet.split(','))
        else:
            l = _py(s_po_sl(a))
    else:
        l = length

    P = PolynomialRing(parent_ring, l, alphabet)
    x = P.gens()
    res = P(0)
    while pointer != NULL:
        exps = _py_vector(s_po_s(pointer))
        res += _py(s_po_k(pointer)) *prod([ x[i]**exps[i] for i in range(min(len(exps),l))])
        pointer = s_po_n(pointer)
    return res

cdef object _op_polynom(object d, OP res):
    late_import()

    poly_ring = d.parent()

    if not isinstance(poly_ring, MPolynomialRing_base):
        raise TypeError("you must pass a multivariate polynomial")
    base_ring = poly_ring.base_ring()

    if not ( base_ring == ZZ or base_ring == QQ):
        raise TypeError("the base ring must be either ZZ or QQ")

    cdef OP c = callocobject(), v = callocobject()
    cdef OP pointer = res
    m = d.monomials()
    exp = d.exponents()
    cdef int n, i
    n = len(exp)

    for i from 0 <= i < n:
        _op_il_vector(exp[i], v)
        _op(d.monomial_coefficient(poly_ring(m[i])), c)
        if i != n-1:
            m_skn_po(v,c, callocobject(), pointer)
        else:
            m_skn_po(v,c, NULL, pointer)
        pointer = s_po_n(pointer)

    freeall(c)
    freeall(v)
    return None


#######################################
#Schur symmetric functions and friends#
#######################################
cdef object _py_schur(OP a):
    late_import()
    z_elt = _py_schur_general(a)
    if len(z_elt) == 0:
        return SymmetricFunctions(ZZ).s()(0)

    #Figure out the parent ring of a coefficient
    R = z_elt[next(iter(z_elt))].parent()

    s = SymmetricFunctions(R).s()
    z = s(0)
    z._monomial_coefficients = z_elt
    return z

cdef void* _op_schur(object d, OP res) noexcept:
    _op_schur_general(d, res)

cdef object _py_monomial(OP a): #Monomial symmetric functions
    late_import()
    z_elt = _py_schur_general(a)
    if len(z_elt) == 0:
        return SymmetricFunctions(ZZ).m()(0)

    R = z_elt[next(iter(z_elt))].parent()

    m = SymmetricFunctions(R).m()
    z = m(0)
    z._monomial_coefficients = z_elt
    return z

cdef void* _op_monomial(object d, OP res) noexcept: #Monomial symmetric functions
    cdef OP pointer = res
    _op_schur_general(d, res)
    while pointer != NULL:
        c_o_k(pointer, MONOMIAL)
        pointer = s_s_n(pointer)

cdef object _py_powsym(OP a):  #Power-sum symmetric functions
    late_import()
    z_elt = _py_schur_general(a)
    if len(z_elt) == 0:
        return SymmetricFunctions(ZZ).p()(0)

    R = z_elt[next(iter(z_elt))].parent()

    p = SymmetricFunctions(R).p()
    z = p(0)
    z._monomial_coefficients = z_elt
    return z

cdef void* _op_powsym(object d, OP res) noexcept: #Power-sum symmetric functions
    cdef OP pointer = res
    _op_schur_general(d, res)
    while pointer != NULL:
        c_o_k(pointer, POWSYM)
        pointer = s_s_n(pointer)


cdef object _py_elmsym(OP a): #Elementary symmetric functions
    late_import()
    z_elt = _py_schur_general(a)
    if len(z_elt) == 0:
        return SymmetricFunctions(ZZ).e()(0)

    R = z_elt[next(iter(z_elt))].parent()

    e = SymmetricFunctions(R).e()
    z = e(0)
    z._monomial_coefficients = z_elt
    return z

cdef void* _op_elmsym(object d, OP res) noexcept: #Elementary symmetric functions
    cdef OP pointer = res
    _op_schur_general(d, res)
    while pointer != NULL:
        c_o_k(pointer, ELMSYM)
        pointer = s_s_n(pointer)


cdef object _py_homsym(OP a):  # Homogeneous symmetric functions
    late_import()
    z_elt = _py_schur_general(a)
    if len(z_elt) == 0:
        return SymmetricFunctions(ZZ).h()(0)

    R = z_elt[next(iter(z_elt))].parent()

    h = SymmetricFunctions(R).h()
    z = h(0)
    z._monomial_coefficients = z_elt
    return z

cdef void* _op_homsym(object d, OP res) noexcept:  # Homogeneous symmetric functions
    cdef OP pointer = res
    _op_schur_general(d, res)
    while pointer != NULL:
        c_o_k(pointer, HOMSYM)
        pointer = s_s_n(pointer)


cdef object _py_schur_general(OP a):
    cdef OP pointer = a
    d = {}
    if a == NULL:
        return d
    while pointer != NULL:
        d[ _py_partition(s_s_s(pointer)) ] = _py(s_s_k(pointer))
        pointer = s_s_n(pointer)
    return d

cdef void* _op_schur_general(object d, OP res) noexcept:
    if isinstance(d, dict):
        _op_schur_general_dict(d, res)
    else:
        _op_schur_general_sf(d, res)

cdef void* _op_schur_general_sf(object f, OP res) noexcept:
    late_import()
    base_ring = f.parent().base_ring()
    if not ( base_ring is QQ or base_ring is ZZ ):
        raise ValueError("the base ring must be either ZZ or QQ")

    _op_schur_general_dict( f.monomial_coefficients(), res)

cdef void* _op_schur_general_dict(object d, OP res) noexcept:
    late_import()

    cdef OP next
    cdef OP pointer = res
    cdef INT n, i


    keys = builtinlist(d)
    n = len(keys)

    if n == 0:
        raise ValueError("the dictionary must be nonempty")

    b_skn_s(callocobject(), callocobject(), NULL, res)
    _op_partition(keys[0], s_s_s(res))
    _op(d[keys[0]], s_s_k(res))


    for i from 0 < i < n:
        next = callocobject()

        b_skn_s(callocobject(), callocobject(), NULL, next)
        _op_partition(keys[i], s_s_s(next))
        _op(d[keys[i]], s_s_k(next))

        insert(next, res, NULL, NULL)


######################
#Schubert Polynomials#
######################
cdef void* _op_schubert_general(object d, OP res) noexcept:
    if isinstance(d, dict):
        _op_schubert_dict(d, res)
    else:
        _op_schubert_sp(d, res)

cdef void* _op_schubert_perm(object a, OP res) noexcept:
    cdef OP caperm = callocobject()
    _op_permutation(a, caperm)
    m_perm_sch(caperm, res)
    freeall(caperm)

cdef void* _op_schubert_sp(object f, OP res) noexcept:
    late_import()
    base_ring = f.parent().base_ring()
    if not ( base_ring is QQ or base_ring is ZZ ):
        raise ValueError("the base ring must be either ZZ or QQ")

    _op_schubert_dict( f.monomial_coefficients(), res)

cdef void* _op_schubert_dict(object d, OP res) noexcept:
    late_import()

    cdef OP next
    cdef OP pointer = res
    cdef INT n, i

    keys = builtinlist(d)
    n = len(keys)

    if n == 0:
        raise ValueError("the dictionary must be nonempty")

    b_skn_sch(callocobject(), callocobject(), NULL, res)
    _op_permutation(keys[0], s_sch_s(res))
    _op(d[keys[0]], s_sch_k(res))


    for i from 0 < i < n:
        next = callocobject()

        b_skn_sch(callocobject(), callocobject(), NULL, next)
        _op_permutation(keys[i], s_sch_s(next))
        _op(d[keys[i]], s_sch_k(next))

        insert(next, res, NULL, NULL)

cdef object _py_schubert(OP a):
    late_import()
    cdef OP pointer = a
    cdef dict z_elt = {}

    # SCHUBERT is (like) a list, so we also need to make sure it is not empty
    if a == NULL or empty_listp(a):
        return SchubertPolynomialRing(ZZ).zero()

    while pointer != NULL:
        z_elt[_py(s_s_s(pointer)).remove_extra_fixed_points() ] = _py(s_sch_k(pointer))
        pointer = s_sch_n(pointer)

    if not z_elt:
        return SchubertPolynomialRing(ZZ).zero()

    R = z_elt[next(iter(z_elt))].parent()
    X = SchubertPolynomialRing(R)
    # If the element constructor ends up copying the input dict in the future,
    #   then this will not be as fast as creating a copy of the zero element
    #   and explicitly setting the _monomial_coefficients.
    return X.element_class(X, z_elt)


##########
#Matrices#
##########
cdef object _py_matrix(OP a):

    late_import()

    cdef INT i,j,rows, cols
    rows = S_M_HI(a)
    cols = S_M_LI(a)

    res = []
    for i from 0 <= i < rows:
        row = []
        for j from 0 <= j < cols:
            row.append( _py(S_M_IJ(a,i,j)) )

        res.append(row)

    #return res
    if res == [] or res is None:
        return res
    else:
        return matrix_constructor(res)


cdef void* _op_matrix(object a, OP res) noexcept:
    #FIXME: only constructs integer matrices

    cdef INT i,j,rows, cols

    rows = a.nrows()
    cols = a.ncols()

    m_ilih_m(rows, cols, res)

    for i from 0 <= i < rows:
        for j from 0 <= j < cols:
            _op_integer( a[(i,j)], S_M_IJ(res,i,j) )

##########
#Tableaux#
##########
cdef object _py_tableau(OP t):

    late_import()

    cdef INT i,j,rows, cols, added, is_skew = 0
    cdef OP a
    a = S_T_S(t)
    rows = S_M_HI(a)
    cols = S_M_LI(a)

    res = []
    for i from 0 <= i < rows:
        row = []
        added = 0
        for j from 0 <= j < cols:
            if s_o_k(S_M_IJ(a,i,j)) == EMPTY:
                if added:
                    break
                else:
                    row.append( None )
                    is_skew = 1
            else:
                added = 1
                row.append( _py(S_M_IJ(a,i,j)) )

        res.append(row)

    #return res
    if is_skew:
        return SkewTableau(res)
    else:
        return Tableau(res)


def start():
    anfang()

    # Automatically clean up when sage exits.
    import atexit
    atexit.register(ende)
