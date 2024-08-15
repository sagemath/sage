cdef extern from 'symmetrica/def.h':
    INT mult_schubert_schubert(OP a, OP b, OP result)
    INT m_perm_sch(OP a, OP b)
    INT t_SCHUBERT_POLYNOM(OP a, OP b)
    INT t_POLYNOM_SCHUBERT(OP a, OP b)
    INT mult_schubert_variable(OP a, OP i, OP r)
    INT divdiff_perm_schubert(OP perm, OP sb, OP res)
    INT scalarproduct_schubert(OP a, OP b, OP c)
    INT divdiff_schubert(OP a, OP schub, OP res)

    INT t_2SCHUBERT_POLYNOM(OP a,OP b)
    INT mult_schubert_polynom(OP a,OP b,OP c)


cdef object _check_schubert(object a, OP ca):
    if a in Permutations():
        if isinstance(a, builtinlist):
            a = Permutation(a)
        _op_schubert_perm(a, ca)
        return max(a.reduced_word()+[0])
    elif isinstance(a, SchubertPolynomial_class):
        br = a.parent().base_ring()
        if (br == QQ or br == ZZ):
            _op_schubert_sp(a, ca)
            return min(max(i.reduced_word(), default=0) for i in a.support())
        else:
            raise ValueError("a must be a Schubert polynomial over ZZ or QQ")
    else:
        raise TypeError("a must be a permutation or a Schubert polynomial")


def mult_schubert_schubert_symmetrica(a, b):
    """
    Multiply the Schubert polynomials `a` and `b`.

    EXAMPLES::

        sage: symmetrica.mult_schubert_schubert([3,2,1], [3,2,1])
        X[5, 3, 1, 2, 4]
    """
    late_import()

    cdef OP ca = callocobject(), cb = callocobject(), cres = callocobject()

    try:
        max_a = _check_schubert(a, ca)
        max_b = _check_schubert(b, cb)
    except (ValueError, TypeError), err:
        freeall(ca)
        freeall(cb)
        freeall(cres)
        raise err


    sig_on()
    mult_schubert_schubert(ca, cb, cres)
    sig_off()

    res = _py(cres)

    freeall(ca)
    freeall(cb)
    freeall(cres)

    return res


def t_SCHUBERT_POLYNOM_symmetrica(a):
    """
    Convert a Schubert polynomial to a 'regular' multivariate
    polynomial.

    EXAMPLES::

        sage: symmetrica.t_SCHUBERT_POLYNOM([3,2,1])
        x0^2*x1
    """
    late_import()

    cdef OP ca = callocobject(), cres = callocobject()

    try:
        max_a = _check_schubert(a, ca)
    except (ValueError, TypeError), err:
        freeall(ca)
        freeall(cres)
        raise err

    sig_on()
    t_SCHUBERT_POLYNOM(ca, cres)
    sig_off()

    res = _py(cres)

    freeall(ca)
    freeall(cres)

    return res


def t_POLYNOM_SCHUBERT_symmetrica(a):
    """
    Convert a multivariate polynomial a to a Schubert polynomial.

    EXAMPLES::

        sage: R.<x1,x2,x3> = QQ[]
        sage: w0 = x1^2*x2
        sage: symmetrica.t_POLYNOM_SCHUBERT(w0)
        X[3, 2, 1]
    """
    late_import()

    cdef OP ca = callocobject(), cres = callocobject()

    if not isinstance(a, MPolynomial):
        freeall(ca)
        freeall(cres)
        raise TypeError("a (= %s) must be a multivariate polynomial")
    else:
        br = a.parent().base_ring()
        if br != QQ and br != ZZ:
            freeall(ca)
            freeall(cres)
            raise ValueError("a's base ring must be either ZZ or QQ")
        else:
            _op_polynom(a, ca)

    sig_on()
    t_POLYNOM_SCHUBERT(ca, cres)
    sig_off()

    res = _py(cres)

    freeall(ca)
    freeall(cres)

    return res


def mult_schubert_variable_symmetrica(a, i):
    """
    Return the product of `a` and `x_i`.  Note that indexing with `i`
    starts at 1.

    EXAMPLES::

        sage: symmetrica.mult_schubert_variable([3,2,1], 2)
        X[3, 2, 4, 1]
        sage: symmetrica.mult_schubert_variable([3,2,1], 4)
        X[3, 2, 1, 4, 6, 5] - X[3, 2, 1, 5, 4]
    """
    late_import()

    cdef OP ca = callocobject(), ci = callocobject(),  cres = callocobject()

    try:
        max_a = _check_schubert(a, ca)
    except (ValueError, TypeError), err:
        freeall(ca)
        freeall(ci)
        freeall(cres)
        raise err

    _op_integer(i, ci)

    sig_on()
    mult_schubert_variable(ca, ci, cres)
    sig_off()

    res = _py(cres)

    freeall(ca)
    freeall(ci)
    freeall(cres)

    return res


def divdiff_perm_schubert_symmetrica(perm, a):
    r"""
    Return the result of applying the divided difference operator
    `\delta_i` to `a` where `a` is either a permutation or a
    Schubert polynomial over `\QQ`.

    EXAMPLES::

       sage: symmetrica.divdiff_perm_schubert([2,3,1], [3,2,1])
       X[2, 1]
       sage: symmetrica.divdiff_perm_schubert([3,1,2], [3,2,1])
       X[1, 3, 2]
       sage: symmetrica.divdiff_perm_schubert([3,2,4,1], [3,2,1])
       Traceback (most recent call last):
       ...
       ValueError: cannot apply \delta_{[3, 2, 4, 1]} to a (= [3, 2, 1])
    """
    late_import()

    cdef OP ca = callocobject(), cperm = callocobject(),  cres = callocobject()

    try:
        max_a = _check_schubert(a, ca)
    except (ValueError, TypeError), err:
        freeall(ca)
        freeall(cperm)
        freeall(cres)
        raise err

    if perm not in Permutations():
        freeall(ca)
        freeall(cperm)
        freeall(cres)
        raise TypeError("perm must be a permutation")
    else:
        perm = Permutation(perm)
        rw = perm.reduced_word()
        max_perm = max(rw)
        _op_permutation(perm, cperm)

    if max_perm > max_a:
        freeall(ca)
        freeall(cperm)
        freeall(cres)
        raise ValueError(r"cannot apply \delta_{%s} to a (= %s)" % (perm, a))

    sig_on()
    divdiff_perm_schubert(cperm, ca, cres)
    sig_off()

    res = _py(cres)

    freeall(ca)
    freeall(cperm)
    freeall(cres)

    return res


def scalarproduct_schubert_symmetrica(a, b):
    """
    EXAMPLES::

        sage: symmetrica.scalarproduct_schubert([3,2,1], [3,2,1])
        X[1, 3, 5, 2, 4]
        sage: symmetrica.scalarproduct_schubert([3,2,1], [2,1,3])
        X[1, 2, 4, 3]
    """
    late_import()

    cdef OP ca = callocobject(), cb = callocobject(), cres = callocobject()

    try:
        max_a = _check_schubert(a, ca)
        max_b = _check_schubert(b, cb)
    except (ValueError, TypeError), err:
        freeall(ca)
        freeall(cb)
        freeall(cres)
        raise err

    sig_on()
    scalarproduct_schubert(ca, cb, cres)
    sig_off()

    if empty_listp(cres):
        res = Integer(0)
    else:
        res = _py(cres)

    freeall(ca)
    freeall(cb)
    freeall(cres)

    return res


def divdiff_schubert_symmetrica(i, a):
    r"""
    Return the result of applying the divided difference operator
    `\delta_i` to `a` where `a` is either a permutation or a
    Schubert polynomial over `\QQ`.

    EXAMPLES::

       sage: symmetrica.divdiff_schubert(1, [3,2,1])
       X[2, 3, 1]
       sage: symmetrica.divdiff_schubert(2, [3,2,1])
       X[3, 1, 2]
       sage: symmetrica.divdiff_schubert(3, [3,2,1])
       Traceback (most recent call last):
       ...
       ValueError: cannot apply \delta_{3} to a (= [3, 2, 1])
    """
    late_import()

    cdef OP ca = callocobject(), ci = callocobject(),  cres = callocobject()

    try:
        max_a = _check_schubert(a, ca)
    except (ValueError, TypeError), err:
        freeall(ca)
        freeall(ci)
        freeall(cres)
        raise err

    if not isinstance(i, (int, Integer)):
        freeall(ca)
        freeall(ci)
        freeall(cres)
        raise TypeError("i must be an integer")
    else:
        _op_integer(i, ci)

    if i > max_a or i <= 0:
        freeall(ca)
        freeall(ci)
        freeall(cres)
        raise ValueError(r"cannot apply \delta_{%s} to a (= %s)" % (i, a))

    sig_on()
    divdiff_schubert(ci, ca, cres)
    sig_off()

    res = _py(cres)

    freeall(ca)
    freeall(ci)
    freeall(cres)

    return res
