cdef extern from 'symmetrica/def.h':
    INT outerproduct_schur(OP parta, OP partb, OP result)
    INT dimension_schur(OP a, OP result)
    INT part_part_skewschur(OP big, OP small, OP result)
    INT newtrans(OP perm, OP schur)
    INT compute_schur_with_alphabet(OP part, OP length, OP poly)
    INT compute_homsym_with_alphabet(OP number, OP length, OP poly)
    INT compute_elmsym_with_alphabet(OP number, OP length, OP poly)
    INT compute_monomial_with_alphabet(OP partition, OP length, OP poly)
    INT compute_powsym_with_alphabet(OP number, OP length, OP poly)
    INT compute_schur_with_alphabet_det(OP part, OP length, OP poly)

    INT part_part_skewschur(OP a, OP b, OP c)

    INT t_SCHUR_MONOMIAL(OP schur, OP result)
    INT t_SCHUR_HOMSYM(OP a, OP b)
    INT t_SCHUR_ELMSYM(OP a, OP b)

    INT t_MONOMIAL_SCHUR(OP a, OP b)
    INT t_MONOMIAL_HOMSYM(OP a, OP b)
    INT t_MONOMIAL_ELMSYM(OP a, OP b)

    INT t_ELMSYM_SCHUR(OP a, OP b)
    INT t_ELMSYM_MONOMIAL(OP a, OP b)
    INT t_ELMSYM_HOMSYM(OP a, OP b)

    INT t_HOMSYM_SCHUR(OP a, OP b)
    INT t_HOMSYM_MONOMIAL(OP a, OP b)
    INT t_HOMSYM_ELMSYM(OP a, OP b)

    INT t_POWSYM_SCHUR(OP a, OP b)
    INT t_SCHUR_POWSYM(OP a, OP b)
    INT t_POWSYM_HOMSYM(OP a, OP b)
    INT t_HOMSYM_POWSYM(OP a, OP b)
    INT t_POWSYM_ELMSYM(OP a, OP b)
    INT t_ELMSYM_POWSYM(OP a, OP b)
    INT t_POWSYM_MONOMIAL(OP a, OP b)
    INT t_MONOMIAL_POWSYM(OP a, OP b)

    INT hall_littlewood(OP part, OP res)

    INT mult_schur_schur(OP s1, OP s2, OP res)
    INT mult_monomial_monomial(OP m1, OP m2, OP res)

    INT t_POLYNOM_POWER(OP a, OP b)
    INT t_POLYNOM_SCHUR(OP a, OP b)
    INT t_POLYNOM_ELMSYM(OP a, OP b)
    INT t_POLYNOM_MONOMIAL(OP a, OP b)

    INT symmetricp(OP a)

    INT scalarproduct_schur(OP a, OP b, OP c)


def outerproduct_schur_symmetrica(parta, partb):
    """
    You enter two PARTITION objects, and the result is
    a SCHUR object, which is the expansion of the product
    of the two schurfunctions, labeled by
    the two PARTITION objects parta and partb.
    Of course this can also be interpreted as the decomposition of the
    outer tensor product of two irreducible representations of the
    symmetric group.

    EXAMPLES::

        sage: symmetrica.outerproduct_schur([2],[2])
        s[2, 2] + s[3, 1] + s[4]
    """
    cdef OP cparta, cpartb, cresult

    cparta  = callocobject()
    cpartb  = callocobject()
    cresult = callocobject()

    _op_partition(parta, cparta)
    _op_partition(partb, cpartb)

    sig_on()
    outerproduct_schur(cparta, cpartb, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cparta)
    freeall(cpartb)
    freeall(cresult)

    return res


def dimension_schur_symmetrica(s):
    """
    You enter a SCHUR object a, and the result is the
    dimension of the corresponding representation of the
    symmetric group sn.
    """
    cdef OP ca, cresult

    cresult = callocobject()
    ca      = callocobject()

    _op_schur(s, ca)
    sig_on()
    dimension_schur(ca, cresult)
    sig_off()
    res = _py(cresult)

    freeall(ca)
    freeall(cresult)

    return res


def newtrans_symmetrica(perm):
    """
    Compute the decomposition of a schubertpolynomial labeled by
    the permutation perm, as a sum of Schurfunction.

    FIXME!
    """
    cdef OP cperm = callocobject(), cresult = callocobject()

    _op_permutation(perm, cperm)

    sig_on()
    newtrans(cperm, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(cperm)

    return res


def compute_schur_with_alphabet_symmetrica(part, length, alphabet='x'):
    """
    Compute the expansion of a schurfunction labeled by a
    partition PART as a POLYNOM erg. The INTEGER length specifies the
    length of the alphabet.

    EXAMPLES::

        sage: symmetrica.compute_schur_with_alphabet(2,2)
        x0^2 + x0*x1 + x1^2
        sage: symmetrica.compute_schur_with_alphabet([2],2)
        x0^2 + x0*x1 + x1^2
        sage: symmetrica.compute_schur_with_alphabet(Partition([2]),2)
        x0^2 + x0*x1 + x1^2
        sage: symmetrica.compute_schur_with_alphabet(Partition([2]),2,'y')
        y0^2 + y0*y1 + y1^2
        sage: symmetrica.compute_schur_with_alphabet(Partition([2]),2,'a,b')
        a^2 + a*b + b^2
        sage: symmetrica.compute_schur_with_alphabet([2,1],1,'x')
        0
    """
    late_import()
    cdef OP cpart = callocobject(), cresult = callocobject(), clength = callocobject()

    if isinstance(part, (int, Integer)):
        _op_partition([part], cpart)
    elif isinstance(part, (builtinlist, Partition)):
        _op_partition(part, cpart)
    else:
        raise NotImplementedError("n must be an integer or partition")
    _op_integer(length, clength)

    sig_on()
    compute_schur_with_alphabet(cpart, clength, cresult)
    sig_off()

    res = _py_polynom_alphabet(cresult, alphabet, length)

    freeall(cresult)
    freeall(cpart)

    return res


def compute_homsym_with_alphabet_symmetrica(n, length, alphabet='x'):
    """
    Compute the expansion of a homogeneous(=complete) symmetric
    function labeled by a INTEGER number as a POLYNOM erg.
    The object number may also be a  PARTITION or a HOM_SYM object.
    The INTEGER laenge specifies the length of the alphabet.
    Both routines are the same.

    EXAMPLES::

        sage: symmetrica.compute_homsym_with_alphabet(3,1,'x')
        x^3
        sage: symmetrica.compute_homsym_with_alphabet([2,1],1,'x')
        x^3
        sage: symmetrica.compute_homsym_with_alphabet([2,1],2,'x')
        x0^3 + 2*x0^2*x1 + 2*x0*x1^2 + x1^3
        sage: symmetrica.compute_homsym_with_alphabet([2,1],2,'a,b')
        a^3 + 2*a^2*b + 2*a*b^2 + b^3
        sage: symmetrica.compute_homsym_with_alphabet([2,1],2,'x').parent()
        Multivariate Polynomial Ring in x0, x1 over Integer Ring
    """
    late_import()
    cdef OP cn = callocobject(), clength = callocobject(), cresult = callocobject()

    if isinstance(n, (int, Integer)):
        _op_integer(n, cn)
    elif isinstance(n, (builtinlist, Partition)):
        _op_partition(n, cn)
    else:
        raise NotImplementedError("n must be an integer or a partition")

    _op_integer(length, clength)

    sig_on()
    compute_homsym_with_alphabet(cn, clength, cresult)
    sig_off()

    res = _py_polynom_alphabet(cresult, alphabet, length)

    freeall(cresult)
    freeall(cn)
    freeall(clength)

    return res


def compute_elmsym_with_alphabet_symmetrica(n, length, alphabet='x'):
    """
    Compute the expansion of a elementary symmetric
    function labeled by a INTEGER number as a POLYNOM erg.
    The object number may also be a  PARTITION or a ELM_SYM object.
    The INTEGER length specifies the length of the alphabet.
    Both routines are the same.

    EXAMPLES::

        sage: a = symmetrica.compute_elmsym_with_alphabet(2,2); a
        x0*x1
        sage: a.parent()
        Multivariate Polynomial Ring in x0, x1 over Integer Ring
        sage: a = symmetrica.compute_elmsym_with_alphabet([2],2); a
        x0*x1
        sage: symmetrica.compute_elmsym_with_alphabet(3,2)
        0
        sage: symmetrica.compute_elmsym_with_alphabet([3,2,1],2)
        0
    """
    late_import()
    cdef OP cn = callocobject(), clength = callocobject(), cresult = callocobject()

    if isinstance(n, (int, Integer)):
        if n > length:
            return 0
        _op_integer(n, cn)
    elif isinstance(n, (builtinlist, Partition)):
        if max(n) > length:
            return 0
        _op_partition(n, cn)
    else:
        raise NotImplementedError("n must be an integer or a partition")

    _op_integer(length, clength)

    sig_on()
    compute_elmsym_with_alphabet(cn, clength, cresult)
    sig_off()

    res = _py_polynom_alphabet(cresult, alphabet, length)

    freeall(cresult)
    freeall(cn)
    freeall(clength)

    return res


def compute_monomial_with_alphabet_symmetrica(n, length, alphabet='x'):
    """
    Compute the expansion of a monomial symmetric
    function labeled by a PARTITION number as a POLYNOM erg.
    The INTEGER laenge specifies the length of the alphabet.

    EXAMPLES::

        sage: symmetrica.compute_monomial_with_alphabet([2,1],2,'x')
        x0^2*x1 + x0*x1^2
        sage: symmetrica.compute_monomial_with_alphabet([1,1,1],2,'x')
        0
        sage: symmetrica.compute_monomial_with_alphabet(2,2,'x')
        x0^2 + x1^2
        sage: symmetrica.compute_monomial_with_alphabet(2,2,'a,b')
        a^2 + b^2
        sage: symmetrica.compute_monomial_with_alphabet(2,2,'x').parent()
        Multivariate Polynomial Ring in x0, x1 over Integer Ring
    """
    late_import()
    cdef OP cn = callocobject(), clength = callocobject(), cresult = callocobject()
    if isinstance(n, (int, Integer)):
        _op_integer(n, cn)
    elif isinstance(n, (builtinlist, Partition)):
        if len(n) > length:
            return 0
        _op_partition(n, cn)
    else:
        raise NotImplementedError("n must be an integer or a partition")

    _op_integer(length, clength)

    sig_on()
    compute_monomial_with_alphabet(cn, clength, cresult)
    sig_off()

    res = _py_polynom_alphabet(cresult, alphabet, length)

    freeall(cresult)
    freeall(cn)
    freeall(clength)

    return res


def compute_powsym_with_alphabet_symmetrica(n, length, alphabet='x'):
    """
    Compute the expansion of a power symmetric
    function labeled by a INTEGER label or by a PARTITION label
    or a POW_SYM label as a POLYNOM erg.
    The INTEGER laenge specifies the length of the alphabet.

    EXAMPLES::

        sage: symmetrica.compute_powsym_with_alphabet(2,2,'x')
        x0^2 + x1^2
        sage: symmetrica.compute_powsym_with_alphabet(2,2,'x').parent()
        Multivariate Polynomial Ring in x0, x1 over Integer Ring
        sage: symmetrica.compute_powsym_with_alphabet([2],2,'x')
        x0^2 + x1^2
        sage: symmetrica.compute_powsym_with_alphabet([2],2,'a,b')
        a^2 + b^2
        sage: symmetrica.compute_powsym_with_alphabet([2,1],2,'a,b')
        a^3 + a^2*b + a*b^2 + b^3
    """
    late_import()
    cdef OP cn = callocobject(), clength = callocobject(), cresult = callocobject()

    if isinstance(n, (int, Integer)):
        _op_integer(n, cn)
    elif isinstance(n, (builtinlist, Partition)):
        _op_partition(n, cn)
    else:
        raise NotImplementedError("need to write code for POW_SYM")

    _op_integer(length, clength)

    sig_on()
    compute_powsym_with_alphabet(cn, clength, cresult)
    sig_off()

    res = _py_polynom_alphabet(cresult, alphabet, length)

    freeall(cresult)
    freeall(cn)
    freeall(clength)

    return res


def compute_schur_with_alphabet_det_symmetrica(part, length, alphabet='x'):
    """
    EXAMPLES::

        sage: symmetrica.compute_schur_with_alphabet_det(2,2)
        x0^2 + x0*x1 + x1^2
        sage: symmetrica.compute_schur_with_alphabet_det([2],2)
        x0^2 + x0*x1 + x1^2
        sage: symmetrica.compute_schur_with_alphabet_det(Partition([2]),2)
        x0^2 + x0*x1 + x1^2
        sage: symmetrica.compute_schur_with_alphabet_det(Partition([2]),2,'y')
        y0^2 + y0*y1 + y1^2
        sage: symmetrica.compute_schur_with_alphabet_det(Partition([2]),2,'a,b')
        a^2 + a*b + b^2
    """
    cdef OP cpart = callocobject(), cresult = callocobject(), clength = callocobject()

    if isinstance(part, (int, Integer)):
        _op_partition([part], cpart)
    elif isinstance(part, (builtinlist, Partition)):
        _op_partition(part, cpart)
    else:
        raise NotImplementedError("n must be an integer or partition")

    _op_integer(length, clength)

    sig_on()
    compute_schur_with_alphabet_det(cpart, clength, cresult)
    sig_off()

    res = _py_polynom_alphabet(cresult, alphabet, length)

    freeall(cresult)
    freeall(cpart)
    freeall(clength)

    return res


def part_part_skewschur_symmetrica(outer, inner):
    """
    Return the skew Schur function s_{outer/inner}.

    EXAMPLES::

        sage: symmetrica.part_part_skewschur([3,2,1],[2,1])
        s[1, 1, 1] + 2*s[2, 1] + s[3]
    """
    cdef OP couter = callocobject(), cinner = callocobject(), cresult = callocobject()

    _op_partition(outer, couter)
    _op_partition(inner, cinner)

    sig_on()
    part_part_skewschur(couter, cinner, cresult)
    sig_off()

    res = _py(cresult)

    freeall(couter)
    freeall(cinner)
    freeall(cresult)

    return res


def hall_littlewood_symmetrica(part):
    """
    Compute the so called Hall Littlewood Polynomials, i.e.
    a SCHUR object, whose coefficient are polynomials in one
    variable. The method, which is used for the computation is described
    in the paper: A.O. Morris The Characters of the group GL(n,q)
    Math Zeitschr 81, 112-123 (1963)
    """

    cdef OP cpart = callocobject(), cresult = callocobject()
    cdef OP pointer

    if len(part) == 0:
        raise TypeError("part must be a partition of a positive integer")

    _op_partition(part, cpart)

    sig_on()
    hall_littlewood(cpart, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(cpart)

    return res


def t_SCHUR_MONOMIAL_symmetrica(schur):
    cdef OP cschur = callocobject(), cresult = callocobject()

    _op_schur(schur, cschur)

    sig_on()
    t_SCHUR_MONOMIAL(cschur, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(cschur)

    return res


def t_SCHUR_HOMSYM_symmetrica(schur):
    cdef OP cschur = callocobject(), cresult = callocobject()

    _op_schur(schur, cschur)

    sig_on()
    t_SCHUR_HOMSYM(cschur, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(cschur)

    return res


def t_SCHUR_ELMSYM_symmetrica(schur):
    cdef OP cschur = callocobject(), cresult = callocobject()

    _op_schur(schur, cschur)

    sig_on()
    t_SCHUR_ELMSYM(cschur, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(cschur)

    return res


def t_SCHUR_POWSYM_symmetrica(schur):
    cdef OP cschur = callocobject(), cresult = callocobject()

    _op_schur(schur, cschur)

    sig_on()
    t_SCHUR_POWSYM(cschur, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(cschur)

    return res


def t_POLYNOM_SCHUR_symmetrica(p):
    r"""
    Convert a symmetric polynomial with base ring `\QQ` or `\ZZ` into a
    symmetric function in the Schur basis.
    """
    cdef OP polynom = callocobject(), cresult = callocobject()

    _op_polynom(p, polynom)

    if not symmetricp(polynom):
        raise ValueError("the polynomial must be symmetric")

    sig_on()
    t_POLYNOM_SCHUR(polynom, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(polynom)

    return res


def t_MONOMIAL_HOMSYM_symmetrica(monomial):
    cdef OP cmonomial = callocobject(), cresult = callocobject()

    _op_monomial(monomial, cmonomial)

    sig_on()
    t_MONOMIAL_HOMSYM(cmonomial, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(cmonomial)

    return res


def t_MONOMIAL_ELMSYM_symmetrica(monomial):
    cdef OP cmonomial = callocobject(), cresult = callocobject()

    _op_monomial(monomial, cmonomial)

    sig_on()
    t_MONOMIAL_ELMSYM(cmonomial, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(cmonomial)

    return res


def t_MONOMIAL_SCHUR_symmetrica(monomial):
    cdef OP cmonomial = callocobject(), cresult = callocobject()

    _op_monomial(monomial, cmonomial)

    sig_on()
    t_MONOMIAL_SCHUR(cmonomial, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(cmonomial)

    return res


def t_MONOMIAL_POWSYM_symmetrica(monomial):
    cdef OP cmonomial = callocobject(), cresult = callocobject()

    _op_monomial(monomial, cmonomial)

    sig_on()
    t_MONOMIAL_POWSYM(cmonomial, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(cmonomial)

    return res


def t_POLYNOM_MONOMIAL_symmetrica(p):
    r"""
    Convert a symmetric polynomial with base ring `\QQ` or `\ZZ` into a
    symmetric function in the monomial basis.
    """
    cdef OP polynom = callocobject(), cresult = callocobject()

    _op_polynom(p, polynom)

    if not symmetricp(polynom):
        raise ValueError("the polynomial must be symmetric")

    sig_on()
    t_POLYNOM_MONOMIAL(polynom, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(polynom)

    return res


def t_ELMSYM_SCHUR_symmetrica(elmsym):
    cdef OP celmsym = callocobject(), cresult = callocobject()

    _op_elmsym(elmsym, celmsym)

    sig_on()
    t_ELMSYM_SCHUR(celmsym, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(celmsym)

    return res


def t_ELMSYM_POWSYM_symmetrica(elmsym):
    cdef OP celmsym = callocobject(), cresult = callocobject()

    _op_elmsym(elmsym, celmsym)

    sig_on()
    t_ELMSYM_POWSYM(celmsym, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(celmsym)

    return res


def t_ELMSYM_MONOMIAL_symmetrica(elmsym):
    cdef OP celmsym = callocobject(), cresult = callocobject()

    _op_elmsym(elmsym, celmsym)

    sig_on()
    t_ELMSYM_MONOMIAL(celmsym, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(celmsym)

    return res


def t_ELMSYM_HOMSYM_symmetrica(elmsym):
    cdef OP celmsym = callocobject(), cresult = callocobject()

    _op_elmsym(elmsym, celmsym)

    sig_on()
    t_ELMSYM_HOMSYM(celmsym, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(celmsym)

    return res


def t_POLYNOM_ELMSYM_symmetrica(p):
    r"""
    Convert a symmetric polynomial with base ring `\QQ` or `\ZZ` into a
    symmetric function in the elementary basis.
    """
    cdef OP polynom = callocobject(), cresult = callocobject()

    _op_polynom(p, polynom)

    if not symmetricp(polynom):
        raise ValueError("the polynomial must be symmetric")

    sig_on()
    t_POLYNOM_ELMSYM(polynom, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(polynom)

    return res


def t_HOMSYM_SCHUR_symmetrica(homsym):
    cdef OP chomsym = callocobject(), cresult = callocobject()

    _op_homsym(homsym, chomsym)

    sig_on()
    t_HOMSYM_SCHUR(chomsym, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(chomsym)

    return res


def t_HOMSYM_POWSYM_symmetrica(homsym):
    cdef OP chomsym = callocobject(), cresult = callocobject()

    _op_homsym(homsym, chomsym)

    sig_on()
    t_HOMSYM_POWSYM(chomsym, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(chomsym)

    return res


def t_HOMSYM_MONOMIAL_symmetrica(homsym):
    cdef OP chomsym = callocobject(), cresult = callocobject()

    _op_homsym(homsym, chomsym)

    sig_on()
    t_HOMSYM_MONOMIAL(chomsym, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(chomsym)

    return res


def t_HOMSYM_ELMSYM_symmetrica(homsym):
    cdef OP chomsym = callocobject(), cresult = callocobject()

    _op_homsym(homsym, chomsym)

    sig_on()
    t_HOMSYM_ELMSYM(chomsym, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(chomsym)

    return res


def t_POWSYM_MONOMIAL_symmetrica(powsym):
    cdef OP cpowsym = callocobject(), cresult = callocobject()

    _op_powsym(powsym, cpowsym)

    sig_on()
    t_POWSYM_MONOMIAL(cpowsym, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(cpowsym)

    return res


def t_POWSYM_SCHUR_symmetrica(powsym):
    cdef OP cpowsym = callocobject(), cresult = callocobject()

    _op_powsym(powsym, cpowsym)

    sig_on()
    t_POWSYM_SCHUR(cpowsym, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(cpowsym)

    return res


def t_POWSYM_ELMSYM_symmetrica(powsym):
    cdef OP cpowsym = callocobject(), cresult = callocobject()

    _op_powsym(powsym, cpowsym)

    sig_on()
    t_POWSYM_ELMSYM(cpowsym, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(cpowsym)

    return res


def t_POWSYM_HOMSYM_symmetrica(powsym):
    cdef OP cpowsym = callocobject(), cresult = callocobject()

    _op_powsym(powsym, cpowsym)

    sig_on()
    t_POWSYM_HOMSYM(cpowsym, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(cpowsym)

    return res


def t_POLYNOM_POWER_symmetrica(p):
    r"""
    Convert a symmetric polynomial with base ring `\QQ` or `\ZZ` into a
    symmetric function in the power sum basis.
    """
    cdef OP polynom = callocobject(), cresult = callocobject()

    _op_polynom(p, polynom)

    if not symmetricp(polynom):
        raise ValueError("the polynomial must be symmetric")

    sig_on()
    t_POLYNOM_POWER(polynom, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cresult)
    freeall(polynom)

    return res


def mult_schur_schur_symmetrica(s1, s2):
    cdef OP cs1 = callocobject(), cs2 = callocobject(), cresult = callocobject()

    _op_schur(s1, cs1)
    _op_schur(s2, cs2)

    sig_on()
    mult_schur_schur(cs1, cs2, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cs1)
    freeall(cs2)
    freeall(cresult)

    return res


def mult_monomial_monomial_symmetrica(m1, m2):
    cdef OP cm1 = callocobject(), cm2 = callocobject(), cresult = callocobject()

    _op_monomial(m1, cm1)
    _op_monomial(m2, cm2)

    sig_on()
    mult_monomial_monomial(cm1, cm2, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cm1)
    freeall(cm2)
    freeall(cresult)

    return res


def scalarproduct_schur_symmetrica(s1, s2):
    cdef OP cs1 = callocobject(), cs2 = callocobject(), cresult = callocobject()

    _op_schur(s1, cs1)
    _op_schur(s2, cs2)

    sig_on()
    scalarproduct_schur(cs1, cs2, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cs1)
    freeall(cs2)
    freeall(cresult)

    return res
