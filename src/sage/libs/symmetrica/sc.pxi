cdef extern from 'symmetrica/def.h':
    INT chartafel(OP degree, OP result)
    INT charvalue(OP irred, OP cls, OP result, OP table)
    INT kranztafel(OP a, OP b, OP res, OP co, OP cl)
    INT c_ijk_sn(OP i, OP j, OP k, OP res)


def chartafel_symmetrica(n):
    """
    You enter the degree of the symmetric group, as INTEGER
    object and the result is a MATRIX object: the charactertable
    of the symmetric group of the given degree.

    EXAMPLES::

        sage: symmetrica.chartafel(3)
        [ 1  1  1]
        [-1  0  2]
        [ 1 -1  1]
        sage: symmetrica.chartafel(4)
        [ 1  1  1  1  1]
        [-1  0 -1  1  3]
        [ 0 -1  2  0  2]
        [ 1  0 -1 -1  3]
        [-1  1  1 -1  1]
    """
    cdef OP cn, cres

    cn   = callocobject()
    cres = callocobject()

    _op_integer(n, cn)

    chartafel(cn, cres)

    res = _py(cres)

    freeall(cn)
    freeall(cres)

    return res


def charvalue_symmetrica(irred, cls, table=None):
    """
    You enter a PARTITION object part, labelling the irreducible
    character, you enter a PARTITION object class, labeling the class
    or class may be a PERMUTATION object, then result becomes the value
    of that character on that class or permutation. Note that the
    table may be NULL, in which case the value is computed, or it may be
    taken from a precalculated charactertable.

    FIXME: add table parameter

    EXAMPLES::

        sage: n = 3
        sage: m = matrix([[symmetrica.charvalue(irred, cls) for cls in Partitions(n)] for irred in Partitions(n)]); m
        [ 1  1  1]
        [-1  0  2]
        [ 1 -1  1]
        sage: m == symmetrica.chartafel(n)
        True
        sage: n = 4
        sage: m = matrix([[symmetrica.charvalue(irred, cls) for cls in Partitions(n)] for irred in Partitions(n)])
        sage: m == symmetrica.chartafel(n)
        True
    """
    cdef OP cirred, cclass, ctable, cresult

    cirred = callocobject()
    cclass = callocobject()
    cresult = callocobject()

    if table is None:
        ctable = NULL
    else:
        ctable = callocobject()
        _op_matrix(table, ctable)

    #FIXME: assume that class is a partition
    _op_partition(cls, cclass)

    _op_partition(irred, cirred)

    charvalue(cirred, cclass, cresult, ctable)

    res = _py(cresult)

    freeall(cirred)
    freeall(cclass)
    freeall(cresult)
    if ctable != NULL:
        freeall(ctable)

    return res


def kranztafel_symmetrica(a, b):
    r"""
    You enter the INTEGER objects, say `a` and `b`, and ``res`` becomes a
    MATRIX object, the charactertable of `S_b \wr S_a`, ``co`` becomes a
    VECTOR object of classorders and ``cl`` becomes a VECTOR object of
    the classlabels.

    EXAMPLES::

        sage: (a,b,c) = symmetrica.kranztafel(2,2)
        sage: a
        [ 1 -1  1 -1  1]
        [ 1  1  1  1  1]
        [-1  1  1 -1  1]
        [ 0  0  2  0 -2]
        [-1 -1  1  1  1]
        sage: b
        [2, 2, 1, 2, 1]
        sage: for m in c: print(m)
        [0 0]
        [0 1]
        [0 0]
        [1 0]
        [0 2]
        [0 0]
        [1 1]
        [0 0]
        [2 0]
        [0 0]
    """
    cdef OP ca, cb, cres, cco, ccl

    ca = callocobject()
    cb = callocobject()
    cres = callocobject()
    cco = callocobject()
    ccl = callocobject()

    _op_integer(a, ca)
    _op_integer(b, cb)

    kranztafel(ca,cb,cres,cco,ccl)

    res = _py(cres)
    co  = _py(cco)
    cl  = _py(ccl)

    freeall(ca)
    freeall(cb)
    freeall(cres)
    freeall(cco)
    freeall(ccl)

    return (res, co, cl)


## def c_ijk_sn_symmetrica(i, j, k):
##     """
##     computes the coefficients of the class multiplication in the
##     group algebra of the S_n. It uses the method described in
##     Curtis/Reiner: Methods of representation theory I p. 216

##     EXAMPLES:

##     """

##     cdef OP ci, cj, ck, cresult


##     ci = callocobject()
##     cj = callocobject()
##     cresult = callocobject()
##     ck = callocobject()


##     _op_partition(i, ci)
##     _op_partition(j, cj)
##     _op_partition(k, ck)

##     c_ijk_sn(ci, cj, ck, cresult)

##     res = _py(cresult)

##     freeall(ci)
##     freeall(cj)
##     freeall(cresult)
##     freeall(ck)

##     return res
