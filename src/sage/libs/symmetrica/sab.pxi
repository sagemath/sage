# sage_setup: distribution = sagemath-combinat
cdef extern from 'symmetrica/def.h':
    INT dimension_symmetrization(OP n, OP part, OP a)
    INT bdg(OP part, OP perm, OP D)
    INT sdg(OP part, OP perm, OP D)
    INT odg(OP part, OP perm, OP D)
    INT ndg(OP part, OP perm, OP D)
    INT specht_dg(OP part, OP perm, OP D)
    INT glmndg(OP m, OP n, OP M, INT VAR)


def dimension_symmetrization_symmetrica(n, part):
    """
    Compute the dimension of the degree of a irreducible
    representation of the GL_n, n is a INTEGER object, labeled
    by the PARTITION object a.
    """
    cdef OP cn, cpart, cres


    cn    = callocobject()
    cpart = callocobject()
    cres  = callocobject()

    _op_partition(part, cpart)
    _op_integer(n, cn)

    dimension_symmetrization(cn, cpart, cres)
    res = _py(cres)

    freeall(cn)
    freeall(cpart)
    freeall(cres)


    return res


def bdg_symmetrica(part, perm):
    """
    Calculate the irreducible matrix representation
    D^part(perm), whose entries are of integral numbers.

    REFERENCE: H. Boerner:
               Darstellungen von Gruppen, Springer 1955.
               pp. 104-107.
    """
    cdef OP cpart, cperm, cD


    cpart = callocobject()
    cperm = callocobject()
    cD    = callocobject()

    _op_partition(part, cpart)
    _op_permutation(perm, cperm)

    bdg(cpart, cperm, cD)
    res = _py_matrix(cD)

    freeall(cpart)
    freeall(cperm)
    freeall(cD)


def sdg_symmetrica(part, perm):
    """
    Calculate the irreducible matrix representation
    D^part(perm), which consists of rational numbers.

    REFERENCE: G. James/ A. Kerber:
               Representation Theory of the Symmetric Group.
               Addison/Wesley 1981.
               pp. 124-126.
    """
    cdef OP cpart, cperm, cD


    cpart = callocobject()
    cperm = callocobject()
    cD    = callocobject()

    _op_partition(part, cpart)
    _op_permutation(perm, cperm)

    sdg(cpart, cperm, cD)
    res = _py_matrix(cD)

    freeall(cpart)
    freeall(cperm)
    freeall(cD)

    return res


def odg_symmetrica(part, perm):
    """
    Calculate the irreducible matrix representation
    D^part(perm), which consists of real numbers.

    REFERENCE: G. James/ A. Kerber:
               Representation Theory of the Symmetric Group.
               Addison/Wesley 1981.
               pp. 127-129.
    """
    cdef OP cpart, cperm, cD


    cpart = callocobject()
    cperm = callocobject()
    cD    = callocobject()

    _op_partition(part, cpart)
    _op_permutation(perm, cperm)

    odg(cpart, cperm, cD)
    res = _py_matrix(cD)

    freeall(cpart)
    freeall(cperm)
    freeall(cD)

    return res


def ndg_symmetrica(part, perm):
    cdef OP cpart, cperm, cD

    cpart = callocobject()
    cperm = callocobject()
    cD    = callocobject()

    _op_partition(part, cpart)
    _op_permutation(perm, cperm)

    ndg(cpart, cperm, cD)
    res = _py_matrix(cD)

    freeall(cpart)
    freeall(cperm)
    freeall(cD)

    return res


def specht_dg_symmetrica(part, perm):
    cdef OP cpart, cperm, cD

    cpart = callocobject()
    cperm = callocobject()
    cD    = callocobject()

    _op_partition(part, cpart)
    _op_permutation(perm, cperm)

    specht_dg(cpart, cperm, cD)
    res = _py_matrix(cD)

    freeall(cpart)
    freeall(cperm)
    freeall(cD)

    return res


## def glmndg_symmetrica(m, n, VAR=0):
##     """
##     If VAR is equal to 0 the orthogonal representation
##     is used for the decomposition, otherwise, if VAR
##     equals 1, the natural representation is considered.

##     The result is the  VECTOR-Object M, consisting of
##     components of type MATRIX, representing the several
##     irreducible matrix representations of GLm(C) with
##     part_1' <= m, where part is a partition of n.

##     """
##     cdef OP cm, cn, cM

##     cm = callocobject()
##     _op_integer(m, cm)

##     cn = callocobject()
##     _op_integer(n, cn)

##     cM = callocobject()

##     glmndg(cm, cn, cM, VAR)
##     res = _py(cM)


##    freeall(cm)
##    freeall(cn)
##    freeall(cM)

##    return res
