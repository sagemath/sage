cdef extern from 'symmetrica/def.h':
    INT strict_to_odd_part(OP s, OP o)
    INT odd_to_strict_part(OP o, OP s)
    INT q_core(OP part, OP d,  OP core)
    INT gupta_nm(OP n, OP m, OP res)
    INT gupta_tafel(OP max, OP res)
    INT random_partition(OP nx, OP res)


def strict_to_odd_part_symmetrica(part):
    """
    Implement the bijection between strict partitions
    and partitions with odd parts. input is a VECTOR type partition, the
    result is a partition of the same weight with only odd parts.
    """

    #Make sure that the partition is strict
    cdef INT i
    for i from 0 <= i < len(part)-1:
        if part[i] == part[i+1]:
            raise ValueError("the partition part (= %s) must be strict" % str(part))

    cdef OP cpart, cres
    anfang()
    cpart = callocobject()
    cres = callocobject()

    _op_partition(part, cpart)

    strict_to_odd_part(cpart, cres)

    res = _py(cres)

    freeall(cpart)
    freeall(cres)
    ende()

    return res


def odd_to_strict_part_symmetrica(part):
    """
    Implement the bijection between partitions with odd parts
    and strict partitions. input is a VECTOR type partition, the
    result is a partition of the same weight with different parts.
    """

    #Make sure that the partition is strict
    cdef INT i
    for i from 0 <= i < len(part):
        if part[i] % 2 == 0:
            raise ValueError("the partition part (= %s) must be odd" % str(part))

    cdef OP cpart, cres
    anfang()
    cpart = callocobject()
    cres = callocobject()

    _op_partition(part, cpart)

    odd_to_strict_part(cpart, cres)

    res = _py(cres)

    freeall(cpart)
    freeall(cres)
    ende()

    return res


def q_core_symmetrica(part, d):
    """
    Compute the q-core of a PARTITION object
    part. This is the remaining partition (=res) after
    removing of all hooks of length d (= INTEGER object).
    The result may be an empty object, if the whole
    partition disappears.
    """


    cdef OP cpart, cres, cd
    anfang()
    cpart = callocobject()
    cd = callocobject()
    cres = callocobject()

    _op_partition(part, cpart)
    _op_integer(d, cd)

    q_core(cpart, cd, cres)

    res = _py(cres)

    freeall(cpart)
    freeall(cres)
    freeall(cd)
    ende()

    return res


def gupta_nm_symmetrica(n, m):
    """
    This routine computes the number of partitions
    of n with maximal part m. The result is erg. The
    input n,m must be INTEGER objects. The result is
    freed first to an empty object. The result must
    be a different from m and n.
    """


    cdef OP cn, cm, cres
    anfang()
    cm = callocobject()
    cn = callocobject()
    cres = callocobject()


    _op_integer(n, cn)
    _op_integer(m, cm)

    gupta_nm(cn, cm, cres)

    res = _py(cres)

    freeall(cn)
    freeall(cres)
    freeall(cm)
    ende()

    return res


def gupta_tafel_symmetrica(max):
    """
    It computes the table of the above values. The entry
    n,m is the result of gupta_nm. mat is freed first.
    max must be an INTEGER object, it is the maximum
    weight for the partitions. max must be different from
    result.
    """


    cdef OP cmax, cres
    anfang()

    cmax = callocobject()
    cres = callocobject()


    _op_integer(max, cmax)

    gupta_tafel(cmax, cres)

    res = _py(cres)

    freeall(cmax)
    freeall(cres)
    ende()

    return res


def random_partition_symmetrica(n):
    """
    Return a random partition p of the entered weight w.

    w must be an INTEGER object, p becomes a PARTITION object.
    Type of partition is VECTOR . It uses the algorithm of
    Nijenhuis and Wilf, p.76
    """


    cdef OP cn, cres
    anfang()

    cn = callocobject()
    cres = callocobject()


    _op_integer(n, cn)

    random_partition(cn, cres)

    res = _py(cres)

    freeall(cn)
    freeall(cres)
    ende()

    return res
