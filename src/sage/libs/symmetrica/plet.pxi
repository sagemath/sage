cdef extern from 'symmetrica/def.h':
    INT plethysm(OP s1, OP s2, OP res)
    INT schur_schur_plet(OP p1, OP p2, OP res)


def plethysm_symmetrica(outer, inner):
    cdef OP couter = callocobject(), cinner = callocobject(), cresult = callocobject()

    _op_schur(outer, couter)
    _op_schur(inner, cinner)

    sig_on()
    plethysm(couter, cinner, cresult)
    sig_off()

    res = _py(cresult)

    freeall(couter)
    freeall(cinner)
    freeall(cresult)

    return res


def schur_schur_plet_symmetrica(outer, inner):
    cdef OP couter = callocobject(), cinner = callocobject(), cresult = callocobject()

    _op_partition(outer, couter)
    _op_partition(inner, cinner)

    sig_on()
    schur_schur_plet(couter, cinner, cresult)
    sig_off()

    res = _py(cresult)

    freeall(couter)
    freeall(cinner)
    freeall(cresult)

    return res
