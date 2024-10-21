# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules

from cpython.object cimport *
from cysignals.signals cimport sig_on, sig_off

cdef extern from 'symmetrica/def.h':
    INT kostka_number(OP shape, OP content, OP result)
    INT kostka_tab(OP shape, OP content, OP result)
    INT kostka_tafel(OP n, OP result)


def kostka_number_symmetrica(shape, content):
    """
    Compute the kostkanumber, i.e. the number of
    tableaux of given shape, which is a PARTITION object, and
    of given content, which also is a PARTITION object, or a VECTOR
    object with INTEGER entries. The
    result is an INTEGER object, which is freed to an empty
    object at the beginning. The shape could also be a
    SKEWPARTITION object, then we compute the number of
    skewtableaux of the given shape.

    EXAMPLES::

        sage: symmetrica.kostka_number([2,1],[1,1,1])
        2
        sage: symmetrica.kostka_number([1,1,1],[1,1,1])
        1
        sage: symmetrica.kostka_number([3],[1,1,1])
        1
    """
    cdef OP cshape = callocobject(), ccontent = callocobject(), result = callocobject()

    if isinstance(shape, <type>builtinlist):
        if isinstance(shape[0], <type>builtinlist):
            shape = SkewPartition(shape)
        else:
            shape = Partition(shape)


    if isinstance(shape, <type>SkewPartition):
        _op_skew_partition(shape, cshape)
    else:
        _op_partition(shape, cshape)

    _op_partition(content, ccontent)

    kostka_number(ccontent, cshape, result)

    res = _py(result)

    freeall(cshape)
    freeall(ccontent)
    freeall(result)


    return res


def kostka_tab_symmetrica(shape, content):
    """
    Compute the list of tableaux of given shape
    and content. shape is a PARTITION object or a
    SKEWPARTITION object and
    content is a PARTITION object or a VECTOR object with
    INTEGER entries, the result becomes a
    LIST object whose entries are the computed TABLEAUX
    object.

    EXAMPLES::

        sage: symmetrica.kostka_tab([3],[1,1,1])
        [[[1, 2, 3]]]
        sage: symmetrica.kostka_tab([2,1],[1,1,1])
        [[[1, 2], [3]], [[1, 3], [2]]]
        sage: symmetrica.kostka_tab([1,1,1],[1,1,1])
        [[[1], [2], [3]]]
        sage: symmetrica.kostka_tab([[2,2,1],[1,1]],[1,1,1])
        [[[None, 1], [None, 2], [3]],
         [[None, 1], [None, 3], [2]],
         [[None, 2], [None, 3], [1]]]
        sage: symmetrica.kostka_tab([[2,2],[1]],[1,1,1])
        [[[None, 1], [2, 3]], [[None, 2], [1, 3]]]
    """
    late_import()

    cdef OP cshape = callocobject(), ccontent = callocobject(), result = callocobject()
    cdef INT err

    if isinstance(shape, <type>builtinlist):
        if isinstance(shape[0], <type>builtinlist):
            shape = SkewPartition(shape)
        else:
            shape = Partition(shape)


    if isinstance(shape, <type>SkewPartition):
        _op_skew_partition(shape, cshape)
    else:
        _op_partition(shape, cshape)


    #Check to make sure the content is compatible with the shape.

    _op_il_vector(content, ccontent)

    err = kostka_tab(cshape, ccontent, result)

    res = _py(result)

    freeall(cshape)
    freeall(ccontent)
    freeall(result)


    return res


def kostka_tafel_symmetrica(n):
    """
    Return the table of Kostka numbers of weight `n`.

    EXAMPLES::

        sage: symmetrica.kostka_tafel(1)
        [1]

        sage: symmetrica.kostka_tafel(2)
        [1 0]
        [1 1]

        sage: symmetrica.kostka_tafel(3)
        [1 0 0]
        [1 1 0]
        [1 2 1]

        sage: symmetrica.kostka_tafel(4)
        [1 0 0 0 0]
        [1 1 0 0 0]
        [1 1 1 0 0]
        [1 2 1 1 0]
        [1 3 2 3 1]

        sage: symmetrica.kostka_tafel(5)
        [1 0 0 0 0 0 0]
        [1 1 0 0 0 0 0]
        [1 1 1 0 0 0 0]
        [1 2 1 1 0 0 0]
        [1 2 2 1 1 0 0]
        [1 3 3 3 2 1 0]
        [1 4 5 6 5 4 1]
    """

    cdef OP cn = callocobject(), cresult = callocobject()

    _op_integer(n, cn)

    sig_on()
    kostka_tafel(cn, cresult)
    sig_off()

    res = _py(cresult)

    freeall(cn)
    freeall(cresult)

    return res
