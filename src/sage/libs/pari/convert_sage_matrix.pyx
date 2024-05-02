# sage.doctest: needs sage.modules

from cypari2.gen cimport Gen
from cypari2.types cimport (GEN, typ, t_INT, t_FRAC, t_REAL, t_COMPLEX,
                            t_INTMOD, t_PADIC, t_INFINITY, t_VEC, t_COL,
                            t_VECSMALL, t_MAT, t_STR,
                            lg, precp)
from cypari2.paridecl cimport *

from sage.matrix.args cimport (MatrixArgs, MA_ENTRIES_SEQ_SEQ,
                               MA_ENTRIES_SEQ_FLAT, MA_ENTRIES_CALLABLE,
                               MA_ENTRIES_UNKNOWN, MA_ENTRIES_SCALAR)

from .convert_sage cimport gen_to_sage


def gen_to_sage_matrix(Gen z, locals=None):
    cdef GEN g = z.g
    nc = lg(g) - 1
    nr = 0 if nc == 0 else lg(gel(g,1)) - 1
    ma = MatrixArgs.__new__(MatrixArgs)
    ma.nrows = nr
    ma.ncols = nc
    ma.entries = [gen_to_sage(z[i,j], locals) for i in range(nr) for j in range(nc)]
    return ma.matrix()


def pari_typ_to_entries_type(MatrixArgs self):
    """
    Determine the ``entries_type`` of a :class:`sage.matrix.args.MatrixArgs`
    with PARI entries.

    This will modify the entries.

    TESTS:

    ``MA_ENTRIES_SEQ_SEQ``::

        sage: from sage.libs.pari.convert_sage import pari_typ_to_entries_type
        sage: from sage.matrix.args import MatrixArgs
        sage: ma = MatrixArgs(QQ, entries=pari("[1,2;3,4]"))
        sage: 0x10_03 == pari_typ_to_entries_type(ma)
        True

    ``MA_ENTRIES_SEQ_FLAT``::

        sage: ma = MatrixArgs(QQ, entries=pari("[1,2]"))
        sage: 0x10_04 == pari_typ_to_entries_type(ma)
        True
        sage: ma = MatrixArgs(QQ, entries=pari(vector([1,2])))
        sage: 0x10_04 == pari_typ_to_entries_type(ma)
        True
        sage: ma = MatrixArgs(QQ, entries=pari(matrix(2, range(4))[0]))
        sage: 0x10_04 == pari_typ_to_entries_type(ma)
        True

    ``MA_ENTRIES_CALLABLE``::

        sage: ma = MatrixArgs(QQ, entries=pari(lambda x: x))
        sage: 0x13_06 == pari_typ_to_entries_type(ma)
        True

    ``MA_ENTRIES_SCALAR``::

        sage: ma = MatrixArgs(QQ, entries=pari(1/2))
        sage: 0x17_02 == pari_typ_to_entries_type(ma)
        True

    ``MA_ENTRIES_UNKNOWN``::

        sage: ma = MatrixArgs(QQ, entries=pari('"2"'))
        sage: 0 == pari_typ_to_entries_type(ma)
        True

    A second call gives an error::

        sage: ma = MatrixArgs(QQ, entries=pari("[1,2]"))
        sage: 0x10_04 == pari_typ_to_entries_type(ma)
        True
        sage: 0x10_04 == pari_typ_to_entries_type(ma)
        Traceback (most recent call last):
        ...
        ValueError: entries are not a PARI generator
    """
    if not isinstance(self.entries, Gen):
        raise ValueError("entries are not a PARI generator")
    cdef long t = typ((<Gen>self.entries).g)
    if t == t_MAT:
        R = self.base
        if R is None:
            self.entries = self.entries.Col().sage()
        else:
            self.entries = [[R(x) for x in v]
                            for v in self.entries.mattranspose()]
        return MA_ENTRIES_SEQ_SEQ
    elif t in [t_VEC, t_COL, t_VECSMALL, t_LIST]:
        self.entries = self.entries.sage()
        return MA_ENTRIES_SEQ_FLAT
    elif t == t_CLOSURE:
        return MA_ENTRIES_CALLABLE
    elif t == t_STR:
        return MA_ENTRIES_UNKNOWN
    else:
        self.entries = self.entries.sage()
        return MA_ENTRIES_SCALAR
