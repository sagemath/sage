# distutils: libraries = gmp pari
"""
Declarations for private functions from PARI

Ideally, we shouldn't use these, but for technical reasons, we have to.
"""

from .types cimport *

cdef extern from "pari/paripriv.h":
    int t_FF_FpXQ, t_FF_Flxq, t_FF_F2xq

    int gpd_QUIET, gpd_TEST, gpd_EMACS, gpd_TEXMACS

    struct pariout_t:
        char format  # e,f,g
        long fieldw  # 0 (ignored) or field width
        long sigd    # -1 (all) or number of significant digits printed
        int sp       # 0 = suppress whitespace from output
        int prettyp  # output style: raw, prettyprint, etc
        int TeXstyle

    struct gp_data:
        pariout_t *fmt
        unsigned long flags

    extern gp_data* GP_DATA

cdef extern from "pari/anal.h":
    char* closure_func_err()
