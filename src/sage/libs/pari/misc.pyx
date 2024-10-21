# sage_setup: distribution = sagemath-pari

from cysignals.signals cimport sig_on
from cypari2.paridecl cimport *
from cypari2.stack cimport new_gen


cdef Gen new_t_POL_from_int_star(int* vals, unsigned long length, long varnum):
    """
    Convert an array of ints to a PARI polynomial.

    Note that ``length = degree + 1``, so that recognizing 0 is easier.
    """
    cdef GEN z
    cdef unsigned long i

    sig_on()
    z = cgetg(length + 2, t_POL)
    if length == 0:
        # Polynomial is zero
        z[1] = evalvarn(varnum) + evalsigne(0)
    else:
        z[1] = evalvarn(varnum) + evalsigne(1)
        for i in range(length):
            set_gel(z, i+2, stoi(vals[i]))

    return new_gen(z)
