# sage_setup: distribution = sagemath-modules
from sage.libs.mpfr.types cimport mpfr_t

cdef extern from "mpc.h":
    ctypedef struct __mpc_struct:
        mpfr_t re
        mpfr_t im
    ctypedef __mpc_struct mpc_t[1]
    ctypedef __mpc_struct* mpc_ptr
    ctypedef __mpc_struct* mpc_srcptr

    ctypedef enum mpc_rnd_t:
        MPC_RNDNN
        MPC_RNDZN
        MPC_RNDUN
        MPC_RNDDN
        MPC_RNDNZ
        MPC_RNDZZ
        MPC_RNDUZ
        MPC_RNDDZ
        MPC_RNDNU
        MPC_RNDZU
        MPC_RNDUU
        MPC_RNDDU
        MPC_RNDND
        MPC_RNDZD
        MPC_RNDUD
        MPC_RNDDD
