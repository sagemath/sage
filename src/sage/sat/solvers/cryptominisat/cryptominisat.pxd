from sage.sat.solvers.satsolver cimport SatSolver
from .decl cimport Solver

cdef class CryptoMiniSat(SatSolver):
    cdef Solver *_solver
