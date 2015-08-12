#*****************************************************************************
#       Copyright (C) 2012 Martin Albrecht <martinralbrecht@googlemail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.sat.solvers.cryptominisat.decl cimport SolverConf as SolverConfC

cdef extern from "sage/sat/solvers/cryptominisat/solverconf_helper.h":
    ctypedef enum sc_type:
        t_int
        t_float
        t_double
        t_Var
        t_bool
        t_uint32_t
        t_uint64_t

    ctypedef struct sc_entry "sc_entry":
        char *name
        sc_type type
        void *target
        char *doc

    size_t setup_map(sc_entry *entries, SolverConfC conf, Py_ssize_t n)

cdef class SolverConf:
    cdef SolverConfC *_conf
    cdef sc_entry[100] _map
    cdef Py_ssize_t _nopts
