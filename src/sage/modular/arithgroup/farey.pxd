#*****************************************************************************
#       Copyright (C) 2011 Hartmut Monien <monien@th.physik.uni-bonn.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cdef extern from 'gmpxx.h':
    cdef cppclass mpz_class:
        mpz_class()
        mpz_class(int i)
        mpz_class(mpz_t z)
        mpz_class(mpz_class)
        mpz_t get_mpz_t()
        mpz_class operator%(mpz_class, mpz_class)
    cdef cppclass mpq_class:
        mpq_class()
        mpz_t get_num_mpz_t()
        mpz_t get_den_mpz_t()

cdef extern from 'sl2z.hpp':
    cppclass cpp_SL2Z "SL2Z":
        mpz_class a, b, c, d
        cpp_SL2Z(int, int, int, int)
        cpp_SL2Z(mpz_class, mpz_class, mpz_class, mpz_class)
        mpz_class a()
        mpz_class b()
        mpz_class c()
        mpz_class d()

cdef extern from 'farey.hpp':
    cppclass is_element_Gamma0:
        is_element_Gamma0(int)
    cppclass is_element_Gamma1:
        is_element_Gamma1(int)
    cppclass is_element_Gamma:
        is_element_Gamma(int)
    cppclass is_element_GammaH:
        is_element_GammaH(int, object)
    cppclass cpp_farey "FareySymbol":
        cpp_farey()
        cpp_farey(object)
        cpp_farey(object, is_element_Gamma*)
        cpp_farey(object, is_element_Gamma0*)
        cpp_farey(object, is_element_Gamma1*)
        cpp_farey(object, is_element_GammaH*)
        size_t genus()
        size_t index()
        size_t level()
        size_t nu2()
        size_t nu3()
        object is_element(mpz_t, mpz_t, mpz_t, mpz_t)
        object word_problem(mpz_t, mpz_t, mpz_t, mpz_t)
        size_t get_cusp_class(mpz_t, mpz_t)
        object get_cusps()
        object get_cusp_widths()
        object get_transformation_to_cusp(mpz_t, mpz_t)
        object get_fractions()
        object get_coset()
        object get_generators()
        object get_pairings()
        object get_paired_sides()
        object get_pairing_matrices()
        object dumps()

