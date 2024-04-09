##############################################################################
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#       Copyright (C) 2016 Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.numerical.backends.generic_backend cimport GenericBackend

cdef class InteractiveLPBackend(GenericBackend):

    cdef object lp
    cdef object row_names
    cdef object prob_name

    cdef object lp_std_form
    cdef object std_form_transformation
    cdef object final_dictionary
    cdef int verbosity

    cpdef int add_variable(self,
                           lower_bound=*,
                           upper_bound=*,
                           binary=*,
                           continuous=*,
                           integer=*,
                           obj=*,
                           name=*,
                           coefficients=*) \
                           except -1

    cpdef dictionary(self)

    cpdef interactive_lp_problem(self)
