# distutils: extra_compile_args = ECL_CFLAGS
# distutils: include_dirs = ECL_INCDIR
# distutils: libraries = ECL_LIBRARIES
# distutils: library_dirs = ECL_LIBDIR
# distutils: extra_link_args = ECL_LIBEXTRA

###############################################################################
#   Sage: Open Source Mathematical Software
#       Copyright (C) 2009 Nils Bruin <nbruin@sfu.ca>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

#ecl's header files export a very large number of definitions, referred to here
#as ecl.h bindings. In reality ecl.h includes a couple of other header files.
#We only include the cython translations of the ones we need here.

#ecl's naming conventions have been followed as much as possible. These
#conventions are not entirely consistently used in ECL itself, though, and
#those cases are corrected in the cython bindings.
#(see cl_cons and ecl_fixnum)
# cl_*  functions are proper common lisp routines. In particular, input and
#       return types are cl_object
# ecl_* functions provide other interfaces. Arguments or return types might
#       include C types.
#
# In addition, type predicates that return a C boolean are named:
# bint_* bindings corresponding to a mix of macros and ecl_* functions in ecl.h

cdef extern from "ecl/ecl.h":

    # Typedefs

    ctypedef long int cl_fixnum
    ctypedef cl_fixnum cl_narg
    ctypedef void *cl_object
    ctypedef unsigned long int cl_index
    ctypedef int ecl_character

    ctypedef enum ecl_option:
        ECL_OPT_INCREMENTAL_GC = 0,
        ECL_OPT_TRAP_SIGSEGV,
        ECL_OPT_TRAP_SIGFPE,
        ECL_OPT_TRAP_SIGINT,
        ECL_OPT_TRAP_SIGILL,
        ECL_OPT_TRAP_SIGBUS,
        ECL_OPT_TRAP_SIGPIPE,
        ECL_OPT_TRAP_INTERRUPT_SIGNAL,
        ECL_OPT_SIGNAL_HANDLING_THREAD,
        ECL_OPT_SIGNAL_QUEUE_SIZE,
        ECL_OPT_BOOTED,
        ECL_OPT_BIND_STACK_SIZE,
        ECL_OPT_BIND_STACK_SAFETY_AREA,
        ECL_OPT_FRAME_STACK_SIZE,
        ECL_OPT_FRAME_STACK_SAFETY_AREA,
        ECL_OPT_LISP_STACK_SIZE,
        ECL_OPT_LISP_STACK_SAFETY_AREA,
        ECL_OPT_C_STACK_SIZE,
        ECL_OPT_C_STACK_SAFETY_AREA,
        ECL_OPT_HEAP_SIZE,
        ECL_OPT_HEAP_SAFETY_AREA,
        ECL_OPT_THREAD_INTERRUPT_SIGNAL,
        ECL_OPT_SET_GMP_MEMORY_FUNCTIONS,
        ECL_OPT_LIMIT

    # boot and shutdown

    cl_fixnum ecl_get_option(int option)
    void ecl_set_option(int option, cl_fixnum value)
    void cl_boot(int argc, char **argv)
    void cl_shutdown()

    # predefined symbols

    cl_object ECL_NIL
    cl_object ECL_T
    cl_fixnum MOST_POSITIVE_FIXNUM
    cl_fixnum MOST_NEGATIVE_FIXNUM

    # Type predicates returning a cl_object

    cl_object cl_symbolp(cl_object o)
    cl_object cl_numberp(cl_object o)
    cl_object cl_integerp(cl_object x)
    cl_object cl_rationalp(cl_object x)
    cl_object cl_floatp(cl_object x)

    # Type predicates returning a C boolean

    bint bint_floatp "floatp" (cl_object x)
    bint bint_numberp "ecl_numberp" (cl_object x)
    bint bint_eql "ecl_eql"(cl_object x, cl_object y)
    bint bint_equal "ecl_equal"(cl_object x, cl_object y)
    bint bint_equalp "ecl_equalp"(cl_object x, cl_object y)
    bint bint_stringp "ecl_stringp"(cl_object x)
    bint bint_fixnump "FIXNUMP"(cl_object o)
    bint bint_characterp "CHARACTERP"(cl_object o)
    bint bint_nullp "Null"(cl_object o)
    bint bint_listp "LISTP" (cl_object o)
    bint bint_consp "CONSP" (cl_object o)
    bint bint_atomp "ATOM" (cl_object o)

    # Equality tests

    cl_object cl_eq(cl_object x, cl_object y)
    cl_object cl_eql(cl_object x, cl_object y)
    cl_object cl_equal(cl_object x, cl_object y)

    # ECL numeric type conversion

    cl_object ecl_make_fixnum(cl_fixnum i)
    cl_object ecl_make_integer(cl_fixnum i)
    cl_object ecl_make_unsigned_integer(cl_index i)
    cl_fixnum ecl_fixint "fixint" (cl_object x)

    cl_object ecl_make_ratio(cl_object num, cl_object den)
    cl_object cl_numerator(cl_object x)
    cl_object cl_denominator(cl_object x)

    cl_object ecl_make_doublefloat(double f)
    double ecl_to_double(cl_object x)

    # list manipulation

    cl_object cl_cons "ecl_cons" (cl_object a, cl_object d)
    cl_object cl_car(cl_object x)
    cl_object cl_cdr(cl_object x)
    cl_object cl_caar(cl_object x)
    cl_object cl_cadr(cl_object x)
    cl_object cl_cdar(cl_object x)
    cl_object cl_cddr(cl_object x)
    cl_object cl_rplaca(cl_object x, cl_object v)
    cl_object cl_rplacd(cl_object x, cl_object v)
    cl_object ecl_list1(cl_object a)

    # string parsing and string IO

    char *ecl_base_string_pointer_safe(cl_object f)
    cl_object ecl_read_from_cstring(char *s)
    cl_object ecl_read_from_cstring_safe(char *s, cl_object err)
    cl_object cl_write_to_string(cl_narg narg, cl_object o)
    cl_object ecl_cstring_to_base_string_or_nil(char *s)
    cl_object si_coerce_to_base_string(cl_object x)
    cl_object si_base_string_p(cl_object x)
    cl_object cl_make_string(cl_narg narg, cl_object o, ...)
    ecl_character ecl_char(cl_object s, cl_index i)
    ecl_character ecl_char_set(cl_object s, cl_index i, ecl_character c)

    # S-expr evaluation and function calls

    cl_object cl_safe_eval(cl_object form, cl_object env, cl_object value)
    cl_object cl_eval(cl_object form)
    cl_object cl_funcall(cl_narg narg, cl_object fun, cl_object arg1,...)
    cl_object cl_apply(cl_narg narg, cl_object fun, cl_object args)
    cl_object cl_set(cl_object var, cl_object val)
    int ecl_nvalues "NVALUES"
    cl_object ecl_values "VALUES"(int n)

    # Common Lisp "EQUAL" compatible hash function

    cl_object cl_sxhash(cl_object key)

    # symbols

    cl_object ecl_make_symbol(const char *name, const char *package)

    # sequences

    cl_fixnum ecl_length(cl_object x)
