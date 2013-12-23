###############################################################################
#   SAGE: Open Source Mathematical Software
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#       Copyright (C) 2008 Burcin Erocal
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

# NOTE: Because of the except+'s below, i.e., C++ exception handling,
# we do *not* have to use sig_on() and sig_off(). We do use it a little
# in the actual pyx code to catch control-c for long running functions.

# distutils: language = c++
# distutils: libraries = pynac gmp

from cpython cimport PyObject

cdef extern from "ginac_wrap.h":
    void ginac_pyinit_Integer(object)
    void ginac_pyinit_Float(object)
    void ginac_pyinit_I(object)

    # forward declaration of GEx
    ctypedef struct GEx "ex"

    ctypedef struct GBasic "basic":
        unsigned int gethash()
        int compare(GBasic other)

    ctypedef struct GConstant "constant":
        unsigned get_serial()

    ctypedef struct GInfinity "infinity":
        bint is_unsigned_infinity()
        bint is_plus_infinity()
        bint is_minus_infinity()
        GEx get_direction()
        GEx conjugate()
        GEx real_part()
        GEx imag_part()

    ctypedef struct GSymbol "symbol":
        unsigned get_domain()
        void set_domain(unsigned d)
        void set_texname(char* t)

    GSymbol* GSymbol_construct_str "Construct_p<symbol, char*>" \
            (void *mem, char* m)

    void GSymbol_destruct "Destruct<symbol>"(GSymbol *mem)

    object GSymbol_to_str "_to_PyString<symbol>"(GSymbol *s)

    ctypedef struct GExPair "std::pair<ex, ex>":
        pass
    ctypedef struct GExMap "exmap":
        void insert(GExPair e)

    ctypedef struct GExListIter "GiNaC::lst::const_iterator":
        void inc "operator++" ()
        GEx obj "operator*" ()
        bint is_not_equal "operator!=" (GExListIter i)

    ctypedef struct GExList "GiNaC::lst":
        GExListIter begin()
        GExListIter end()
        GExList append_sym "append" (GSymbol e)

    ctypedef struct GEx "ex":
        unsigned int gethash()        except +
        int compare(GEx other)        except +
        GEx expand(unsigned int opt)  except +
        GEx collect(GEx s, bint dist) except +
        GEx diff(GSymbol s, int d)    except +
        GEx series(GEx s, int order, unsigned options) except +
        bint is_equal(GEx s)          except +
        bint is_zero()                except +
        bint is_polynomial(GEx vars)  except +
        bint match(GEx pattern, GExList s) except +
        bint find(GEx pattern, GExList s) except +
        bint has(GEx pattern)         except +
        GEx subs(GEx expr)            except +
        GEx subs_map "subs" (GExMap map) except +
        GEx coeff(GEx expr, int n)    except +
        GEx lcoeff(GEx expr)          except +
        GEx tcoeff(GEx expr)          except +
        GEx normal()                  except +
        GEx numer()                   except +
        GEx denom()                   except +
        GEx numer_denom()             except +
        int degree(GEx expr)          except +
        int ldegree(GEx expr)         except +
        GEx unit(GEx expr)            except +
        GEx content(GEx expr)         except +
        GEx primpart(GEx expr)        except +
        void unitcontprim(GEx expr, GEx unit, GEx cont, GEx prim) except +
        GEx rhs()                     except +
        GEx lhs()                     except +
        int nops()                    except +
        GEx op "sorted_op" (int i)    except +
        GEx eval(int level)           except +
        GEx evalf(int level, object parent) except +
        GEx conjugate()               except +
        GEx real_part()               except +
        GEx imag_part()               except +
        bint info(unsigned)           except +
        void dbgprint()
        void dbgprinttree()

    GExPair make_pair "std::make_pair" (GEx, GEx)

    ctypedef struct GNumeric "numeric":
        bint is_positive() except +
        bint is_negative() except +

    # Numericals
    bint is_a_numeric "is_a<numeric>" (GEx e)
    GNumeric ex_to_numeric "ex_to<numeric>" (GEx e)
    # given a GEx that is known to be a numeric, return reference to
    # the underlying PyObject*.
    object py_object_from_numeric(GEx e)     except +

    # Algorithms
    GEx g_gcd "gcd"(GEx a, GEx b) except +

    # Pattern matching wildcards
    GEx g_wild "wild"(unsigned int label) except +

    # Series back to poly
    GEx series_to_poly(GEx e) except +
    bint is_a_series "is_a<pseries>" (GEx e)
    # you must ensure that is_a_series(e) is true before calling this:
    bint g_is_a_terminating_series(GEx e) except +

    # Relations
    ctypedef enum operators "relational::operators":
        equal               "GiNaC::relational::equal"
        not_equal           "GiNaC::relational::not_equal"
        less                "GiNaC::relational::less"
        less_or_equal       "GiNaC::relational::less_or_equal"
        greater             "GiNaC::relational::greater"
        greater_or_equal    "GiNaC::relational::greater_or_equal"

    bint is_negative(GEx x)                  except +
    bint is_a_relational "is_a<relational>" (GEx e)
    bint relational_to_bool(GEx e)
    operators relational_operator(GEx e)
    operators switch_operator(operators op)
    GEx relational(GEx lhs, GEx rhs, operators o)
    GEx g_lt "LT_WRAP" (GEx left, GEx right) except +
    GEx g_eq "EQ_WRAP" (GEx left, GEx right) except +
    GEx g_gt "GT_WRAP" (GEx left, GEx right) except +
    GEx g_le "LE_WRAP" (GEx left, GEx right) except +
    GEx g_ne "NE_WRAP" (GEx left, GEx right) except +
    GEx g_ge "GE_WRAP" (GEx left, GEx right) except +

    #Domains
    unsigned domain_complex "GiNaC::domain::complex"
    unsigned domain_real "GiNaC::domain::real"
    unsigned domain_positive "GiNaC::domain::positive"
    unsigned domain_infinity "GiNaC::domain::infinity"

    # info flags
    unsigned info_real          "GiNaC::info_flags::real"
    unsigned info_rational      "GiNaC::info_flags::rational"
    unsigned info_integer       "GiNaC::info_flags::integer"
    unsigned info_positive      "GiNaC::info_flags::positive"
    unsigned info_negative      "GiNaC::info_flags::negative"
    unsigned info_nonnegative   "GiNaC::info_flags::nonnegative"
    unsigned info_posint        "GiNaC::info_flags::posint"
    unsigned info_negint        "GiNaC::info_flags::negint"
    unsigned info_nonnegint     "GiNaC::info_flags::nonnegint"
    unsigned info_even          "GiNaC::info_flags::even"
    unsigned info_odd           "GiNaC::info_flags::odd"
    unsigned info_rational_function "GiNaC::info_flags::rational_function"

    # Constants
    GEx g_Pi "Pi"
    GEx g_Catalan "Catalan"
    GEx g_Euler "Euler"

    GConstant* GConstant_construct(void *mem, char* name, char* texname, unsigned domain)
    bint is_a_constant "is_a<constant>" (GEx e)
    void GConstant_destruct "Destruct<constant>"(GConstant *mem) except +
    GConstant* GConstant_construct_str "Construct_p<constant, char*>" \
            (void *mem, char* name) except +

    # Infinities
    bint is_a_infinity "is_a<GiNaC::infinity>" (GEx e)
    GEx g_UnsignedInfinity "UnsignedInfinity"
    GEx g_Infinity "Infinity"
    GEx g_mInfinity "-Infinity"
    GInfinity ex_to_infinity "ex_to<GiNaC::infinity>" (GEx e)

    # I is not a constant, but a numeric object
    # we declare it here for easy reference
    GEx g_I "I"

    # Destructor and constructor
    void GEx_destruct "Destruct<ex>"(GEx *mem) except +
    GEx* GEx_construct_symbol "Construct_p<ex, symbol>" \
            (void *mem, GSymbol m) except +
    GEx* GEx_construct_ex "Construct_p<ex, ex>" (void *mem, GEx m) except +
    GEx* GEx_construct_long "Construct_p<ex, long>" (void *mem, long n) except +
    GEx* GEx_construct_double "Construct_p<ex, double>" \
            (void *mem, double d) except +

    GEx* GEx_construct_pyobject "ASSIGN_WRAP" (GEx mem, object n)

    # Conversions
    double GEx_to_double(GEx e, int* success) except +
    object GEx_to_str "_to_PyString<ex>"(GEx *s) except +
    object GEx_to_str_latex "_to_PyString_latex<ex>"(GEx *s) except +

    bint is_a_symbol "is_a<symbol>" (GEx e)
    GSymbol ex_to_symbol "ex_to<symbol>" (GEx e)

    ctypedef struct GParamSetIter "paramset::const_iterator":
        void inc "operator++" ()
        unsigned obj "operator*" ()
        bint is_not_equal "operator!=" (GParamSetIter i)

    ctypedef struct GParamSet "paramset":
        GParamSetIter begin()
        GParamSetIter end()
        int size()

    ctypedef struct GExVector "exvector":
        void push_back(GEx)
        int size()
        GEx at(int i)

    ctypedef struct GExSetIter "std::set<ex, ex_is_less>::const_iterator":
        void inc "operator++" ()
        GEx obj "operator*" ()
        bint is_not_equal "operator!=" (GExSetIter i)

    ctypedef struct GExSet "std::set<ex, ex_is_less>":
        GExSetIter begin()
        GExSetIter end()

    void g_list_symbols "list_symbols" (GEx e, GExSet s)

    # more is_a tests
    bint is_a_add "is_a<add>" (GEx e)
    bint is_a_mul "is_a<mul>" (GEx e)
    bint is_a_power "is_a<power>" (GEx e)
    bint is_a_fderivative "is_a<fderivative>" (GEx e)
    bint is_a_function "is_a<function>" (GEx e)
    bint is_exactly_a_function "is_exactly_a<function>" (GEx e)
    bint is_a_ncmul "is_a<ncmul>" (GEx e)

    # Arithmetic
    int ginac_error()
    GEx gadd "ADD_WRAP" (GEx left, GEx right) except +
    GEx gsub "SUB_WRAP" (GEx left, GEx right) except +
    GEx gmul "MUL_WRAP" (GEx left, GEx right) except +
    GEx gdiv "DIV_WRAP" (GEx left, GEx right) except +
    GEx g_pow "pow" (GEx left, GEx exp)      except +

    GEx g_hold_wrapper "HOLD" (GEx (GEx) except+, GEx ex, bint) except +
    GEx g_hold_wrapper_vec "HOLD" (GEx (GExVector) except+, GExVector vec,
            bint) except +
    GEx g_hold2_wrapper "HOLD2" (GEx (GEx, GEx) except+, GEx ex, GEx ex,
            bint) except +

    #GSymbol get_symbol(char* s)              except +
    GSymbol ginac_symbol "GiNaC::symbol" (char* s, char* t, unsigned d) except +
    GSymbol ginac_new_symbol "GiNaC::symbol" () except +
    GEx g_collect_common_factors "collect_common_factors" (GEx e) except +

    # standard library string
    ctypedef struct stdstring "std::string":
        stdstring assign(char* s, Py_ssize_t l)
        char* c_str()
        unsigned int size()
        char at(unsigned int ind)

    stdstring* stdstring_construct_cstr \
            "new std::string" (char* s, unsigned int l)
    void stdstring_delete "Delete<std::string>"(stdstring* s)

    # Archive
    ctypedef struct GArchive "archive":
        void archive_ex(GEx e, char* name) except +
        GEx unarchive_ex(GExList sym_lst, unsigned ind) except +
        void printraw "printraw(std::cout); " (int t)

    object GArchive_to_str "_to_PyString<archive>"(GArchive *s)
    void GArchive_from_str "_from_str_len<archive>"(GArchive *ar, char* s,
            unsigned int l)


    GEx g_abs "GiNaC::abs" (GEx x)                      except + # absolute value
    GEx g_step "GiNaC::step" (GEx x)                    except + # step function
    GEx g_csgn "GiNaC::csgn" (GEx x)                    except + # complex sign
    GEx g_conjugate "GiNaC::conjugate_function" (GEx x) except + # complex conjugation
    GEx g_real_part "GiNaC::real_part_function" (GEx x) except + # real part
    GEx g_imag_part "GiNaC::imag_part_function" (GEx x) except + # imaginary part
    GEx g_sqrt "GiNaC::sqrt" (GEx x)                    except + # square root (not a GiNaC function, rather an alias for pow(x, numeric(1, 2)))
    GEx g_sin "GiNaC::sin" (GEx x)                      except + # sine
    GEx g_cos "GiNaC::cos" (GEx x)                      except + # cosine
    GEx g_tan "GiNaC::tan" (GEx x)                      except + # tangent
    GEx g_asin "GiNaC::asin" (GEx x)                    except + # inverse sine
    GEx g_acos "GiNaC::acos" (GEx x)                    except + # inverse cosine
    GEx g_atan "GiNaC::atan" (GEx x)                    except + # inverse tangent
    GEx g_atan2 "GiNaC::atan2" (GEx y, GEx x)           except + # inverse tangent with two arguments
    GEx g_sinh "GiNaC::sinh" (GEx x)                    except + # hyperbolic sine
    GEx g_cosh "GiNaC::cosh" (GEx x)                    except + # hyperbolic cosine
    GEx g_tanh "GiNaC::tanh" (GEx x)                    except + # hyperbolic tangent
    GEx g_asinh "GiNaC::asinh" (GEx x)                  except + # inverse hyperbolic sine
    GEx g_acosh "GiNaC::acosh" (GEx x)                  except + # inverse hyperbolic cosine
    GEx g_atanh "GiNaC::atanh" (GEx x)                  except + # inverse hyperbolic tangent
    GEx g_exp "GiNaC::exp" (GEx x)                      except + # exponential function
    GEx g_log "GiNaC::log" (GEx x)                      except + # natural logarithm
    GEx g_Li2 "GiNaC::Li2" (GEx x)                      except + # dilogarithm
    GEx g_Li "GiNaC::Li" (GEx m, GEx x)                 except + # classical polylogarithm as well as multiple polylogarithm
    GEx g_G "GiNaC::G" (GEx a, GEx y)                   except + # multiple polylogarithm
    GEx g_G2 "GiNaC::G" (GEx a, GEx s, GEx y)           except + # multiple polylogarithm with explicit signs for the imaginary parts
    GEx g_S "GiNaC::S" (GEx n, GEx p, GEx x)            except + # Nielsen's generalized polylogarithm
    GEx g_H "GiNaC::H" (GEx m, GEx x)                   except + # harmonic polylogarithm
    GEx g_zeta "GiNaC::zeta" (GEx m)                    except + # Riemann's zeta function as well as multiple zeta value
    GEx g_zeta2 "GiNaC::zeta" (GEx m, GEx s)            except + # alternating Euler sum
    GEx g_zetaderiv "GiNaC::zetaderiv" (GEx n, GEx x)   except + # derivatives of Riemann's zeta function
    GEx g_tgamma "GiNaC::tgamma" (GEx x)                except + # gamma function
    GEx g_lgamma "GiNaC::lgamma" (GEx x)                except + # logarithm of gamma function
    GEx g_beta "GiNaC::beta" (GEx x, GEx y)             except + # beta function (tgamma*tgamma(y)/tgamma(x+y))
    GEx g_psi "GiNaC::psi" (GEx x)                      except + # psi (digamma) function
    GEx g_psi2 "GiNaC::psi" (GEx n, GEx x)              except + # derivatives of psi function (polygamma functions)
    GEx g_factorial "GiNaC::factorial" (GEx n)          except + # factorial function n!
    GEx g_binomial "GiNaC::binomial" (GEx n, GEx k)     except + # binomial coefficients
    GEx g_Order "GiNaC::Order" (GEx x)                  except + # order term function in truncated power series

    # wrapper around arithmetic to allow "hold"ing results
    GEx g_power_construct "GiNaC::power" (GEx b, GEx p) except +
    GEx g_add_construct "GiNaC::add" (GExVector, bint)  except +
    GEx g_mul_construct "GiNaC::mul" (GExVector, bint)  except +

    GEx g_ex1_2 "GiNaC::_ex1_2"

    ctypedef struct GFunction "function":
        unsigned get_serial()
        char* get_name "get_name().c_str" ()

    ctypedef struct GFDerivative "fderivative":
        GParamSet get_parameter_set()

    GFunction ex_to_function "ex_to<function>" (GEx ex)
    GFDerivative ex_to_fderivative "ex_to<fderivative>" (GEx ex)

    GEx g_function_evalv(unsigned int serial, GExVector, bint) except +
    GEx g_function_eval0(unsigned int serial, bint) except +
    GEx g_function_eval1(unsigned int serial, GEx, bint) except +
    GEx g_function_eval2(unsigned int serial, GEx, GEx, bint) except +
    GEx g_function_eval3(unsigned int serial, GEx, GEx, GEx, bint) except +

    ctypedef struct GFunctionOpt "function_options":
        unsigned get_nparams()
        void set_python_func()
        GFunctionOpt eval_func(object f)
        GFunctionOpt evalf_func(object f)
        GFunctionOpt conjugate_func(object f)
        GFunctionOpt real_part_func(object f)
        GFunctionOpt imag_part_func(object f)
        GFunctionOpt derivative_func(object f)
        GFunctionOpt power_func(object f)
        GFunctionOpt series_func(object f)
        GFunctionOpt latex_name(char* name)
        GFunctionOpt do_not_apply_chain_rule()
        GFunctionOpt do_not_evalf_params()
        void set_print_latex_func(object f)
        void set_print_dflt_func(object f)
        char* get_name()
        char* get_latex_name()


    ctypedef struct GFunctionOptVector "vector<function_options>":
        int size()
        GFunctionOpt index "operator[]" (int ind)

    void g_foptions_assign "ASSIGN_WRAP" (GFunctionOpt, GFunctionOpt)

    GFunctionOpt g_function_options "GiNaC::function_options" \
            (char *m)
    GFunctionOpt g_function_options_args "GiNaC::function_options" \
            (char *m, unsigned nargs)
    unsigned g_register_new "GiNaC::function::register_new" (GFunctionOpt opt)

    unsigned find_function "GiNaC::function::find_function" (char* name,
            unsigned nargs) except +ValueError

    GFunctionOptVector g_registered_functions \
            "GiNaC::function::registered_functions" ()

    # these serials allow us to map pynac function objects to
    # Sage special functions for the .operator() method of expressions
    unsigned abs_serial "GiNaC::abs_SERIAL::serial"
    unsigned step_serial "GiNaC::step_SERIAL::serial"# step function
    unsigned csgn_serial "GiNaC::csgn_SERIAL::serial"# complex sign
    unsigned conjugate_serial "GiNaC::conjugate_SERIAL::serial"# complex conjugation
    unsigned real_part_serial "GiNaC::real_part_SERIAL::serial" # real part
    unsigned imag_part_serial "GiNaC::imag_part_SERIAL::serial" # imaginary part
    unsigned sin_serial "GiNaC::sin_SERIAL::serial" # sine
    unsigned cos_serial "GiNaC::cos_SERIAL::serial" # cosine
    unsigned tan_serial "GiNaC::tan_SERIAL::serial" # tangent
    unsigned asin_serial "GiNaC::asin_SERIAL::serial" # inverse sine
    unsigned acos_serial "GiNaC::acos_SERIAL::serial" # inverse cosine
    unsigned atan_serial "GiNaC::atan_SERIAL::serial" # inverse tangent
    unsigned atan2_serial "GiNaC::atan2_SERIAL::serial" # inverse tangent with two arguments
    unsigned sinh_serial "GiNaC::sinh_SERIAL::serial" # hyperbolic sine
    unsigned cosh_serial "GiNaC::cosh_SERIAL::serial" # hyperbolic cosine
    unsigned tanh_serial "GiNaC::tanh_SERIAL::serial" # hyperbolic tangent
    unsigned asinh_serial "GiNaC::asinh_SERIAL::serial" # inverse hyperbolic sine
    unsigned acosh_serial "GiNaC::acosh_SERIAL::serial" # inverse hyperbolic cosine
    unsigned atanh_serial "GiNaC::atanh_SERIAL::serial" # inverse hyperbolic tangent
    unsigned exp_serial "GiNaC::exp_SERIAL::serial" # exponential function
    unsigned log_serial "GiNaC::log_SERIAL::serial" # natural logarithm
    unsigned Li2_serial "GiNaC::Li2_SERIAL::serial" # dilogarithm
    unsigned Li_serial "GiNaC::Li_SERIAL::serial" # classical polylogarithm as well as multiple polylogarithm
    unsigned G_serial "GiNaC::G_SERIAL::serial" # multiple polylogarithm
    #unsigned G2_serial "GiNaC::G_SERIAL::serial" # multiple polylogarithm with explicit signs for the imaginary parts
    unsigned S_serial "GiNaC::S_SERIAL::serial" # Nielsen's generalized polylogarithm
    unsigned H_serial "GiNaC::H_SERIAL::serial" # harmonic polylogarithm
    unsigned zeta1_serial "GiNaC::zeta1_SERIAL::serial" # Riemann's zeta function as well as multiple zeta value
    unsigned zeta2_serial "GiNaC::zeta2_SERIAL::serial" # alternating Euler sum
    unsigned zetaderiv_serial "GiNaC::zetaderiv_SERIAL::serial" # derivatives of Riemann's zeta function
    unsigned tgamma_serial "GiNaC::tgamma_SERIAL::serial" # gamma function
    unsigned lgamma_serial "GiNaC::lgamma_SERIAL::serial" # logarithm of gamma function
    unsigned beta_serial "GiNaC::beta_SERIAL::serial" # beta function (tgamma*tgamma(y)/tgamma(x+y))
    unsigned psi_serial "GiNaC::psi_SERIAL::serial" # psi (digamma) function
    #unsigned psi2_serial "GiNaC::psi_SERIAL::serial" # derivatives of psi function (polygamma functions)
    unsigned factorial_serial "GiNaC::factorial_SERIAL::serial" # factorial function n!
    unsigned binomial_serial "GiNaC::binomial_SERIAL::serial" # binomial coefficients
    unsigned Order_serial "GiNaC::Order_SERIAL::serial" # order term function in truncated power series

    cdef extern from *:
        ctypedef char* const_char_ptr "const char*"
        ctypedef GParamSet const_paramset_ref "const GiNaC::paramset&"

    ctypedef struct py_funcs_struct:
        object (*py_binomial)(object a, object b)
        object (*py_binomial_int)(int n, unsigned int k) except +
        object (*py_gcd)(object a ,object b)
        object (*py_lcm)(object a ,object b)
        object (*py_real)(object a)
        object (*py_imag)(object a)
        object (*py_numer)(object a)
        object (*py_denom)(object a)
        object (*py_conjugate)(object a)
        bint (*py_is_rational)(object a)  except +
        bint (*py_is_crational)(object a)  except +
        bint (*py_is_real)(object a)  except +
        bint (*py_is_integer)(object a)  except +
        bint (*py_is_equal)(object a, object b)  except +
        bint (*py_is_even)(object a)  except +
        bint (*py_is_cinteger)(object a)  except +
        bint (*py_is_prime)(object n)  except +

        object (*py_integer_from_long)(long int x) except +
        object (*py_integer_from_python_obj)(object x) except +

        object (*py_float)(object a, PyObject* parent) except +
        object (*py_RDF_from_double)(double x)


        object (*py_factorial)(object x) except +
        object (*py_doublefactorial)(object x) except +
        object (*py_fibonacci)(object x) except +
        object (*py_step)(object x) except +
        object (*py_bernoulli)(object x) except +
        object (*py_sin)(object x) except +
        object (*py_cos)(object x) except +
        object (*py_zeta)(object x) except +
        object (*py_exp)(object x) except +
        object (*py_log)(object x) except +
        object (*py_tan)(object x) except +
        object (*py_asin)(object x) except +
        object (*py_acos)(object x) except +
        object (*py_atan)(object x) except +
        object (*py_atan2)(object x, object y) except +
        object (*py_sinh)(object x) except +
        object (*py_cosh)(object x) except +
        object (*py_tanh)(object x) except +
        object (*py_asinh)(object x) except +
        object (*py_acosh)(object x) except +
        object (*py_atanh)(object x) except +
        object (*py_tgamma)(object x) except +
        object (*py_lgamma)(object x) except +
        object (*py_isqrt)(object x) except +
        object (*py_sqrt)(object x) except +
        object (*py_abs)(object x) except +
        object (*py_mod)(object x, object y) except +
        object (*py_smod)(object x, object y) except +
        object (*py_irem)(object x, object y) except +
        object (*py_iquo)(object x, object y) except +
        object (*py_iquo2)(object x, object y) except +
        object (*py_li)(object x, object n, object parent) except +
        object (*py_li2)(object x) except +
        object (*py_psi)(object x) except +
        object (*py_psi2)(object n, object x) except +

        int    (*py_int_length)(object x) except -1

        object (*py_eval_constant)(unsigned serial, object parent)
        object (*py_eval_unsigned_infinity)()
        object (*py_eval_infinity)()
        object (*py_eval_neg_infinity)()

        int    (*py_get_parent_char)(object o) except -1

        stdstring* (*py_latex)(object o, int level) except +
        stdstring* (*py_repr)(object o, int level) except +

        stdstring* (*py_dumps)(object o) except +
        object (*py_loads)(object s) except +

        object exvector_to_PyTuple(GExVector seq)
        GEx pyExpression_to_ex(object res) except *
        object ex_to_pyExpression(GEx juice)
        int py_get_ginac_serial()

        object py_get_sfunction_from_serial(unsigned s) except +
        unsigned py_get_serial_from_sfunction(object f) except +
        unsigned py_get_serial_for_new_sfunction(stdstring &s, unsigned nargs) \
                except +

        stdstring* py_print_function(unsigned id, object args) except +
        stdstring* py_latex_function(unsigned id, object args) except +


        GConstant py_get_constant(const_char_ptr name) except +

        stdstring* (*py_print_fderivative)(unsigned id, object params, object args) except +
        stdstring* (*py_latex_fderivative)(unsigned id, object params, object args) except +
        object (*paramset_to_PyTuple)(const_paramset_ref s)

        object (*py_rational_power_parts)(object basis, object exp)

    py_funcs_struct py_funcs "GiNaC::py_funcs"

cdef extern from "pynac/order.h":
    bint print_order_compare "GiNaC::print_order().compare" \
            (GEx left, GEx right) except +
    bint print_order_compare_mul "GiNaC::print_order_mul().compare" \
            (GEx left, GEx right) except +
    bint print_order "GiNaC::print_order()" \
            (GEx left, GEx right) except +
    bint print_order_mul "GiNaC::print_order_mul()" \
            (GEx left, GEx right) except +
