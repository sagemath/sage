# sage_setup: distribution = sagemath-symbolics
"""
Declarations for pynac, a Python frontend for ginac

Check that we can externally cimport this (:issue:`18825`)::

    sage: cython(                                                                       # needs sage.misc.cython
    ....: '''
    ....: cimport sage.symbolic.expression
    ....: ''')
"""

# ****************************************************************************
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#       Copyright (C) 2008 Burcin Erocal
#       Copyright (C) 2017 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cpython.object cimport PyObject
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.string cimport string as stdstring
from sage.libs.gmp.types cimport mpz_t, mpq_t, mpz_ptr, mpq_ptr

cdef extern from "pynac_wrap.h":
    void ginac_pyinit_Integer(object)
    void ginac_pyinit_Float(object)
    void ginac_pyinit_I(object)

    # forward declarations
    cdef cppclass GEx "ex"
    cdef cppclass GExprSeq "exprseq"

    cdef cppclass GBasic "basic":
        long gethash()
        int compare(GBasic other)

    cdef cppclass GConstant "constant":
        GConstant()
        GConstant(char* name, void* evalf, char* texname, unsigned domain)
        unsigned get_serial()

    # This is actually a function, but we declare it as void* for
    # simplicity.
    void* ConstantEvalf

    cdef cppclass GInfinity "infinity":
        bint is_unsigned_infinity()
        bint is_plus_infinity()
        bint is_minus_infinity()
        GEx get_direction()
        GEx conjugate()
        GEx real_part()
        GEx imag_part()

    cdef cppclass GSymbol "symbol":
        unsigned get_domain()
        void set_domain(unsigned d)
        void set_texname(char* t)
        const char* get_name "get_name().c_str" ()

    cdef cppclass GExPair "std::pair<ex, ex>":
        pass

    cdef cppclass GExMap "exmap":
        void insert(GExPair e)

    cdef cppclass GExListIter "GiNaC::lst::const_iterator":
        void inc "operator++" ()
        GEx obj "operator*" ()
        bint operator!=(GExListIter i)

    cdef cppclass GExList "GiNaC::lst":
        GExListIter begin()
        GExListIter end()
        GExList append_sym "append" (GSymbol e)

    cdef cppclass GSymbolSetIter "GiNaC::symbolset::const_iterator":
        void inc "operator++" ()
        GEx obj "operator*" ()
        bint operator!=(GSymbolSetIter i)

    cdef cppclass GSymbolSet "GiNaC::symbolset":
        GSymbolSetIter begin()
        GSymbolSetIter end()

    cdef cppclass GEx "ex":
        GEx()
        GEx(GNumeric o)
        GEx(GSymbol m)
        GEx(GEx m)
        GEx(long n)
        GEx(double d)
        GEx(GExprSeq s)
        GEx operator=(object o)
        long gethash()                except +
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
        GSymbolSet free_symbols()     except +
        bint has(GEx pattern)         except +
        GEx subs(GEx expr)            except +
        GEx subs_map "subs" (GExMap map, unsigned options) except +
        GEx coeff(GEx expr, GEx n)    except +
        GEx lcoeff(GEx expr)          except +
        GEx tcoeff(GEx expr)          except +
        void coefficients(GEx s, vector[pair[GEx,GEx]]) except +
        GEx combine_fractions(bint deep) except +
        GEx normal(int level, bint noexpand_combined, bint noexpand_frac) except +
        GEx numer()                   except +
        GEx denom()                   except +
        GEx numer_denom()             except +
        GNumeric degree(GEx expr)          except +
        GNumeric ldegree(GEx expr)         except +
        GEx unit(GEx expr)            except +
        GEx content(GEx expr)         except +
        GEx primpart(GEx expr)        except +
        void unitcontprim(GEx expr, GEx unit, GEx cont, GEx prim) except +
        GEx rhs()                     except +
        GEx lhs()                     except +
        int nops()                    except +
        GEx op "sorted_op" (int i)    except +
        GEx eval(int level)           except +
        GEx evalf(int level, parent)  except +
        GEx conjugate()               except +
        GEx real_part()               except +
        GEx imag_part()               except +
        bint info(unsigned)           except +
        void dbgprint()
        void dbgprinttree()

    GExPair make_pair "std::make_pair" (GEx, GEx)

    cdef cppclass GNumeric "numeric":
        bint is_positive() except +
        bint is_negative() except +

    # Numericals
    bint is_a_numeric "is_a<numeric>" (GEx e)
    GNumeric ex_to_numeric "ex_to<numeric>" (GEx e)
    # given a GEx that is known to be a numeric, return reference to
    # the underlying PyObject*.
    py_object_from_numeric(GEx e)     except +

    # Algorithms
    GEx g_gcd "gcd"(GEx a, GEx b) except +
    bint g_factor "factor"(GEx expr, GEx res) except +
    GEx g_gosper_term "gosper_term"(GEx the_ex, GEx n) except +
    GEx g_gosper_sum_definite "gosper_sum_definite"(GEx the_ex,
            GEx n, GEx a, GEx b, int* p) except +
    GEx g_gosper_sum_indefinite "gosper_sum_indefinite"(GEx the_ex,
            GEx n, int* p) except +
    GEx to_gamma(GEx expr)          except +
    GEx gamma_normalize(GEx expr)   except +
    GEx g_resultant "resultant"(GEx a, GEx b, GEx v) except +

    # Pattern matching wildcards
    GEx g_wild "wild"(unsigned int label) except +
    bint haswild(GEx x) except +

    # Series back to poly
    GEx series_to_poly(GEx e) except +
    bint is_a_series "is_a<pseries>" (GEx e)
    # you must ensure that is_a_series(e) is true before calling this:
    bint g_is_a_terminating_series(GEx e) except +
    GEx g_series_var(GEx e) except +

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
    unsigned decide_relational(GEx e) except +
    operators relational_operator(GEx e)
    GEx relational(GEx lhs, GEx rhs, operators o)
    GEx operator<(GEx left, GEx right) except +
    GEx operator==(GEx left, GEx right) except +
    GEx operator>(GEx left, GEx right) except +
    GEx operator<=(GEx left, GEx right) except +
    GEx operator!=(GEx left, GEx right) except +
    GEx operator>=(GEx left, GEx right) except +

    # Domains
    unsigned domain_complex "GiNaC::domain::complex"
    unsigned domain_real "GiNaC::domain::real"
    unsigned domain_positive "GiNaC::domain::positive"
    unsigned domain_infinity "GiNaC::domain::infinity"
    unsigned domain_integer "GiNaC::domain::integer"

    # relational outcomes
    unsigned relational_true "GiNaC::relational::result::True"
    unsigned relational_false "GiNaC::relational::result::False"
    unsigned relational_undecidable "GiNaC::relational::result::undecidable"
    unsigned relational_notimplemented "GiNaC::relational::result::notimplemented"

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

    # assumptions
    void pynac_assume_rel "GiNaC::assume" (GEx rel)
    void pynac_assume_gdecl "GiNaC::assume" (GEx x, char*)
    void pynac_forget_rel "GiNaC::forget" (GEx rel)
    void pynac_forget_gdecl "GiNaC::forget" (GEx x, char*)

    # Constants
    GEx g_Pi "Pi"
    GEx g_Catalan "Catalan"
    GEx g_Euler "Euler"
    GEx g_NaN "NaN"

    bint is_a_constant "is_a<constant>" (GEx e)

    # Infinities
    bint is_a_infinity "is_a<GiNaC::infinity>" (GEx e)
    GEx g_UnsignedInfinity "UnsignedInfinity"
    GEx g_Infinity "Infinity"
    GEx g_mInfinity "-Infinity"
    GInfinity ex_to_infinity "ex_to<GiNaC::infinity>" (GEx e)

    # I is not a constant, but a numeric object
    # we declare it here for easy reference
    GEx g_I "I"

    # Conversions
    double GEx_to_double(GEx e, int* success) except +
    GEx_to_str_latex "_to_PyString_latex<ex>"(GEx *s) except +

    bint is_a_symbol "is_a<symbol>" (GEx e)
    GSymbol ex_to_symbol "ex_to<symbol>" (GEx e)

    cdef cppclass GParamSetIter "paramset::const_iterator":
        void inc "operator++" ()
        unsigned obj "operator*" ()
        bint operator!=(GParamSetIter i)

    cdef cppclass GParamSet "paramset":
        GParamSetIter begin()
        GParamSetIter end()
        int size()

    cdef cppclass GExVector "exvector":
        void push_back(GEx)
        int size()
        GEx at(int i)

    cdef cppclass GExprSeq "exprseq":
        GExprSeq()
        GExprSeq(GExVector m)

    bint is_a_exprseq "is_a<exprseq>" (GEx e)
    bint is_exactly_a_exprseq "is_exactly_a<exprseq>" (GEx e)

    cdef cppclass GExSetIter "std::set<ex, ex_is_less>::const_iterator":
        void inc "operator++" ()
        GEx obj "operator*" ()
        bint operator!=(GExSetIter i)

    cdef cppclass GExSet "std::set<ex, ex_is_less>":
        GExSetIter begin()
        GExSetIter end()

    void g_list_symbols "list_symbols" (GEx e, GExSet s)

    # more is_a tests
    bint is_a_add "is_a<GiNaC::add>" (GEx e)
    bint is_a_mul "is_a<GiNaC::mul>" (GEx e)
    bint is_a_power "is_a<GiNaC::power>" (GEx e)
    bint is_a_fderivative "is_a<GiNaC::fderivative>" (GEx e)
    bint is_a_function "is_a<GiNaC::function>" (GEx e)
    bint is_exactly_a_function "is_exactly_a<GiNaC::function>" (GEx e)

    # Arithmetic
    int ginac_error()
    GEx operator+(GEx left, GEx right) except +
    GEx operator-(GEx left, GEx right) except +
    GEx operator*(GEx left, GEx right) except +
    GEx operator/(GEx left, GEx right) except +
    GEx g_pow "pow" (GEx left, GEx exp) except +

    GEx g_hold_wrapper "HOLD" (GEx (GEx) except+, GEx ex, bint) except +
    GEx g_hold_wrapper_vec "HOLD" (GEx (GExVector) except+, GExVector vec,
            bint) except +
    GEx g_hold2_wrapper "HOLD2" (GEx (GEx, GEx) except+, GEx ex, GEx ex,
            bint) except +
    void g_set_state "GiNaC::set_state" (stdstring & s, bint b) except +

    GSymbol ginac_symbol "GiNaC::symbol" (char* s, char* t, unsigned d) except +
    GSymbol ginac_new_symbol "GiNaC::symbol" () except +
    GEx g_collect_common_factors "collect_common_factors" (GEx e) except +

    # Archive
    cdef cppclass GArchive "archive":
        void archive_ex(GEx e, char* name) except +
        GEx unarchive_ex(GExList sym_lst, unsigned ind) except +
        void printraw "printraw(std::cout); " (int t)

    GEx g_abs "GiNaC::abs" (GEx x)                      except + # absolute value
    GEx g_step "GiNaC::unit_step" (GEx x)               except + # step function
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
    GEx g_stieltjes "GiNaC::stieltjes" (GEx m)          except + # Stieltjes constants
    GEx g_zetaderiv "GiNaC::zetaderiv" (GEx n, GEx x)   except + # derivatives of Riemann's zeta function
    GEx g_gamma "GiNaC::gamma" (GEx x)                  except + # gamma function
    GEx g_lgamma "GiNaC::lgamma" (GEx x)                except + # logarithm of gamma function
    GEx g_beta "GiNaC::beta" (GEx x, GEx y)             except + # beta function (gamma(x)*gamma(y)/gamma(x+y))
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

    cdef cppclass GFunction "function":
        unsigned get_serial()
        char* get_name "get_name().c_str" ()

    cdef cppclass GFDerivative "fderivative":
        GParamSet get_parameter_set()

    GFunction ex_to_function "ex_to<function>" (GEx ex)
    GFDerivative ex_to_fderivative "ex_to<fderivative>" (GEx ex)

    GEx g_function_evalv(unsigned int serial, GExVector, bint) except +
    GEx g_function_eval0(unsigned int serial, bint) except +
    GEx g_function_eval1(unsigned int serial, GEx, bint) except +
    GEx g_function_eval2(unsigned int serial, GEx, GEx, bint) except +
    GEx g_function_eval3(unsigned int serial, GEx, GEx, GEx, bint) except +

    cdef cppclass GFunctionOpt "function_options":
        GFunctionOpt operator=(GFunctionOpt)
        unsigned get_nparams()
        void set_python_func()
        GFunctionOpt eval_func(f)
        GFunctionOpt subs_func(f)
        GFunctionOpt evalf_func(f)
        GFunctionOpt conjugate_func(f)
        GFunctionOpt real_part_func(f)
        GFunctionOpt imag_part_func(f)
        GFunctionOpt derivative_func(f)
        GFunctionOpt power_func(f)
        GFunctionOpt series_func(f)
        GFunctionOpt latex_name(char* name)
        GFunctionOpt do_not_apply_chain_rule()
        GFunctionOpt do_not_evalf_params()
        void set_print_latex_func(f)
        void set_print_dflt_func(f)
        char* get_name()
        char* get_latex_name()

    cdef cppclass GFunctionOptVector "vector<function_options>":
        int size()
        GFunctionOpt index "operator[]" (int ind)

    void g_foptions_assign "ASSIGN_WRAP" (GFunctionOpt, GFunctionOpt)

    GFunctionOpt g_function_options "GiNaC::function_options" (char *m)
    GFunctionOpt g_function_options_args "GiNaC::function_options" \
        (char *m, unsigned nargs)
    unsigned g_register_new "GiNaC::function::register_new" (GFunctionOpt opt)

    unsigned find_function "GiNaC::function::find_function" (char* name,
            unsigned nargs) except +

    bint has_symbol "GiNaC::has_symbol" (GEx ex)
    bint has_symbol_or_function "GiNaC::has_symbol_or_function" (GEx ex)

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
    # unsigned G2_serial "GiNaC::G_SERIAL::serial" # multiple polylogarithm with explicit signs for the imaginary parts
    unsigned S_serial "GiNaC::S_SERIAL::serial" # Nielsen's generalized polylogarithm
    unsigned H_serial "GiNaC::H_SERIAL::serial" # harmonic polylogarithm
    unsigned zeta1_serial "GiNaC::zeta1_SERIAL::serial" # Riemann's zeta function as well as multiple zeta value
    unsigned zeta2_serial "GiNaC::zeta2_SERIAL::serial" # alternating Euler sum
    unsigned stieltjes1_serial "GiNaC::stieltjes1_SERIAL::serial" # Stieltjes constants
    unsigned zetaderiv_serial "GiNaC::zetaderiv_SERIAL::serial" # derivatives of Riemann's zeta function
    unsigned tgamma_serial "GiNaC::tgamma_SERIAL::serial" # gamma function
    unsigned lgamma_serial "GiNaC::lgamma_SERIAL::serial" # logarithm of gamma function
    unsigned beta_serial "GiNaC::beta_SERIAL::serial" # beta function (tgamma(x)*tgamma(y)/tgamma(x+y))
    unsigned psi_serial "GiNaC::psi_SERIAL::serial" # psi (digamma) function
    # unsigned psi2_serial "GiNaC::psi_SERIAL::serial" # derivatives of psi function (polygamma functions)
    unsigned factorial_serial "GiNaC::factorial_SERIAL::serial" # factorial function n!
    unsigned binomial_serial "GiNaC::binomial_SERIAL::serial" # binomial coefficients
    unsigned Order_serial "GiNaC::Order_SERIAL::serial" # order term function in truncated power series

    ctypedef GParamSet const_paramset_ref "const GiNaC::paramset&"

    ctypedef struct py_funcs_struct:
        py_gcd(a, b)
        py_lcm(a, b)
        py_real(a)
        py_imag(a)
        py_numer(a)
        py_denom(a)
        bint py_is_rational(a)
        bint py_is_real(a)
        bint py_is_integer(a)
        bint py_is_equal(a, b)
        bint py_is_even(a)
        bint py_is_prime(n)
        bint py_is_exact(x)

        py_integer_from_long(long int x)
        py_integer_from_python_obj(x)
        py_integer_from_mpz(mpz_t)
        py_rational_from_mpq(mpq_t)
        bint py_is_Integer(x)
        bint py_is_Rational(x)
        mpz_ptr py_mpz_from_integer(x)
        mpq_ptr py_mpq_from_rational(x)

        py_float(a, PyObject* kwds)
        py_RDF_from_double(double x)

        py_factorial(x)
        py_doublefactorial(x)
        py_fibonacci(x)
        py_step(x)
        py_bernoulli(x)
        py_sin(x)
        py_cos(x)
        py_stieltjes(x)
        py_zeta(x)
        py_exp(x)
        py_log(x)
        py_tan(x)
        py_asin(x)
        py_acos(x)
        py_atan(x)
        py_atan2(x, y)
        py_sinh(x)
        py_cosh(x)
        py_tanh(x)
        py_asinh(x)
        py_acosh(x)
        py_atanh(x)
        py_isqrt(x)
        py_sqrt(x)
        py_mod(x, y)
        py_smod(x, y)
        py_irem(x, y)
        py_psi(x)
        py_psi2(n, x)

        py_eval_constant(unsigned serial, parent)
        py_eval_unsigned_infinity()
        py_eval_infinity()
        py_eval_neg_infinity()

        int py_get_parent_char(o) except -1

        stdstring* py_latex(o, int level)
        stdstring* py_repr(o, int level)

        stdstring* py_dumps(o)
        py_loads(s)

        exvector_to_PyTuple(GExVector seq)
        GEx pyExpression_to_ex(res) except *
        ex_to_pyExpression(GEx juice)
        int py_get_ginac_serial()
        subs_args_to_PyTuple(GExMap map, unsigned options, GExVector seq)

        py_get_sfunction_from_serial(unsigned s)
        unsigned py_get_serial_from_sfunction(f)
        unsigned py_get_serial_for_new_sfunction(stdstring &s, unsigned nargs)

        stdstring* py_print_function(unsigned id, args)
        stdstring* py_latex_function(unsigned id, args)

        GConstant py_get_constant(const char* name)

        stdstring* py_print_fderivative(unsigned id, params, args)
        stdstring* py_latex_fderivative(unsigned id, params, args)
        paramset_to_PyTuple(const_paramset_ref s)

    py_funcs_struct py_funcs "GiNaC::py_funcs"

cdef extern from "ginac/order.h":
    bint print_order_compare "GiNaC::print_order().compare" \
        (GEx left, GEx right) except +
    bint print_order_compare_mul "GiNaC::print_order_mul().compare" \
        (GEx left, GEx right) except +
