# sage_setup: distribution = sagemath-giac
# distutils: language = c++
# ****************************************************************************
#       Copyright (C) 2012, Frederic Han <frederic.han@imj-prg.fr>
#                     2020, Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.gmp.mpz cimport mpz_t, mpz_set
from libcpp.string cimport string

cdef extern from "giac/giac.h" namespace "giac":
    cdef cppclass context:
         context()

    cdef struct ref_mpz_t:
         pass
    cdef struct ref_real_object:
         pass
    cdef struct ref_complex:
         pass
    cdef struct ref_identificateur:
         pass
    cdef struct ref_symbolic:
         pass
    cdef struct ref_modulo:
         pass
    cdef struct ref_algext:
         pass
    cdef struct giac_float:
         pass
    cdef cppclass Tref_fraction[T]:
         pass
    cdef cppclass Tref_tensor[T]:
         pass
    cdef cppclass vecteur:
         vecteur(int)
         vecteur()
         void push_back(gen &)
         int size()
    cdef struct ref_vecteur:
         pass
    cdef struct ref_sparse_poly1:
         pass
    cdef struct ref_string:
         pass
    cdef struct ref_gen_user:
         pass
    cdef struct ref_gen_map:
         pass
    cdef struct ref_eqwdata:
         pass
    cdef struct ref_grob:
         pass
    cdef struct ref_void_pointer:
         pass

    cdef cppclass gen:
         gen()  except +
         gen(char *, context *) except +
         gen(string , context *) except +
         gen(int )  except +
         gen(long long )  except +
         gen(double ) except +
         #gen(ref_mpz_t * ) except +
         gen(mpz_t & ) except +
         gen(void *ptr,short int subt) except +
         gen(gen )  except +
         gen (vecteur & v,short int s) except +
         gen (ref_vecteur * vptr,short int s) except +

         mpz_t * ref_ZINTptr() except +
         gen * ref_MODptr() except +
         vecteur * ref_VECTptr() except +

         #
         unsigned char type
         signed char subtype
         # (the meaning of types from dispatch.h)
         #    // immediate type (without mem allocation) should be < _ZINT
         #    _INT_= 0, // int val
         #    _DOUBLE_= 1, // double _DOUBLE_val
         #    // all type below or equal to _DOUBLE_ must be non pointers
         #    _ZINT= 2, // mpz_t * _ZINTptr
         #    _REAL= 3, // mpf_t * _REALptr
         #    // all type strictly below _CPLX must be real types
         #    _CPLX= 4, // gen * _CPLXptr
         #    _POLY= 5, // polynome * _POLYptr
         #    _IDNT= 6, // identificateur * _IDNTptr
         #    _VECT= 7, // vecteur * _VECTptr
         #    _SYMB= 8, // symbolic * _SYMBptr
         #    _SPOL1= 9, // sparse_poly1 * _SPOL1ptr
         #    _FRAC= 10, // fraction * _FRACptr
         #    _EXT= 11, // gen * _EXTptr
         #    _STRNG= 12, // string * _STRNGptr
         #    _FUNC= 13, // unary_fonction_ptr * _FUNCptr
         #    _ROOT= 14, // real_complex_rootof *_ROOTptr
         #    _MOD= 15, // gen * _MODptr
         #    _USER= 16, // gen_user * _USERptr
         #    _MAP=17, // map<gen.gen> * _MAPptr
         #    _EQW=18, // eqwdata * _EQWptr
         #    _GROB=19, // grob * _GROBptr
         #    _POINTER_=20, // void * _POINTER_val
         #    _FLOAT_=21 // immediate, _FLOAT_val


         # immediate types
         int val # immediate int (type _INT_)
         double _DOUBLE_val # immediate float (type _DOUBLE_)
         giac_float _FLOAT_val

         # pointer types
         ref_mpz_t * __ZINTptr # long int (type _ZINT)
         ref_real_object * __REALptr # extended double (type _REAL)
         ref_complex * __CPLXptr  # complex as an gen[2] array (type _CPLX)
         ref_identificateur * __IDNTptr # global name identifier (type _IDNT)
         ref_symbolic * __SYMBptr # for symbolic objects (type _SYMB)
         ref_modulo * __MODptr
         ref_algext * __EXTptr # 2 gens for alg. extension (type ext)
         # alg ext: 1st gen is a std::vector or a fraction, 2nd gen is
         # a/ a std::vector, the minimal monic polynomial (the roots are permutable)
         # b/ a real_complex_rootof given by it's min poly and
         # c/ another type meaning that the root is expressed in terms
         #    of another rootof, in this case ext_reduce should be called
         # For 2nd order extension, X^2=d is used if d!=1 mod 4
         # X is the positive solution
         # if d=1 mod 4 the equation is X^2-X=(d-1)/4
         Tref_fraction[gen] * __FRACptr # fraction (type _FRAC)
         Tref_tensor[gen] * __POLYptr  # multidim. sparse polynomials (type poly)
         # _VECTosite types (std::vector<>)
         ref_vecteur * __VECTptr  # vecteur: std::vectors & dense_POLY1 (type _VECT)
         ref_sparse_poly1 * __SPOL1ptr  # std::vector<monome>: sparse 1-d poly (type _SPOL1)
         ref_string * __STRNGptr
         unsigned _FUNC_       # ref_unary_function_ptr * __FUNCptr;
         ref_gen_user * __USERptr
         ref_gen_map * __MAPptr
         ref_eqwdata * __EQWptr
         ref_grob * __GROBptr
         ref_void_pointer * __POINTERptr

         #operators
         gen operator[](int i) except +
         gen operator[](gen & i) except +
         gen operator()(gen & i,context * contextptr) except +
         gen operator()(gen & i,gen & progname,context * contextptr) except +

         gen operator+(gen & b) except +
         gen operator-(gen & b) except +
         gen operator*(gen & b) except +
         gen operator/(gen & b) except +


    gen GIAC_rdiv "rdiv"(gen & a,gen & b) except + # rational division
    gen GIAC_eval "eval" (gen &,int , context *)  except +
    gen GIAC_protecteval "protecteval" (gen , int, context *)  except +
    gen GIAC_pow "pow"(gen & ,gen & , context * ) except +
    gen GIAC_neg "operator-"(gen & ) except +
    gen GIAC_pos "operator+"(gen & ) except +
    gen GIAC_factor "_factor" (gen &, context *) except +
    gen GIAC_factors "_factors" (gen &, context *) except +
    gen GIAC_normal "normal" (gen &, context *)  except +
    gen GIAC_gcd "_gcd" (gen & args, context *) except +
    gen GIAC_smod "_smod" (gen & args, context * ) except +
    gen GIAC_mods "_mods" (gen & args, context * ) except +
    gen GIAC_makemod "_makemod" (gen & , context * ) except +
    string GIAC_print "print" (gen &, context *) except +
    string GIAC_gen2tex "gen2tex" (gen &, context *) except +
    ref_vecteur * GIAC_makenewvecteur "makenewvecteur"(gen & a,gen & b) except +
    gen GIAC_size "_size"(gen & , context *) except +
    gen GIAC_pari_unlock "_pari_unlock"(gen & , context *) except +

    unsigned int GIAC_taille "taille"(gen & , unsigned int) except +
    void GIAC_try_parse_i "try_parse_i"(bool , context *) except +

    string GIAC_giac_aide_dir "giac_aide_dir"() except +
    string GIAC_set_langage "_set_langage"(int , context *) except +

    gen GIAC_sto "sto" (gen &, gen &, bool, context *) except +

    int GIACctrl_c "ctrl_c"

    #test
    gen GIAC_Airy_Ai "_Airy_Ai" (gen &, context *) except +
    gen GIAC_ifactor "_ifactor" (gen &, context *) except +
        

cdef extern from "misc.h":
     void ressetctrl_c() except +
     int testctrl_c() except +
     int giacgencmp( gen & , gen & , context *) except +
     int giacgenrichcmp( gen & , gen & , int, context *) except +
     #NB: we use the following multiplication otherwise some  giac errors make python quit:
     #l=giac([1,2]); l.transpose()*l
     gen GIAC_giacmul "giacmul"( gen & , gen & , context *) except +
     gen GIAC_giacdiv "giacdiv"( gen & , gen & , context *) except +
     gen GIAC_giacmod "giacmod"( gen & , gen & , context *) except +
     #
     string browser_help(gen & , int lang) except +

     void GIAC_archive "archivegen"( string , gen & , context *) except +
     gen GIAC_unarchive "unarchivegen"( string , context *) except +
