//
//  farey.cpp
//
//  Implementation of FareySymbol
//
//****************************************************************************
//       Copyright (C) 2011 Hartmut Monien <monien@th.physik.uni-bonn.de>
//                          Stefan Krämer <skraemer@th.physik.uni-bonn.de>
//
//  Distributed under the terms of the GNU General Public License (GPL)
//
//    This code is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    General Public License for more details.
//
//  The full text of the GPL is available at:
//
//                  http://www.gnu.org/licenses/
//****************************************************************************

#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <functional>

#include <gmpxx.h>
#include <Python.h>

#include "farey.hpp"
#include "sage/modular/arithgroup/farey_symbol.h"


using namespace std;

inline
ostream& tab(ostream& os) { os << "\t"; return os; }

template <class T>
inline
ostream& operator<<(ostream& os, const vector<T>& v) {
  os << v.size() << " ";
  for(typename vector<T>::const_iterator i=v.begin(); i!=v.end(); i++) {
    os << *i << " ";
  }
  return os;
}

template <class T>
inline
istream& operator>>(istream& is, vector<T>& v) {
  size_t n;
  is >> n;
  for(size_t i=0; i<n; i++) {
    T tmp;
    is >> tmp;
    v.push_back(tmp);
  }
  return is;
}

template <>
inline
istream& operator>>(istream& is, vector<SL2Z>& v) {
  size_t n;
  is >> n;
  for(size_t i=0; i<n; i++) {
    SL2Z tmp(1, 0, 0, 1);
    is >> tmp;
    v.push_back(tmp);
  }
  return is;
}

inline
istream& operator>>(istream& is, FareySymbol& F) {
  is >> F.pairing_max
     >> F.pairing
     >> F.cusp_classes
     >> F.a
     >> F.b
     >> F.x
     >> F.coset
     >> F.generators
     >> F.cusps
     >> F.cusp_widths
     >> F.reductions
     >> F.even
     >> F.pairing_in_group;
  return is;
}

inline
ostream& operator<<(ostream& os, const FareySymbol& F) {
  os << F.pairing_max << " "
     << F.pairing
     << F.cusp_classes
     << F.a
     << F.b
     << F.x
     << F.coset
     << F.generators
     << F.cusps
     << F.cusp_widths
     << F.reductions
     << F.even << " "
     << F.pairing_in_group;
  return os;
}

inline
mpq_class operator/(const mpz_class& a, const mpz_class& b) {
  mpq_class result(a, b);
  result.canonicalize();
  return result;
}

inline
mpq_class
operator*(const SL2Z& M, const mpq_class& z) {
  mpz_class p = z.get_num(), q = z.get_den();
  if( M.c()*p+M.d()*q == 0 ) {
    throw string(__FUNCTION__) + ": division by zero.";
  }
  mpq_class result(M.a()*p+M.b()*q, M.c()*p+M.d()*q);
  result.canonicalize();
  return result;
}

inline
vector<mpq_class>
operator*(const SL2Z& M, const vector<mpq_class>& v) {
  vector<mpq_class> result;
  for(size_t j=0; j<v.size(); j++) result.push_back(M*v[j]);
  return result;
}

inline
mpz_class
floor(const mpq_class r) {
  mpz_class result = mpz_class(r.get_num()/r.get_den());
  if( r >= 0 ) {
    return result;
  } else {
    return result - 1;
  }
}

inline
mpz_class lcm(const mpz_class& a, const mpz_class& b) {
  mpz_class result;
  mpz_lcm(result.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
  return result;
}

inline
mpz_class lcm(const vector<mpz_class>& v) {
  mpz_class q(1);
  for(size_t i=0; i<v.size(); i++) {
    q = lcm(q, v[i]);
  }
  return q;
}

inline
mpz_class gcd(const mpz_class& a, const mpz_class& b) {
  mpz_class result;
  mpz_gcd( result.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
  return result;
}

inline
void
gcd_ext(mpz_class& g, mpz_class& r, mpz_class& s, const mpz_class& a, const mpz_class& b) {
  mpz_gcdext(g.get_mpz_t(),
	     r.get_mpz_t(), s.get_mpz_t(),
	     a.get_mpz_t(), b.get_mpz_t());
}

inline
/* cf. Rademacher, The Fourier Coefficients of the Modular Invariant
   J(\tau), p. 502 */
SL2Z rademacher_matrix(const mpq_class& q) {
  mpz_class h = q.get_num(), k = q.get_den(), j;
  mpz_class g, r, s;
  gcd_ext(g, r, s, h, k);
  if( r < 0 ) {
    j = -r;
  } else {
    j = k-r;
  }
  return SL2Z(j, -(h*j+1)/k, k, -h);
}

//--- Helper class for checking membership of SL2Z matrix in group GammaH ----

is_element_GammaH::is_element_GammaH(int p_, PyObject* gen_list) : p(p_) {
  typedef vector<long>::const_iterator const_iterator;
  // list of generators
  vector<long> gen;
  Py_ssize_t ngen = PyList_Size(gen_list);
  for(Py_ssize_t i=0; i<ngen; i++) {
    PyObject* item = PyList_GetItem(gen_list, i);
    gen.push_back(convert_to_long(item));
  }
  // generate H from generators
  H = gen;
  for(;;) {
    vector<long> m;
    for(const_iterator i=gen.begin(); i!=gen.end(); i++) {
      for(const_iterator j=H.begin(); j!=H.end(); j++) {
        long q = ((*i)*(*j))%p;
        if( find(H.begin(), H.end(), q) == H.end() and
            find(m.begin(), m.end(), q) == m.end() ) {
          m.push_back(q);
        }
      }
    }
    if( m.size() == 0 ) break;
    for(const_iterator i=m.begin(); i!=m.end(); i++) H.push_back(*i);
  }
  // sort for binary searches
  sort(H.begin(), H.end());
}

is_element_GammaH::~is_element_GammaH() {
}

bool is_element_GammaH::is_member(const SL2Z& m) const {
  mpz_class a = m.a()%p; if(a < 0) a+=p;
  mpz_class d = m.d()%p; if(d < 0) d+=p;
  if( m.c()%p != 0 ) return false;
  if( not binary_search(H.begin(), H.end(), a.get_si()) ) return false;
  if( not binary_search(H.begin(), H.end(), d.get_si()) ) return false;
  return true;
}

//--- Helper class for checking membership of SL2Z matrix in group -------------
// group defined by the python object.

is_element_general::is_element_general(PyObject* group_) : group(group_) {
  if( PyObject_HasAttrString(group, "__contains__") ) {
    method = PyObject_GetAttrString(group, "__contains__");
  } else {
    cerr << "group has to define __contains__" << endl;
    throw string(__FUNCTION__) + ": error.";
  }
}

is_element_general::~is_element_general() {
  Py_DECREF(method);
}

bool is_element_general::is_member(const SL2Z& m) const {
  PyObject* arg = convert_to_SL2Z(m);
  PyObject* tuple = PyTuple_New(1);
  PyTuple_SetItem(tuple, 0, arg);
  PyObject *result = PyObject_CallObject(method, tuple);
  Py_DECREF(tuple);
  if( not PyBool_Check(result) ) {
    cerr << "__contains__ does not return bool." << endl;
    throw string(__FUNCTION__) + ": error.";
  }
  bool value = (result == Py_True);
  Py_DECREF(result);
  return value;
}

//--- FareySymbol ------------------------------------------------------------

// SL2Z

FareySymbol::FareySymbol() {
  pairing = vector<int>(2);
  pairing[0] = EVEN;
  pairing[1] = ODD;
  pairing_max = NO;
  a.push_back(0);
  b.push_back(1);
  cusp_widths.push_back(1);
  coset.push_back(SL2Z::E);
  generators.push_back(SL2Z::S);
  generators.push_back(SL2Z::S*SL2Z::R);
  cusp_classes.push_back(0);
  even = true;
  pairing_in_group.push_back(true);
  pairing_in_group.push_back(true);
  for(size_t i=0; i<a.size(); i++) x.push_back(a[i]/b[i]);
}

FareySymbol::FareySymbol(istream& is) {
  is >> (*this);
}

// User defined group and restoring from pickle

FareySymbol::FareySymbol(PyObject* o) {
  if( PyBytes_Check(o) ) {
    // restoration from data
    istringstream is(PyBytes_AsString(o));
    is >> (*this);
  } else {
    // init with user defined group
    is_element_general *group = new is_element_general(o);
    // check for user defined SL2Z
    if( group->is_member(SL2Z::S) and
        group->is_member(SL2Z::T) ) {
      pairing = vector<int>(2);
      pairing[0] = EVEN;
      pairing[1] = ODD;
      pairing_max = NO;
      a.push_back(0);
      b.push_back(1);
      cusp_widths.push_back(1);
      coset.push_back(SL2Z::E);
      generators.push_back(SL2Z::S);
      generators.push_back(SL2Z::S*SL2Z::R);
      cusp_classes.push_back(0);
      even = true;
      pairing_in_group.push_back(true);
      pairing_in_group.push_back(true);
      x.push_back(a[0]/b[0]);
      reductions.push_back(SL2Z::E);
    }
    // check for index two subgroup
    else if( group->is_member(SL2Z( 0,  1, -1, -1)) and
             group->is_member(SL2Z(-1,  1, -1,  0)) ) {
      pairing = vector<int>(2);
      pairing[0] = ODD;
      pairing[1] = ODD;
      pairing_max = NO;
      a.push_back(0);
      b.push_back(1);
      x.push_back(a[0]/b[0]);
      coset.push_back(SL2Z::E);
      if ( group->is_member(SL2Z(0, -1, 1, 1)) ) {
        // index 2 even subgroup
        generators.push_back(SL2Z( 0, 1, -1, -1));
        generators.push_back(SL2Z(-1, 1, -1, 0));
        coset.push_back(SL2Z(0, 1, -1, 0));
      } else {
        // index 4 odd subgroup
        generators.push_back(SL2Z(0, 1, -1, -1));
        coset.push_back(SL2Z(0, -1,  1,  0));
        coset.push_back(SL2Z(1, -1,  1,  0));
        coset.push_back(SL2Z(1,  0, -1,  1));
        generators.push_back(SL2Z(-1, 1, -1, 0));
      }
      cusp_classes.push_back(0);
      reductions.push_back(SL2Z::E);
      if (group->is_member(SL2Z::I)) even = true;
      else even = false;
      pairing_in_group = init_sl2z_lift(group);
    } else {
      // everything else
      init_pairing(group);
      cusp_widths  = init_cusp_widths();
      coset        = init_coset_reps();
      generators   = init_generators(group);
      cusp_classes = init_cusp_classes();
      for(size_t i=0; i<a.size(); i++) x.push_back(a[i]/b[i]);
      cusps        = init_cusps();
      reductions   = init_reductions();
      if (group->is_member(SL2Z::I)) even = true;
      else even = false;
      pairing_in_group = init_sl2z_lift(group);
    }
    delete group;
  }
}

// Predefined subgroups of SL2Z

FareySymbol::FareySymbol(PyObject* o, const is_element_group* group) {
  init_pairing(group);
  cusp_widths  = init_cusp_widths();
  coset        = init_coset_reps();
  generators   = init_generators(group);
  cusp_classes = init_cusp_classes();
  for(size_t i=0; i<a.size(); i++) x.push_back(a[i]/b[i]);
  cusps        = init_cusps();
  reductions   = init_reductions();
  if (group->is_member(SL2Z::I)) even = true;
  else even = false;
  pairing_in_group = init_sl2z_lift(group);
}

FareySymbol::~FareySymbol() {
}

// for debugging purposes

void FareySymbol::dump(ostream& os) const {
  os << "Dumping FareySymbol:" << endl
     << tab << "pairing_max: " << pairing_max << endl
     << tab << "pairing: " << pairing << endl
     << tab << "a: " << a << endl
     << tab << "b: " << b << endl
     << tab << "x: " << x << endl
     << tab << "coset: " << coset << endl
     << tab << "generators: " << generators << endl
     << tab << "cusps: " << cusps << endl
     << tab << "cusp classes: " << cusp_classes << endl
     << tab << "cusp widths: " << cusp_widths << endl
     << tab << "reductions: " << reductions << endl;
}

void FareySymbol::init_pairing(const is_element_group* group) {
  pairing = vector<int>(3, NO);
  const mpq_class infinity(10000000);
  pairing_max = NO;
  if( group->is_member(SL2Z(-1, 1, -1, 0)) ) {
    a.push_back(-1);
    a.push_back(0);
    b.push_back(1);
    b.push_back(1);
  } else {
    a.push_back(0);
    a.push_back(1);
    b.push_back(1);
    b.push_back(1);
  }
  check_pair(group, 0);
  check_pair(group, 1);
  check_pair(group, 2);
  for(;;) {
    int missing_pair(-1);
    mpq_class largest_diameter(0);
    for(size_t i=0; i<pairing.size(); i++) {
      if( pairing[i] == NO ) {
        if( i+1 != pairing.size() ) {
          if( i != 0 ) {
            mpq_class d = a[i]/b[i] - a[i-1]/b[i-1];
            if( d > largest_diameter ) {
              largest_diameter = d;
              missing_pair = (int)(i);
            }
          } else {
            largest_diameter = infinity;
            missing_pair = 0;
          }
        } else {
          largest_diameter = infinity;
          missing_pair = (int)(pairing.size()-1);
          break;
        }
      }
    }
    if( missing_pair == -1 ) {
      break;
    } else {
      mpz_class A, B;
      if( missing_pair+1 == pairing.size() ) {
        A = a[missing_pair-1] + 1;
        B = b[missing_pair-1] + 0;
      } else {
        if( missing_pair == 0 ) {
          A = a[0] - 1;
          B = b[0] + 0;
        } else {
          A = a[missing_pair-1]+a[missing_pair];
          B = b[missing_pair-1]+b[missing_pair];
        }
      }
      add_term(missing_pair, A/B);
    }
    check_pair(group, missing_pair);
    check_pair(group, missing_pair+1);
  }
}

void FareySymbol::add_term(const int i, const mpq_class& r) {
  a.insert(a.begin()+i, r.get_num());
  b.insert(b.begin()+i, r.get_den());
  pairing.insert(pairing.begin()+i, NO);
}

void FareySymbol::check_pair(const is_element_group* group, const int i) {
  if( pairing[i] == NO ) {
    vector<int> even(pairing), odd(pairing);
    even[i] = EVEN;
    odd [i] = ODD;
    SL2Z A = pairing_matrix(even, i);
    SL2Z B = pairing_matrix(odd , i);
    if( group->is_member(A) or group->is_member(-A) ) {
      pairing[i] = EVEN;
      return;
    } else if( group->is_member(B) or group->is_member(-B) ) {
      pairing[i] = ODD;
      return;
    }
  }
  if( pairing[i] == NO ) {
    for(size_t j=0; j<pairing.size(); j++) {
      if( pairing[j] == NO and i != j ) {
        vector<int> p(pairing);
        p[i] = pairing_max+1;
        p[j] = pairing_max+1;
        SL2Z C = pairing_matrix(p, i);
        if( group->is_member(C) or group->is_member(-C) ) {
          pairing_max++;
          pairing[i] = pairing_max;
          pairing[j] = pairing_max;
          return;
        }
      }
    }
  }
}

SL2Z FareySymbol::pairing_matrix(const vector<int>& p, const size_t i) const {
  mpz_class ai, ai1, bi, bi1, aj, aj1, bj, bj1;
  if( i == 0 ) {
    ai = -1; bi = 0; ai1 = a[0], bi1 = b[0];
  } else if( i+1 == p.size() ) {
    ai = a[i-1]; bi = b[i-1]; ai1 = 1; bi1 = 0;
  } else {
    ai = a[i-1]; bi = b[i-1]; ai1 = a[i]; bi1 = b[i];
  }
  if( p[i] == NO ) {
    throw string(__FUNCTION__)+string(": error");
  } else if( p[i] == EVEN ) {
    return SL2Z(ai1*bi1+ai*bi, -ai*ai-ai1*ai1,
                bi*bi+bi1*bi1, -ai1*bi1-ai*bi);
  } else if( p[i] == ODD ) {
    return SL2Z(ai1*bi1+ai*bi1+ai*bi, -ai*ai-ai*ai1-ai1*ai1,
                bi*bi+bi*bi1+bi1*bi1, -ai1*bi1-ai1*bi-ai*bi);
  } else if( p[i] > NO ) {
    const size_t j = paired_side(p, i);
    if( j == 0 ) {
      aj = -1; bj = 0; aj1 = a[0]; bj1 = b[0];
    } else if( j == a.size() ) {
      aj = a[j-1]; bj = b[j-1]; aj1 = 1; bj1 = 0;
    } else {
      aj = a[j-1]; bj = b[j-1]; aj1 = a[j]; bj1 = b[j];
    }
    return SL2Z(aj1*bi1+aj*bi, -aj*ai-aj1*ai1,
                bj*bi+bj1*bi1, -ai1*bj1-ai*bj);
  }
  return SL2Z::E;
}

inline
SL2Z FareySymbol::pairing_matrix(const size_t i) const {
  return pairing_matrix(pairing, i);
}

SL2Z FareySymbol::pairing_matrix_in_group(const size_t i) const {
  if (pairing_in_group[i])
    return pairing_matrix(i);
  else
    return SL2Z::I*pairing_matrix(i);
}

size_t FareySymbol::paired_side(const vector<int>& p, const size_t n) const {
  if( p[n] == EVEN or p[n] == ODD ) {
    return n;
  } else if( p[n] > NO ) {
    vector<int>::const_iterator i = find(p.begin(), p.end(), p[n]);
    if( i-p.begin() != n ) {
      return i-p.begin();
    } else {
      vector<int>::const_iterator j = find(i+1, p.end(), p[n]);
      return j-p.begin();
    }
  }
  throw string(__FUNCTION__)+string(": error");
  return 0;
}

vector<SL2Z> FareySymbol::init_generators(const is_element_group *group) const {
  const SL2Z I(-1, 0, 0, -1);
  vector<SL2Z> gen;
  vector<int> p;
  for(size_t i=0; i<pairing.size(); i++) {
    if( find(p.begin(), p.end(), pairing[i]) == p.end() ) {
      SL2Z m = pairing_matrix(i);
      if( not group->is_member(m) ) m = I*m;
      if( pairing[i] == ODD and group->is_member(I) ) m = I*m;
      gen.push_back(m);
      if( pairing[i] > NO ) p.push_back(pairing[i]);
    }
  }
  if( nu2() == 0 and nu3() == 0 and group->is_member(I) ) gen.push_back(I);
  return gen;
}

vector<mpq_class> FareySymbol::init_cusp_widths() const {
  static const mpq_class one_half(1, 2);
  vector<mpz_class> A(a), B(b);
  A.push_back(1);
  B.push_back(0);
  vector<mpq_class> w(A.size(), 0);
  for(size_t i=0; i<w.size(); i++) {
    size_t im = (i==0 ? A.size()-1 : i-1);
    size_t ip = (i+1==A.size() ? 0 : i+1);
    w[i] = abs(A[im]*B[ip]-A[ip]*B[im]);
    if( pairing[i ] == ODD ) w[i] += one_half;
    if( pairing[ip] == ODD ) w[i] += one_half;
  }
  return w;
}

vector<SL2Z> FareySymbol::init_coset_reps() const {
  static const mpq_class one_half(1, 2);
  vector<mpz_class> A(a), B(b);
  A.insert(A.begin(), -1);
  B.insert(B.begin(),  0);
  vector<mpq_class> cw(cusp_widths);
  rotate(cw.begin(), cw.end()-1, cw.end());
  vector<int> p(pairing);
  rotate(p.begin(), p.end()-1, p.end());
  vector<SL2Z> reps;
  for(size_t i=0; i<p.size(); i++) {
    size_t j = (i+1) % p.size();
    mpz_class c(A[j]), d(B[j]);
    if( d == 0 ) c = 1;
    mpq_class upper_bound(cw[i]);
    if( p[i] == ODD ) upper_bound += one_half;
    if( p[j] == ODD ) upper_bound -= one_half;
    for(size_t k=0; k<upper_bound; k++) {
      reps.push_back(SL2Z(1, -(int)(k), 0, 1)/SL2Z(-A[i], c, -B[i], d));
    }
  }
  return reps;
}

vector<bool> FareySymbol::init_sl2z_lift(const is_element_group *group) const {
  vector<bool> result;
  for (size_t i=0; i<pairing.size(); i++)
    if (group->is_member(pairing_matrix(i)) )
      result.push_back(true);
    else
      result.push_back(false);
  return result;
}

// init cusp classes is a class of the pairing alone !!!

vector<int> FareySymbol::init_cusp_classes() const {
  vector<int> c(pairing.size(), 0);
  int cusp_number(1);
  for(size_t m=0; m<c.size(); m++) {
    if( c[m] != 0 ) {
      continue;
    }
    c[m] = cusp_number;
    size_t i(m), I, J;
    for(;;) {
      if( pairing[i] == NO ) {
        I = i;
        J = (i==0? pairing.size()-1 : (i-1)%c.size());
      } else {
        I = (i+1)%c.size();
        J = I;
      }
      if( pairing[I] == ODD or pairing[I] == EVEN ) {
        if( c[I] == cusp_number ) {
          cusp_number++;
          break;
        }
        c[J] = cusp_number;
        i = J;
        continue;
      } else if( pairing[I] > NO ) {
        size_t j;
        for(size_t k=0; k<c.size(); k++) {
          if( pairing[k] == pairing[I] and k != I ) j = k;
        }
        if( I != i ) {
          if( c[j] == cusp_number ) {
            cusp_number++;
            break;
          }
          c[j] = cusp_number;
          i = j;
          continue;
        } else {
          if( c[j-1] == cusp_number ) {
            cusp_number++;
            break;
          }
          c[j-1] = cusp_number;
          i = j-1;
          continue;
        }
      }
    }
  }
  for(size_t j=0; j<c.size(); j++) c[j]--;
  return c;
}

vector<mpq_class> FareySymbol::init_cusps() const {
  // initialize cusps by identifying fractions using the cusp_class number
  vector<mpq_class> c;
  for(int i=0; i<number_of_cusps(); i++) {
    for(size_t j=0; j<cusp_classes.size(); j++) {
      if( cusp_classes[j] == i and cusp_classes[j] != cusp_classes.back() ) {
        c.push_back(x[j]);
        break;
      }
    }
  }
  // in earlier version: shift negative cusps to positive ones
  sort(c.begin(), c.end());
  return c;
}

vector<SL2Z> FareySymbol::init_reductions() const {
  vector<SL2Z> reduction(x.size(), SL2Z::E);
  for(size_t j=0; j<x.size(); j++) {
    if( binary_search(cusps.begin(), cusps.end(), x[j]) ) continue;
    size_t k = j;
    mpq_class q = x[j];
    for(;;) {
      SL2Z p = pairing_matrix(k);
      reduction[j] = p*reduction[j];
      if( p.c()*q + p.d() == 0 ) break; // cusp ~ infinity
      q = p*q;
      if( binary_search(cusps.begin(), cusps.end(), q) ) break;
      k = lower_bound(x.begin(), x.end(), q) - x.begin();
    }
  }
  return reduction;
}

size_t FareySymbol::nu2() const {
  return count(pairing.begin(), pairing.end(), EVEN);
}

size_t FareySymbol::nu3() const {
  return count(pairing.begin(), pairing.end(), ODD);
}

size_t FareySymbol::rank_pi() const {
  using namespace std::placeholders;
  if( index() == 2 ) return 1;
  return count_if(pairing.begin(), pairing.end(),
                  bind(greater<int>(), placeholders::_1, 0))/2;
}

size_t FareySymbol::number_of_cusps() const {
  return size_t(*max_element(cusp_classes.begin(), cusp_classes.end()))+1;
}

size_t FareySymbol::genus() const {
  return (rank_pi()-number_of_cusps()+1)/2;
}

size_t FareySymbol::level() const {
  if( index() == 1 ) return 1;
  if( index() == 2 ) return 2;
  vector<mpz_class> A(a), B(b);
  A.push_back(1);
  B.push_back(0);
  vector<mpz_class> width;
  for(size_t i=0; i<number_of_cusps(); i++) {
    mpq_class cusp_width(0);
    for(size_t j=0; j<cusp_widths.size(); j++) {
      if( cusp_classes[j] == i ) {
        cusp_width += cusp_widths[j];
      }
    }
    width.push_back(cusp_width.get_num());
  }
  return lcm(width).get_ui();
}

long FareySymbol::side_index(const mpz_class& a0, const mpz_class& b0,
                             const mpz_class& a1, const mpz_class& b1) const {
  if( b0 == 0) {
    if( ( a1 == a[0] and b1 == b[0]) or
        (-a1 == a[0] and -b1 == b[0])) return 0;
  } else if( b1 == 0 ) {
        if( ( a0 == a.back() and  b0 == b.back()) or
            (-a0 == a.back() and -b0 == b.back())) return a.size();
  } else {
    mpq_class p = a1/b1, q = a0/b0;
    for(size_t j=1; j<a.size(); j++) {
      if( x[j-1] == q and x[j] == p ) return j;
    }
  }
  return -1;
}

void
FareySymbol::LLT_algorithm(const SL2Z& M, vector<int>& p, SL2Z& beta) const {
  // Based on: Kurth/Long - Computations with finite index subgroups
  // of PSL_2(ZZ) using Farey Symbols, p.13f
  beta = M;
  p.clear();
  bool found = false;
  mpq_class m;
  for(;;) {
    size_t k;
    found = false;
    mpz_class A=beta.a(), B=beta.b(), C=beta.c(), D=beta.d();
    if( D == 0 ) {
      if( A/C < x[0] ) {
        found = true;
        k = 0;
      } else if( x.back() < A/C ) {
        found = true;
        k = pairing.size()-1;
      }
    } else if( C == 0 ) {
      if( x.back() < B/D ) {
        found = true;
        k = pairing.size()-1;
      } else if( B/D < x[0] ) {
        found = true;
        k = 0;
      }
    } else if( A/C <= x[0] and B/D <= x[0] ) {
      found = true;
      k = 0;
    } else if( x.back() <= B/D and x.back() <= A/C ) {
      found = true;
      k = pairing.size()-1;
    } else {
      for(size_t i=0; i+1<x.size(); i++) {
        if( (x[i] <  B/D and B/D < A/C and A/C <= x[i+1]) or
            (x[i] <= B/D and B/D < A/C and A/C <  x[i+1]) or
            (x[i] <  A/C and A/C < B/D and B/D <= x[i+1]) or
            (x[i] <= A/C and A/C < B/D and B/D <  x[i+1]) ) {
          found = true;
          k = i+1;
          break;
        }
      }
    }
    if( not found ) break;
    // If the pairing, we found is ODD, we need to split the interval.
    // If the arc is in the right part, we choose the pairing matrix
    // If it is in the left part, we choose its inverse.
    // Note, that this choice differs from the article.
    if( pairing[k] == ODD ) {
      if( k == 0 ) {
        m = mpq_class(a[0]-1, b[0]);
      } else if ( k == pairing.size()-1 ) {
        m = mpq_class(a[0]+1, b[0]);
      } else {
        m = mpq_class(a[k-1]+a[k],b[k-1]+b[k]);
      }
      if ( C == 0 and B/D <=m) {
        beta = pairing_matrix_in_group(k).inverse()*beta;
        p.push_back(-(int)k - 1);
      } else if( C == 0 and B/D >= m) {
        beta = pairing_matrix_in_group(k) * beta;
        p.push_back((int)k + 1);
      } else if( D == 0 and A/C <= m) {
        beta = pairing_matrix_in_group(k).inverse()*beta;
        p.push_back(-(int)k - 1);
      } else if( D == 0 and A/C >=m) {
        beta = pairing_matrix_in_group(k)*beta;
        p.push_back((int)k + 1);
      } else if(A/C <= m and B/D <= m) {
        beta = pairing_matrix_in_group(k).inverse()*beta;
        p.push_back(-(int)k - 1);
      } else if(A/C >= m and B/D >= m) {
        beta = pairing_matrix_in_group(k)*beta;
        p.push_back((int)k + 1);
      } else {
        //Based on Lemma 4 of the article by Kurth/Long,
        //this case can not happen.
        throw string("Mathematical complications in ") +
              __FUNCTION__;
	return;
      }
    } else { // case of EVEN or FREE pairing
      beta = pairing_matrix_in_group(k)*beta;
      p.push_back((int)k+1);
    }
  }
}

bool FareySymbol::is_element(const SL2Z& M) const {
  // based on the same article as LLT_algorithm:
  // Kurth/Long - Computations with finite index ...
  vector<int> p;
  SL2Z beta = SL2Z::E;
  size_t k = 0;
  mpq_class smaller;
  LLT_algorithm(M, p, beta);

  // Case (1) of the article (p14)
  if( even ) {
    if( beta == SL2Z::E or beta == SL2Z::I ) return true;
  } else if( beta == SL2Z::E ) return true;

  // Case (3) of the article
  if( beta == SL2Z::S or beta == SL2Z::U ) {
    // If 0 and infty are adjacent vertices, they only can appear
    // at the beginning or the end.
    if( x[0] == 0 and pairing[0] == EVEN ) return true;
    if( x.back() == 0 and pairing.back() == EVEN ) return true;
  }

  // Case (2) of the article
  if( beta.c() != 0 and beta.d() != 0 ) {
    // Find index of the "edge beta"
    smaller = !((beta.b()/beta.d())<(beta.a()/beta.c())) ?
      (beta.a()/beta.c()) : (beta.b()/beta.d());
    for(size_t i=0; i<x.size(); i++)
      if( x[i] == smaller ) {
	k = i;
	break;
      }
    // Is it a free pairing?
      if( pairing[k+1] > NO ) {
	size_t s = paired_side(pairing, k+1);
	if ( s == 0 and x[0] == 0 and beta.a()/beta.c() > beta.b()/beta.d() )
	  // paired with (0,infty) at the beginning?
	  if( even ) {
	    return true;
	  } else {
	    beta = pairing_matrix_in_group(k)*beta;

	    if ( beta == SL2Z::E ) return true;
	    else return false;
	  }

	if( s == pairing.size() -1 and
	    x.back() == 0  and
	    beta.a()/beta.c() > beta.b()/beta.d() ) {
	  // paired with (0,infty) at the end?
	  return true;
	}
      }
  } else if( beta.c() == 0 and pairing.back() > NO ) {
    // In this case (a,infty) (for some a>0) is paired with (0,infty)
    size_t s = paired_side(pairing, pairing.size()-1);
    if( s == 0 and x.back() == beta.b()/beta.d() ) {
      if( even ) {
	return true;
      } else {
	beta = pairing_matrix_in_group(pairing.size()-1)*beta;
	if( beta == SL2Z::E ) return true;
	else return false;
      }
    }
  }
  return false;
}

SL2Z FareySymbol::reduce_to_fraction(const mpq_class& q) const {
  SL2Z M = rademacher_matrix(q);
  for(size_t j=0; j<coset.size(); j++) {
    SL2Z T = coset[j].inverse()*M;
    if( is_element(T) ) return T;
  }
  return SL2Z::E;
}

SL2Z FareySymbol::reduce_to_elementary_cusp(const mpq_class& q) const {
  SL2Z M = reduce_to_fraction(q);
  if( M.c()*q+M.d() == 0 ) return M;
  mpq_class Q = M*q;
  vector<mpq_class>::const_iterator k = find(x.begin(), x.end(), Q);
  if( k != x.end() ) {
    return reductions[k-x.begin()]*M;
  } else {
    return M;
  }
}

size_t FareySymbol::cusp_class(const mpq_class& q) const {
  typedef vector<int>::const_iterator const_iterator;
  SL2Z M = FareySymbol::reduce_to_elementary_cusp(q);
  if( M.c()*q + M.d() == 0 ) return cusp_classes.back();
  mpq_class Q = M*q;
  size_t k = lower_bound(x.begin(), x.end(), Q) - x.begin();
  return cusp_classes[k];
}

//--- communication with sage ------------------------------------------------

PyObject* FareySymbol::get_transformation_to_cusp(const mpz_t a,
						  const mpz_t b) const {
  mpz_class p(a), q(b);
  if( q == 0 ) return convert_to_SL2Z(SL2Z::E);
  mpq_class r(p, q);
  r.canonicalize();
  PyObject* M = convert_to_SL2Z(reduce_to_elementary_cusp(r));
  return M;
}

size_t FareySymbol::index() const {
  return coset.size();
}

PyObject* FareySymbol::is_element(const mpz_t a, const mpz_t b,
				  const mpz_t c, const mpz_t d) const {
  const SL2Z M = SL2Z(mpz_class(a), mpz_class(b), mpz_class(c), mpz_class(d));
  if( is_element(M) ) {
    Py_RETURN_TRUE;
  } else {
    Py_RETURN_FALSE;
  }
}

PyObject* FareySymbol::word_problem(const mpz_t a, const mpz_t b,
				    const mpz_t c, const mpz_t d, SL2Z *beta) const {
  const SL2Z M = SL2Z(mpz_class(a), mpz_class(b), mpz_class(c), mpz_class(d));
  vector<int> p;
  PyObject* wd;
  SL2Z beta1 = SL2Z::E;
  size_t i;
  LLT_algorithm(M, p, beta1);
  wd = PyList_New(p.size());
  for(i=0; i<p.size(); i++) {
    PyList_SetItem(wd, i, PyLong_FromLong(p[i]));
  }
  *beta = beta1;
  return wd;
}


PyObject* FareySymbol::get_coset() const {
  PyObject* coset_list = PyList_New(coset.size());
  for(size_t i=0; i<coset.size(); i++) {
    PyObject* m = convert_to_SL2Z(coset[i]);
    PyList_SetItem(coset_list, i, m);
  }
  return coset_list;
}

PyObject* FareySymbol::get_generators() const {
  PyObject* generators_list = PyList_New(generators.size());
  for(size_t i=0; i<generators.size(); i++) {
    PyObject* m = convert_to_SL2Z(generators[i]);
    PyList_SetItem(generators_list, i, m);
  }
  return generators_list;
}

PyObject* FareySymbol::get_cusps() const {
  PyObject* cusps_list = PyList_New(cusps.size());
  for(size_t i=0; i<cusps.size(); i++) {
    PyObject* m = convert_to_cusp(cusps[i]);
    PyList_SetItem(cusps_list, i, m);
  }
  return cusps_list;
}

PyObject* FareySymbol::get_cusp_widths() const {
  vector<mpz_class> width;
  for(size_t i=0; i<number_of_cusps(); i++) {
    mpq_class cusp_width(0);
    for(size_t j=0; j<cusp_widths.size(); j++) {
      if( cusp_classes[j] == i ) {
        cusp_width += cusp_widths[j];
      }
    }
    width.push_back(cusp_width.get_num());
  }
  PyObject* cusp_widths_list = PyList_New(width.size());
  for(size_t i=0; i<width.size(); i++) {
    PyObject* m = convert_to_rational(width[i]);
    PyList_SetItem(cusp_widths_list, i, m);
  }
  return cusp_widths_list;
}


size_t
FareySymbol::get_cusp_class(const mpz_t p, const mpz_t q) const {
  mpz_class a(p), b(q);
  if( a != 0 and b == 0 ) return cusp_classes.back();
  mpq_class c(a, b);
  c.canonicalize();
  return cusp_class(c);
}

PyObject* FareySymbol::get_fractions() const {
  PyObject* x_list = PyList_New(x.size());
  for(size_t i=0; i<x.size(); i++) {
    PyObject* m = convert_to_rational(x[i]);
    PyList_SetItem(x_list, i, m);
  }
  return x_list;
}

PyObject* FareySymbol::get_pairings() const {
  PyObject* pairing_list = PyList_New(pairing.size());
  for(size_t i=0; i<pairing.size(); i++) {
    PyObject* m = PyLong_FromLong(long(pairing[i]));
    PyList_SetItem(pairing_list, i, m);
  }
  return pairing_list;
}

PyObject* FareySymbol::get_paired_sides() const {
  vector<int> p;
  for(size_t i=0; i<pairing.size(); i++) {
    if( pairing[i] > NO and
        p.end() == find(p.begin(), p.end(), pairing[i]) ) {
      p.push_back(pairing[i]);
    }
  }
  sort(p.begin(), p.end());
  PyObject* pairing_list = PyList_New(p.size());
  for(vector<int>::const_iterator i=p.begin(); i!=p.end(); i++) {
    vector<int>::const_iterator j = find(pairing.begin(), pairing.end(), *i);
    vector<int>::const_iterator k = find(j+1, pairing.end(), *i);
    PyObject* J = PyLong_FromLong(long(j-pairing.begin()));
    PyObject* K = PyLong_FromLong(long(k-pairing.begin()));
    PyObject* tuple = PyTuple_New(2);
    PyTuple_SetItem(tuple, 0, J);
    PyTuple_SetItem(tuple, 1, K);
    PyList_SetItem(pairing_list, i-p.begin(), tuple);
  }
  return pairing_list;
}

PyObject* FareySymbol::get_pairing_matrices() const {
  PyObject* pairing_matrix_list = PyList_New(pairing.size());
  for(size_t i=0; i<pairing.size(); i++) {
    PyObject* pm = convert_to_SL2Z(pairing_matrix(i));
    PyList_SetItem(pairing_matrix_list, i, pm);
  }
  return pairing_matrix_list;
}

PyObject* FareySymbol::dumps() const {
  std::ostringstream os(ostringstream::out|ostringstream::binary);
  os << (*this);
  PyObject* d = PyBytes_FromString(os.str().c_str());
  return d;
}

