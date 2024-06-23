// sage_setup: distribution = sagemath-schemes
//
//  farey.hpp
//  FareySymbol
//
//*************************************************************************
//       Copyright (C) 2011 Hartmut Monien <monien@th.physik.uni-bonn.de>
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
//*************************************************************************


#ifndef FAREY_SYMBOL_HPP_
#define FAREY_SYMBOL_HPP_

#include <iostream>
#include <vector>
#include <string>

#include <Python.h>
#include <gmpxx.h>
#include "sl2z.hpp"


//--- pure virtual base class for helper class for membership test --------

class is_element_group {
public:
  virtual bool is_member(const SL2Z&) const = 0;
};

class is_element_Gamma0 : public is_element_group {
  const int p;
public:
  is_element_Gamma0(int p_) : p(p_) {
  }
  bool is_member(const SL2Z& V) const {
    return (V.c() % p == 0);
  }
};

class is_element_Gamma1 : public is_element_group {
  const int p;
public:
  is_element_Gamma1(int p_) : p(p_) {
  }
  bool is_member(const SL2Z& V) const {
    return ((V.a()-1) % p == 0 &&
            V.c() % p == 0 &&
            (V.d()-1) % p == 0);
  }
};

class is_element_Gamma : public is_element_group {
  const int p;
public:
  is_element_Gamma(int p_) : p(p_) {
  }
  bool is_member(const SL2Z& V) const {
    return ((V.a()-1) % p == 0 &&
	    V.b() % p == 0 &&
	    V.c() % p == 0 &&
	    (V.d()-1) % p == 0);
  }
};

class is_element_GammaH : public is_element_group {
  const int p;
  std::vector<long> H;
public:
  is_element_GammaH(int p_, PyObject*);
  ~is_element_GammaH();
  bool is_member(const SL2Z&) const;
};

class is_element_general : public is_element_group {
protected:
  PyObject* group;
  PyObject* method;
public:
  is_element_general(PyObject*);
  virtual ~is_element_general();
  bool is_member(const SL2Z&) const;
};

class FareySymbol {
protected:
  enum PAIRING { EVEN=-2, ODD=-3, NO=0, FREE=1 };
  size_t pairing_max;
  std::vector<int> pairing;
  std::vector<int> cusp_classes;
  std::vector<mpz_class> a, b;
  std::vector<mpq_class> x;
  std::vector<SL2Z> coset;
  std::vector<SL2Z> generators;
  std::vector<mpq_class> cusps;
  std::vector<mpq_class> cusp_widths;
  std::vector<SL2Z> reductions;
  bool even;
  std::vector<bool> pairing_in_group; // For membership test:
  //Is the i-th pairing matrix in the group?
  size_t rank_pi() const;
  long side_index(const mpz_class& a0, const mpz_class& b0,
		  const mpz_class& a1, const mpz_class& b1) const;
  void LLT_algorithm(const SL2Z& M, std::vector<int>& p, SL2Z& beta) const;
  void dump(std::ostream& os) const;
  std::vector<SL2Z> init_reductions() const;
private:
  void add_term(const int, const mpq_class&);
  void check_pair(const is_element_group*, const int);
  size_t paired_side(const std::vector<int>& p, const size_t i) const;
  SL2Z pairing_matrix(const std::vector<int>&, const size_t i) const;
  void init_pairing(const is_element_group*);
  std::vector<SL2Z> init_generators(const is_element_group*) const;
  std::vector<SL2Z> init_coset_reps() const;
  std::vector<int> init_cusp_classes() const;
  std::vector<mpq_class> init_cusps() const;
  std::vector<mpq_class> init_cusp_widths() const;
  SL2Z pairing_matrix(const size_t) const;
  SL2Z pairing_matrix_in_group(const size_t) const;
  std::vector<bool> init_sl2z_lift(const is_element_group*) const;
public:
  FareySymbol();
  FareySymbol(std::istream& is);
  FareySymbol(PyObject*);
  FareySymbol(PyObject*, const is_element_group*);
  ~FareySymbol();
  const size_t size() const;
  size_t nu2() const;
  size_t nu3() const;
  size_t index() const;
  size_t number_of_cusps() const;
  size_t level() const;
  size_t genus() const;
  SL2Z reduce_to_fraction(const mpq_class& q) const;
  SL2Z reduce_to_elementary_cusp(const mpq_class& q) const;
  size_t cusp_class(const mpq_class& q) const;
  bool is_element(const SL2Z& M) const;
  friend std::ostream& operator<<(std::ostream&, const FareySymbol&);
  friend std::istream& operator>>(std::istream&,       FareySymbol&);
  //--- communication with sage -------------------------------------------
  PyObject* is_element(const mpz_t, const mpz_t, const mpz_t, const mpz_t) const;
  PyObject* word_problem(const mpz_t, const mpz_t, const mpz_t, const mpz_t, SL2Z *) const;
  PyObject* get_transformation_to_cusp(const mpz_t, const mpz_t) const;
  PyObject* get_cusps() const;
  PyObject* get_cusp_widths() const;
  size_t    get_cusp_class(const mpz_t, const mpz_t) const;
  PyObject* get_fractions() const;
  PyObject* get_coset() const;
  PyObject* get_generators() const;
  PyObject* get_pairings() const;
  PyObject* get_paired_sides() const;
  PyObject* get_pairing_matrices() const;
  PyObject* dumps() const;
};

#endif // FAREY_SYMBOL_HPP_
