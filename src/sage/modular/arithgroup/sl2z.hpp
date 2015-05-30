//****************************************************************************
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
//****************************************************************************

#ifndef SL2Z_HPP_
#define SL2Z_HPP_

#include <iostream>
#include <iomanip>
#include <cassert>
#include <gmpxx.h>

class SL2Z {
public:
  typedef mpz_class ElementType;
protected:
  ElementType M[2][2];
public:
  const static SL2Z E, R, T, S, U, I;
  SL2Z(int a_, int b_, int c_, int d_);
  SL2Z(const ElementType& a_, const ElementType& b_,
       const ElementType& c_, const ElementType& d_);
  SL2Z(const SL2Z&);
  ElementType a() const { return M[0][0]; };
  ElementType b() const { return M[0][1]; };
  ElementType c() const { return M[1][0]; };
  ElementType d() const { return M[1][1]; };
  SL2Z inverse() const;
  SL2Z operator-() const;
  SL2Z operator*=(const SL2Z& x);
  SL2Z operator/=(const SL2Z& x);
  SL2Z mod(const size_t) const;
  friend bool operator==(const SL2Z&, const SL2Z&);
  friend bool operator!=(const SL2Z&, const SL2Z&);
  friend SL2Z operator*(const SL2Z&, const SL2Z&);
  friend SL2Z operator/(const SL2Z&, const SL2Z&);
  friend std::ostream& operator<<(std::ostream&, const SL2Z&);
  friend std::istream& operator>>(std::istream&,       SL2Z&);
};

inline
SL2Z::SL2Z(const SL2Z::ElementType& a_, const SL2Z::ElementType& b_,
           const SL2Z::ElementType& c_, const SL2Z::ElementType& d_) {
  M[0][0] = a_;
  M[0][1] = b_;
  M[1][0] = c_;
  M[1][1] = d_;
  assert(M[0][0]*M[1][1] - M[0][1]*M[1][0] == 1);
}

inline
SL2Z::SL2Z(int a_, int b_, int c_, int d_) {
  M[0][0] = a_;
  M[0][1] = b_;
  M[1][0] = c_;
  M[1][1] = d_;
  assert(M[0][0]*M[1][1] - M[0][1]*M[1][0] == 1);
}

inline
SL2Z::SL2Z(const SL2Z& x) {
  M[0][0] = x.M[0][0];
  M[0][1] = x.M[0][1];
  M[1][0] = x.M[1][0];
  M[1][1] = x.M[1][1];
}

inline
SL2Z SL2Z::operator-() const {
  return SL2Z(-M[0][0],-M[0][1],-M[1][0],-M[1][1]);
}

inline
SL2Z SL2Z::operator*=(const SL2Z& x) {
  SL2Z result(M[0][0]*x.M[0][0] + M[0][1]*x.M[1][0],
	      M[0][0]*x.M[0][1] + M[0][1]*x.M[1][1],
	      M[1][0]*x.M[0][0] + M[1][1]*x.M[1][0],
	      M[1][0]*x.M[0][1] + M[1][1]*x.M[1][1]);
  M[0][0] = result.M[0][0];
  M[0][1] = result.M[0][1];
  M[1][0] = result.M[1][0];
  M[1][1] = result.M[1][1];
  return *this;
}

inline
SL2Z SL2Z::operator/=(const SL2Z& x) {
  SL2Z result( M[0][0]*x.M[1][1] - M[0][1]*x.M[1][0],
	      -M[0][0]*x.M[0][1] + M[0][1]*x.M[0][0],
	       M[1][0]*x.M[1][1] - M[1][1]*x.M[1][0],
	      -M[1][0]*x.M[0][1] + M[1][1]*x.M[0][0]);
  M[0][0] = result.M[0][0];
  M[0][1] = result.M[0][1];
  M[1][0] = result.M[1][0];
  M[1][1] = result.M[1][1];
  return *this;
}

inline
SL2Z SL2Z::inverse() const {
  return SL2Z(M[1][1], -M[0][1], -M[1][0], M[0][0]);
}

inline
SL2Z SL2Z::mod(const size_t n) const {
  return SL2Z(M[0][0]%n, M[0][1]%n, M[1][0]%n, M[1][1]%n);
}

inline
SL2Z operator*(const SL2Z& x, const SL2Z& y) {
  return SL2Z(x.M[0][0]*y.M[0][0] + x.M[0][1]*y.M[1][0],
              x.M[0][0]*y.M[0][1] + x.M[0][1]*y.M[1][1],
              x.M[1][0]*y.M[0][0] + x.M[1][1]*y.M[1][0],
              x.M[1][0]*y.M[0][1] + x.M[1][1]*y.M[1][1]);
}

inline
SL2Z operator/(const SL2Z& x, const SL2Z& y) {
  return SL2Z( x.M[0][0]*y.M[1][1] - x.M[0][1]*y.M[1][0],
              -x.M[0][0]*y.M[0][1] + x.M[0][1]*y.M[0][0],
               x.M[1][0]*y.M[1][1] - x.M[1][1]*y.M[1][0],
              -x.M[1][0]*y.M[0][1] + x.M[1][1]*y.M[0][0]);
}

inline
bool operator==(const SL2Z& x, const SL2Z& y) {
  return (x.M[0][0] == y.M[0][0] and
          x.M[0][1] == y.M[0][1] and
          x.M[1][0] == y.M[1][0] and
          x.M[1][1] == y.M[1][1]);
}

inline
bool operator!=(const SL2Z& x, const SL2Z& y) {
  return not (x == y);
}

inline
std::ostream& operator<<(std::ostream& os, const SL2Z& x) {
  os << "["
     << x.M[0][0] << ", "
     << x.M[0][1] << "; "
     << x.M[1][0] << ", "
     << x.M[1][1]
     << "]";
  return os;
}

inline
std::istream& operator>>(std::istream& is, SL2Z& x) {
  char c;
  is >> c;
  if( c == '[' ) {
    is >> x.M[0][0] >> c;
    if( c == ',' ) {
      is >> x.M[0][1] >> c;
      if( c == ';' ) {
        is >> x.M[1][0] >> c;
        if( c == ',' ) {
          is >> x.M[1][1] >> c;
          if( c != ']' ) is.clear(std::ios_base::badbit);
        } else {
          is.clear(std::ios_base::badbit);
        }
      } else {
        is.clear(std::ios_base::badbit);
      }
    } else {
      is.clear(std::ios_base::badbit);
    }
  } else {
    is.clear(std::ios_base::badbit);
  }
  return is;
}

#endif // SL2Z_HPP_
