/** @file matrix.cpp
 *
 *  Implementation of symbolic matrices */

/*
 *  GiNaC Copyright (C) 1999-2008 Johannes Gutenberg University Mainz, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "matrix.h"
#include "numeric.h"
#include "lst.h"
#include "add.h"
#include "power.h"
#include "symbol.h"
#include "operators.h"
#include "normal.h"
#include "archive.h"
#include "utils.h"
#include "upoly.h"

#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <map>
#include <stdexcept>

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(matrix, basic,
  print_func<print_context>(&matrix::do_print).
  print_func<print_latex>(&matrix::do_print_latex).
  print_func<print_tree>(&matrix::do_print_tree).
  print_func<print_python_repr>(&matrix::do_print_python_repr))

//////////
// default constructor
//////////

/** Default ctor.  Initializes to 1 x 1-dimensional zero-matrix. */
matrix::matrix() : inherited(&matrix::tinfo_static), row(1), col(1), m(1, _ex0)
{
	setflag(status_flags::not_shareable);
}

//////////
// other constructors
//////////

// public

/** Very common ctor.  Initializes to r x c-dimensional zero-matrix.
 *
 *  @param r number of rows
 *  @param c number of cols */
matrix::matrix(unsigned r, unsigned c)
  : inherited(&matrix::tinfo_static), row(r), col(c), m(r*c, _ex0)
{
	setflag(status_flags::not_shareable);
}

// protected

/** Ctor from representation, for internal use only. */
matrix::matrix(unsigned r, unsigned c, exvector  m2)
  : inherited(&matrix::tinfo_static), row(r), col(c), m(std::move(m2))
{
	setflag(status_flags::not_shareable);
}

/** Construct matrix from (flat) list of elements. If the list has fewer
 *  elements than the matrix, the remaining matrix elements are set to zero.
 *  If the list has more elements than the matrix, the excessive elements are
 *  thrown away. */
matrix::matrix(unsigned r, unsigned c, const lst & l)
  : inherited(&matrix::tinfo_static), row(r), col(c), m(r*c, _ex0)
{
	setflag(status_flags::not_shareable);

	size_t i = 0;
	for (auto it = l.begin(); it != l.end(); ++it, ++i) {
		size_t x = i % c;
		size_t y = i / c;
		if (y >= r)
			break; // matrix smaller than list: throw away excessive elements
		m[y*c+x] = *it;
	}
}

//////////
// archiving
//////////

matrix::matrix(const archive_node &n, lst &sym_lst) : inherited(n, sym_lst)
{
	setflag(status_flags::not_shareable);

	if (!(n.find_unsigned("row", row)) || !(n.find_unsigned("col", col)))
		throw (std::runtime_error("unknown matrix dimensions in archive"));
	m.reserve(row * col);
	auto first = n.find_first("m");
	auto last = n.find_last("m");
	++last;
	for (auto i=first; i<last; ++i) {
		ex e;
		n.find_ex_by_loc(i, e, sym_lst);
		m.push_back(e);
	}
}

void matrix::archive(archive_node &n) const
{
	inherited::archive(n);
	n.add_unsigned("row", row);
	n.add_unsigned("col", col);
	auto i = m.begin(), iend = m.end();
	while (i != iend) {
		n.add_ex("m", *i);
		++i;
	}
}

//////////
// functions overriding virtual functions from base classes
//////////

// public

void matrix::print_elements(const print_context & c, const char *row_start, const char *row_end, const char *row_sep, const char *col_sep) const
{
	for (unsigned ro=0; ro<row; ++ro) {
		c.s << row_start;
		for (unsigned co=0; co<col; ++co) {
			m[ro*col+co].print(c);
			if (co < col-1)
				c.s << col_sep;
			else
				c.s << row_end;
		}
		if (ro < row-1)
			c.s << row_sep;
	}
}

void matrix::do_print(const print_context & c, unsigned level) const
{
	c.s << "[";
	print_elements(c, "[", "]", ",", ",");
	c.s << "]";
}

void matrix::do_print_latex(const print_latex & c, unsigned level) const
{
	c.s << "\\left(\\begin{array}{" << std::string(col,'c') << "}";
	print_elements(c, "", "", "\\\\", "&");
	c.s << "\\end{array}\\right)";
}

void matrix::do_print_python_repr(const print_python_repr & c, unsigned level) const
{
	c.s << class_name() << '(';
	print_elements(c, "[", "]", ",", ",");
	c.s << ')';
}

/** nops is defined to be rows x columns. */
size_t matrix::nops() const
{
	return static_cast<size_t>(row) * static_cast<size_t>(col);
}

/** returns matrix entry at position (i/col, i%col). */
const ex matrix::op(size_t i) const
{
	GINAC_ASSERT(i<nops());
	
	return m[i];
}

/** returns writable matrix entry at position (i/col, i%col). */
ex & matrix::let_op(size_t i)
{
	GINAC_ASSERT(i<nops());
	
	ensure_if_modifiable();
	return m[i];
}

/** Evaluate matrix entry by entry. */
ex matrix::eval(int level) const
{
	// check if we have to do anything at all
	if (level == 1 and is_evaluated())
		return *this;
	
	// emergency break
	if (level == -max_recursion_level)
		throw (std::runtime_error("matrix::eval(): recursion limit exceeded"));
	
	// eval() entry by entry
	exvector m2(row*col);
	--level;
	for (unsigned r=0; r<row; ++r)
		for (unsigned c=0; c<col; ++c)
			m2[r*col+c] = m[r*col+c].eval(level);
	
	return (new matrix(row, col, m2))->setflag(status_flags::dynallocated |
	                                           status_flags::evaluated);
}

ex matrix::subs(const exmap & mp, unsigned options) const
{
	exvector m2(row * col);
	for (unsigned r=0; r<row; ++r)
		for (unsigned c=0; c<col; ++c)
			m2[r*col+c] = m[r*col+c].subs(mp, options);

	return matrix(row, col, m2).subs_one_level(mp, options);
}

/** Complex conjugate every matrix entry. */
ex matrix::conjugate() const
{
	std::unique_ptr<exvector> ev(nullptr);
	for (auto i=m.begin(); i!=m.end(); ++i) {
		ex x = i->conjugate();
		if (ev) {
			ev->push_back(x);
			continue;
		}
		if (are_ex_trivially_equal(x, *i)) {
			continue;
		}
		ev.reset(new exvector);
		ev->reserve(m.size());
		for (auto j=m.begin(); j!=i; ++j) {
			ev->push_back(*j);
		}
		ev->push_back(x);
	}
	if (ev) {
		ex result = matrix(row, col, *ev);
		return result;
	}
	return *this;
}

ex matrix::real_part() const
{
	exvector v;
	v.reserve(m.size());
	for (const auto & elem : m)
		v.push_back(elem.real_part());
	return matrix(row, col, v);
}

ex matrix::imag_part() const
{
	exvector v;
	v.reserve(m.size());
	for (const auto & elem : m)
		v.push_back(elem.imag_part());
	return matrix(row, col, v);
}

// protected

int matrix::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_exactly_a<matrix>(other));
	const matrix &o = static_cast<const matrix &>(other);
	
	// compare number of rows
	if (row != o.rows())
		return row < o.rows() ? -1 : 1;
	
	// compare number of columns
	if (col != o.cols())
		return col < o.cols() ? -1 : 1;
	
	// equal number of rows and columns, compare individual elements
	int cmpval;
	for (unsigned r=0; r<row; ++r) {
		for (unsigned c=0; c<col; ++c) {
			cmpval = ((*this)(r,c)).compare(o(r,c));
			if (cmpval!=0) return cmpval;
		}
	}
	// all elements are equal => matrices are equal;
	return 0;
}

bool matrix::match_same_type(const basic & other) const
{
	GINAC_ASSERT(is_exactly_a<matrix>(other));
	const matrix & o = static_cast<const matrix &>(other);
	
	// The number of rows and columns must be the same. This is necessary to
	// prevent a 2x3 matrix from matching a 3x2 one.
	return row == o.rows() && col == o.cols();
}

#if 0
/** Automatic symbolic evaluation of an indexed matrix. */
ex matrix::eval_indexed(const basic & i) const
{
	GINAC_ASSERT(is_a<indexed>(i));
	GINAC_ASSERT(is_a<matrix>(i.op(0)));

	bool all_indices_unsigned = static_cast<const indexed &>(i).all_index_values_are(info_flags::nonnegint);

	// Check indices
	if (i.nops() == 2) {

		// One index, must be one-dimensional vector
		if (row != 1 && col != 1)
			throw (std::runtime_error("matrix::eval_indexed(): vector must have exactly 1 index"));

		const idx & i1 = ex_to<idx>(i.op(1));

		if (col == 1) {

			// Column vector
			if (!i1.get_dim().is_equal(row))
				throw (std::runtime_error("matrix::eval_indexed(): dimension of index must match number of vector elements"));

			// Index numeric -> return vector element
			if (all_indices_unsigned) {
				unsigned n1 = ex_to<numeric>(i1.get_value()).to_int();
				if (n1 >= row)
					throw (std::runtime_error("matrix::eval_indexed(): value of index exceeds number of vector elements"));
				return (*this)(n1, 0);
			}

		} else {

			// Row vector
			if (!i1.get_dim().is_equal(col))
				throw (std::runtime_error("matrix::eval_indexed(): dimension of index must match number of vector elements"));

			// Index numeric -> return vector element
			if (all_indices_unsigned) {
				unsigned n1 = ex_to<numeric>(i1.get_value()).to_int();
				if (n1 >= col)
					throw (std::runtime_error("matrix::eval_indexed(): value of index exceeds number of vector elements"));
				return (*this)(0, n1);
			}
		}

	} else if (i.nops() == 3) {

		// Two indices
		const idx & i1 = ex_to<idx>(i.op(1));
		const idx & i2 = ex_to<idx>(i.op(2));

		if (!i1.get_dim().is_equal(row))
			throw (std::runtime_error("matrix::eval_indexed(): dimension of first index must match number of rows"));
		if (!i2.get_dim().is_equal(col))
			throw (std::runtime_error("matrix::eval_indexed(): dimension of second index must match number of columns"));

		// Pair of dummy indices -> compute trace
		if (is_dummy_pair(i1, i2))
			return trace();

		// Both indices numeric -> return matrix element
		if (all_indices_unsigned) {
			unsigned n1 = ex_to<numeric>(i1.get_value()).to_int(), n2 = ex_to<numeric>(i2.get_value()).to_int();
			if (n1 >= row)
				throw (std::runtime_error("matrix::eval_indexed(): value of first index exceeds number of rows"));
			if (n2 >= col)
				throw (std::runtime_error("matrix::eval_indexed(): value of second index exceeds number of columns"));
			return (*this)(n1, n2);
		}

	} else
		throw (std::runtime_error("matrix::eval_indexed(): matrix must have exactly 2 indices"));

	return i.hold();
}

/** Sum of two indexed matrices. */
ex matrix::add_indexed(const ex & self, const ex & other) const
{
	GINAC_ASSERT(is_a<indexed>(self));
	GINAC_ASSERT(is_a<matrix>(self.op(0)));
	GINAC_ASSERT(is_a<indexed>(other));
	GINAC_ASSERT(self.nops() == 2 || self.nops() == 3);

	// Only add two matrices
	if (is_a<matrix>(other.op(0))) {
		GINAC_ASSERT(other.nops() == 2 || other.nops() == 3);

		const matrix &self_matrix = ex_to<matrix>(self.op(0));
		const matrix &other_matrix = ex_to<matrix>(other.op(0));

		if (self.nops() == 2 && other.nops() == 2) { // vector + vector

			if (self_matrix.row == other_matrix.row)
				return indexed(self_matrix.add(other_matrix), self.op(1));
			else if (self_matrix.row == other_matrix.col)
				return indexed(self_matrix.add(other_matrix.transpose()), self.op(1));

		} else if (self.nops() == 3 && other.nops() == 3) { // matrix + matrix

			if (self.op(1).is_equal(other.op(1)) && self.op(2).is_equal(other.op(2)))
				return indexed(self_matrix.add(other_matrix), self.op(1), self.op(2));
			else if (self.op(1).is_equal(other.op(2)) && self.op(2).is_equal(other.op(1)))
				return indexed(self_matrix.add(other_matrix.transpose()), self.op(1), self.op(2));

		}
	}

	// Don't know what to do, return unevaluated sum
	return self + other;
}

/** Product of an indexed matrix with a number. */
ex matrix::scalar_mul_indexed(const ex & self, const numeric & other) const
{
	GINAC_ASSERT(is_a<indexed>(self));
	GINAC_ASSERT(is_a<matrix>(self.op(0)));
	GINAC_ASSERT(self.nops() == 2 || self.nops() == 3);

	const matrix &self_matrix = ex_to<matrix>(self.op(0));

	if (self.nops() == 2)
		return indexed(self_matrix.mul(other), self.op(1));
	else // self.nops() == 3
		return indexed(self_matrix.mul(other), self.op(1), self.op(2));
}

/** Contraction of an indexed matrix with something else. */
bool matrix::contract_with(exvector::iterator self, exvector::iterator other, exvector & v) const
{
	GINAC_ASSERT(is_a<indexed>(*self));
	GINAC_ASSERT(is_a<indexed>(*other));
	GINAC_ASSERT(self->nops() == 2 || self->nops() == 3);
	GINAC_ASSERT(is_a<matrix>(self->op(0)));

	// Only contract with other matrices
	if (!is_a<matrix>(other->op(0)))
		return false;

	GINAC_ASSERT(other->nops() == 2 || other->nops() == 3);

	const matrix &self_matrix = ex_to<matrix>(self->op(0));
	const matrix &other_matrix = ex_to<matrix>(other->op(0));

	if (self->nops() == 2) {

		if (other->nops() == 2) { // vector * vector (scalar product)

			if (self_matrix.col == 1) {
				if (other_matrix.col == 1) {
					// Column vector * column vector, transpose first vector
					*self = self_matrix.transpose().mul(other_matrix)(0, 0);
				} else {
					// Column vector * row vector, swap factors
					*self = other_matrix.mul(self_matrix)(0, 0);
				}
			} else {
				if (other_matrix.col == 1) {
					// Row vector * column vector, perfect
					*self = self_matrix.mul(other_matrix)(0, 0);
				} else {
					// Row vector * row vector, transpose second vector
					*self = self_matrix.mul(other_matrix.transpose())(0, 0);
				}
			}
			*other = _ex1;
			return true;

		} else { // vector * matrix

			// B_i * A_ij = (B*A)_j (B is row vector)
			if (is_dummy_pair(self->op(1), other->op(1))) {
				if (self_matrix.row == 1)
					*self = indexed(self_matrix.mul(other_matrix), other->op(2));
				else
					*self = indexed(self_matrix.transpose().mul(other_matrix), other->op(2));
				*other = _ex1;
				return true;
			}

			// B_j * A_ij = (A*B)_i (B is column vector)
			if (is_dummy_pair(self->op(1), other->op(2))) {
				if (self_matrix.col == 1)
					*self = indexed(other_matrix.mul(self_matrix), other->op(1));
				else
					*self = indexed(other_matrix.mul(self_matrix.transpose()), other->op(1));
				*other = _ex1;
				return true;
			}
		}

	} else if (other->nops() == 3) { // matrix * matrix

		// A_ij * B_jk = (A*B)_ik
		if (is_dummy_pair(self->op(2), other->op(1))) {
			*self = indexed(self_matrix.mul(other_matrix), self->op(1), other->op(2));
			*other = _ex1;
			return true;
		}

		// A_ij * B_kj = (A*Btrans)_ik
		if (is_dummy_pair(self->op(2), other->op(2))) {
			*self = indexed(self_matrix.mul(other_matrix.transpose()), self->op(1), other->op(1));
			*other = _ex1;
			return true;
		}

		// A_ji * B_jk = (Atrans*B)_ik
		if (is_dummy_pair(self->op(1), other->op(1))) {
			*self = indexed(self_matrix.transpose().mul(other_matrix), self->op(2), other->op(2));
			*other = _ex1;
			return true;
		}

		// A_ji * B_kj = (B*A)_ki
		if (is_dummy_pair(self->op(1), other->op(2))) {
			*self = indexed(other_matrix.mul(self_matrix), other->op(1), self->op(2));
			*other = _ex1;
			return true;
		}
	}

	return false;
}
#endif

//////////
// non-virtual functions in this class
//////////

// public

/** Sum of matrices.
 *
 *  @exception logic_error (incompatible matrices) */
matrix matrix::add(const matrix & other) const
{
	if (col != other.col || row != other.row)
		throw std::logic_error("matrix::add(): incompatible matrices");
	
	exvector sum(this->m);
	auto i = sum.begin(), end = sum.end();
	auto ci = other.m.begin();
	while (i != end)
		*i++ += *ci++;
	
	return matrix(row,col,sum);
}


/** Difference of matrices.
 *
 *  @exception logic_error (incompatible matrices) */
matrix matrix::sub(const matrix & other) const
{
	if (col != other.col || row != other.row)
		throw std::logic_error("matrix::sub(): incompatible matrices");
	
	exvector dif(this->m);
	auto i = dif.begin(), end = dif.end();
	auto ci = other.m.begin();
	while (i != end)
		*i++ -= *ci++;
	
	return matrix(row,col,dif);
}


/** Product of matrices.
 *
 *  @exception logic_error (incompatible matrices) */
matrix matrix::mul(const matrix & other) const
{
	if (this->cols() != other.rows())
		throw std::logic_error("matrix::mul(): incompatible matrices");
	
	exvector prod(this->rows()*other.cols());
	
	for (unsigned r1=0; r1<this->rows(); ++r1) {
		for (unsigned c=0; c<this->cols(); ++c) {
			// Quick test: can we shortcut?
			if (m[r1*col+c].is_zero())
				continue;
			for (unsigned r2=0; r2<other.cols(); ++r2)
				prod[r1*other.col+r2] += (m[r1*col+c] * other.m[c*other.col+r2]);
		}
	}
	return matrix(row, other.col, prod);
}


/** Product of matrix and scalar. */
matrix matrix::mul(const numeric & other) const
{
	exvector prod(row * col);

	for (unsigned r=0; r<row; ++r)
		for (unsigned c=0; c<col; ++c)
			prod[r*col+c] = m[r*col+c] * other;

	return matrix(row, col, prod);
}


/** Product of matrix and scalar expression. */
matrix matrix::mul_scalar(const ex & other) const
{
	if (other.return_type() != return_types::commutative)
		throw std::runtime_error("matrix::mul_scalar(): non-commutative scalar");

	exvector prod(row * col);

	for (unsigned r=0; r<row; ++r)
		for (unsigned c=0; c<col; ++c)
			prod[r*col+c] = m[r*col+c] * other;

	return matrix(row, col, prod);
}


/** Power of a matrix.  Currently handles integer exponents only. */
matrix matrix::pow(const ex & expn) const
{
	if (col!=row)
		throw (std::logic_error("matrix::pow(): matrix not square"));
	
	if (is_exactly_a<numeric>(expn)) {
		// Integer cases are computed by successive multiplication, using the
		// obvious shortcut of storing temporaries, like A^4 == (A*A)*(A*A).
		if (expn.is_integer()) {
			numeric b = ex_to<numeric>(expn);
			matrix A(row,col);
			if (expn.info(info_flags::negative)) {
				b *= -1;
				A = this->inverse();
			} else {
				A = *this;
			}
			matrix C(row,col);
			for (unsigned r=0; r<row; ++r)
				C(r,r) = _ex1;
			if (b.is_zero())
				return C;
			// This loop computes the representation of b in base 2 from right
			// to left and multiplies the factors whenever needed.  Note
			// that this is not entirely optimal but close to optimal and
			// "better" algorithms are much harder to implement.  (See Knuth,
			// TAoCP2, section "Evaluation of Powers" for a good discussion.)
			while (b!=*_num1_p) {
				if (b.is_odd()) {
					C = C.mul(A);
					--b;
				}
				b /= *_num2_p;  // still integer.
				A = A.mul(A);
			}
			return A.mul(C);
		}
	}
	throw (std::runtime_error("matrix::pow(): don't know how to handle exponent"));
}


/** operator() to access elements for reading.
 *
 *  @param ro row of element
 *  @param co column of element
 *  @exception range_error (index out of range) */
const ex & matrix::operator() (unsigned ro, unsigned co) const
{
	if (ro>=row || co>=col)
		throw (std::range_error("matrix::operator(): index out of range"));

	return m[ro*col+co];
}


/** operator() to access elements for writing.
 *
 *  @param ro row of element
 *  @param co column of element
 *  @exception range_error (index out of range) */
ex & matrix::operator() (unsigned ro, unsigned co)
{
	if (ro>=row || co>=col)
		throw (std::range_error("matrix::operator(): index out of range"));

	ensure_if_modifiable();
	return m[ro*col+co];
}


/** Transposed of an m x n matrix, producing a new n x m matrix object that
 *  represents the transposed. */
matrix matrix::transpose() const
{
	exvector trans(this->cols()*this->rows());
	
	for (unsigned r=0; r<this->cols(); ++r)
		for (unsigned c=0; c<this->rows(); ++c)
			trans[r*this->rows()+c] = m[c*this->cols()+r];
	
	return matrix(this->cols(),this->rows(),trans);
}

/** Determinant of square matrix.  This routine doesn't actually calculate the
 *  determinant, it only implements some heuristics about which algorithm to
 *  run.  If all the elements of the matrix are elements of an integral domain
 *  the determinant is also in that integral domain and the result is expanded
 *  only.  If one or more elements are from a quotient field the determinant is
 *  usually also in that quotient field and the result is normalized before it
 *  is returned.  This implies that the determinant of the symbolic 2x2 matrix
 *  [[a/(a-b),1],[b/(a-b),1]] is returned as unity.  (In this respect, it
 *  behaves like MapleV and unlike Mathematica.)
 *
 *  @param     algo allows to chose an algorithm
 *  @return    the determinant as a new expression
 *  @exception logic_error (matrix not square)
 *  @see       determinant_algo */
ex matrix::determinant(unsigned algo) const
{
	if (row!=col)
		throw (std::logic_error("matrix::determinant(): matrix not square"));
	GINAC_ASSERT(row*col==m.capacity());
	
	// Gather some statistical information about this matrix:
	bool numeric_flag = true;
	bool normal_flag = false;
	unsigned sparse_count = 0;  // counts non-zero elements
        for (const auto& elem : m) {
		if (not elem.info(info_flags::numeric))
			numeric_flag = false;
		exmap srl;  // symbol replacement list
		ex rtest = elem.to_rational(srl);
		if (not rtest.is_zero())
			++sparse_count;
		if ((not rtest.info(info_flags::crational_polynomial)) and
		        (not rtest.info(info_flags::rational_function)))
			normal_flag = true;
	}
	
	// Here is the heuristics in case this routine has to decide:
	if (algo == determinant_algo::automatic) {
		// Minor expansion is generally a good guess:
		algo = determinant_algo::laplace;
		// Does anybody know when a matrix is really sparse?
		// Maybe <~row/2.236 nonzero elements average in a row?
		if (row>3 && 5*sparse_count<=row*col)
			algo = determinant_algo::bareiss;
		// Purely numeric matrix can be handled by Gauss elimination.
		// This overrides any prior decisions.
		if (numeric_flag)
			algo = determinant_algo::gauss;
	}
	
	// Trap the trivial case here, since some algorithms don't like it
	if (this->row==1) {
		// for consistency with non-trivial determinants...
		if (normal_flag)
			return m[0].normal();
                return m[0].expand();
	}

	// Compute the determinant
        switch (algo) {
        case determinant_algo::gauss: {
                ex det = 1;
                matrix tmp(*this);
                int sign = tmp.gauss_elimination(true);
                for (unsigned d = 0; d < row; ++d)
                        det *= tmp.m[d * col + d];
                if (normal_flag)
                        return (sign * det).normal();
                return (sign * det).normal().expand();
        }
        case determinant_algo::bareiss: {
                matrix tmp(*this);
                int sign;
                sign = tmp.fraction_free_elimination(true);
                if (normal_flag)
                        return (sign * tmp.m[row * col - 1]).normal();
                return (sign * tmp.m[row * col - 1]).expand();
        }
        case determinant_algo::divfree: {
                matrix tmp(*this);
                int sign;
                sign = tmp.division_free_elimination(true);
                if (sign == 0)
                        return _ex0;
                ex det = tmp.m[row * col - 1];
                // factor out accumulated bogus slag
                for (unsigned d = 0; d < row - 2; ++d)
                        for (unsigned j = 0; j < row - d - 2; ++j)
                                det = (det / tmp.m[d * col + d]).normal();
                return (sign * det);
        }
        case determinant_algo::laplace:
        default: {
                // This is the minor expansion scheme.  We always develop such
                // that the smallest minors (i.e, the trivial 1x1 ones) are on the
                // rightmost column.  For this to be efficient, empirical tests
                // have shown that the emptiest columns (i.e. the ones with most
                // zeros) should be the ones on the right hand side -- although
                // this might seem counter-intuitive (and in contradiction to some
                // literature like the FORM manual).  Please go ahead and test it
                // if you don't believe me!  Therefore we presort the columns of
                // the matrix:
                typedef std::pair<unsigned, unsigned> uintpair;
                std::vector<uintpair> c_zeros; // number of zeros in column
                for (unsigned c = 0; c < col; ++c) {
                        unsigned acc = 0;
                        for (unsigned rr = 0; rr < row; ++rr)
                                if (m[rr * col + c].is_zero())
                                        ++acc;
                        c_zeros.emplace_back(acc, c);
                }
                std::sort(c_zeros.begin(), c_zeros.end());
                std::vector<unsigned> pre_sort;
                for (std::vector<uintpair>::const_iterator i = c_zeros.begin(); 
                                i != c_zeros.end(); ++i)
                        pre_sort.push_back(i->second);
                std::vector<unsigned> pre_sort_test(pre_sort); // permutation_sign() modifies the vector so we make a copy here
                int sign = permutation_sign(pre_sort_test.begin(),
                                pre_sort_test.end());
                exvector result(row * col); // represents sorted matrix
                unsigned c = 0;
                for (const auto &elem : pre_sort) {
                        for (unsigned rr = 0; rr < row; ++rr)
                                result[rr * col + c] = m[rr * col + elem];
                        ++c;
                }

                if (normal_flag)
                        return (sign * matrix(row, col, result).
                                        determinant_minor()).normal();
                return sign * matrix(row, col, result).determinant_minor();
        }
        }
}


/** Trace of a matrix.  The result is normalized if it is in some quotient
 *  field and expanded only otherwise.  This implies that the trace of the
 *  symbolic 2x2 matrix [[a/(a-b),x],[y,b/(b-a)]] is recognized to be unity.
 *
 *  @return    the sum of diagonal elements
 *  @exception logic_error (matrix not square) */
ex matrix::trace() const
{
	if (row != col)
		throw (std::logic_error("matrix::trace(): matrix not square"));
	
	ex tr;
	for (unsigned r=0; r<col; ++r)
		tr += m[r*col+r];
	
	if (tr.info(info_flags::rational_function) &&
	   !tr.info(info_flags::crational_polynomial))
		return tr.normal();
	return tr.expand();
}


/** Characteristic Polynomial.  Following mathematica notation the
 *  characteristic polynomial of a matrix M is defined as the determiant of
 *  (M - lambda * 1) where 1 stands for the unit matrix of the same dimension
 *  as M.  Note that some CASs define it with a sign inside the determinant
 *  which gives rise to an overall sign if the dimension is odd.  This method
 *  returns the characteristic polynomial collected in powers of lambda as a
 *  new expression.
 *
 *  @return    characteristic polynomial as new expression
 *  @exception logic_error (matrix not square)
 *  @see       matrix::determinant() */
ex matrix::charpoly(const ex & lambda) const
{
	if (row != col)
		throw (std::logic_error("matrix::charpoly(): matrix not square"));
	
	bool numeric_flag = true;
        for (const auto & elem : m)
		if (!elem.info(info_flags::numeric)) {
			numeric_flag = false;
                        break;
                }
	
	// The pure numeric case is traditionally rather common.  Hence, it is
	// trapped and we use Leverrier's algorithm which goes as row^3 for
	// every coefficient.  The expensive part is the matrix multiplication.
	if (numeric_flag) {

		matrix B(*this);
		ex c = B.trace();
		ex poly = power(lambda, row) - c*power(lambda, row-1);
		for (unsigned i=1; i<row; ++i) {
			for (unsigned j=0; j<row; ++j)
				B.m[j*col+j] -= c;
			B = this->mul(B);
			c = B.trace() / ex(i+1);
			poly -= c*power(lambda, row-i-1);
		}
		if ((row%2) != 0u)
			return -poly;
		
			return poly;

	} else {
	
		matrix M(*this);
		for (unsigned r=0; r<col; ++r)
			M.m[r*col+r] -= lambda;
	
		return M.determinant().collect(lambda);
	}
}


/** Inverse of this matrix.
 *
 *  @return    the inverted matrix
 *  @exception logic_error (matrix not square)
 *  @exception runtime_error (singular matrix) */
matrix matrix::inverse() const
{
	if (row != col)
		throw (std::logic_error("matrix::inverse(): matrix not square"));
	
	// This routine actually doesn't do anything fancy at all.  We compute the
	// inverse of the matrix A by solving the system A * A^{-1} == Id.
	
	// First populate the identity matrix supposed to become the right hand side.
	matrix identity(row,col);
	for (unsigned i=0; i<row; ++i)
		identity(i,i) = _ex1;
	
	// Populate a dummy matrix of variables, just because of compatibility with
	// matrix::solve() which wants this (for compatibility with under-determined
	// systems of equations).
	matrix vars(row,col);
	for (unsigned r=0; r<row; ++r)
		for (unsigned c=0; c<col; ++c)
			vars(r,c) = symbol();
	
	matrix sol(row,col);
	try {
		sol = this->solve(vars,identity);
	} catch (const std::runtime_error & e) {
	    if (e.what()==std::string("matrix::solve(): inconsistent linear system"))
			throw (std::runtime_error("matrix::inverse(): singular matrix"));
		else
			throw;
	}
	return sol;
}


/** Solve a linear system consisting of a m x n matrix and a m x p right hand
 *  side by applying an elimination scheme to the augmented matrix.
 *
 *  @param vars n x p matrix, all elements must be symbols 
 *  @param rhs m x p matrix
 *  @param algo selects the solving algorithm
 *  @return n x p solution matrix
 *  @exception logic_error (incompatible matrices)
 *  @exception invalid_argument (1st argument must be matrix of symbols)
 *  @exception runtime_error (inconsistent linear system)
 *  @see       solve_algo */
matrix matrix::solve(const matrix & vars,
                     const matrix & rhs,
                     unsigned algo) const
{
	const unsigned mm = this->rows();
	const unsigned n = this->cols();
	const unsigned p = rhs.cols();
	
	// syntax checks    
	if ((rhs.rows() != mm) || (vars.rows() != n) || (vars.col != p))
		throw (std::logic_error("matrix::solve(): incompatible matrices"));
	for (unsigned ro=0; ro<n; ++ro)
		for (unsigned co=0; co<p; ++co)
			if (!vars(ro,co).info(info_flags::symbol))
				throw (std::invalid_argument("matrix::solve(): 1st argument must be matrix of symbols"));
	
	// build the augmented matrix of *this with rhs attached to the right
	matrix aug(mm,n+p);
	for (unsigned r=0; r<mm; ++r) {
		for (unsigned c=0; c<n; ++c)
			aug.m[r*(n+p)+c] = this->m[r*n+c];
		for (unsigned c=0; c<p; ++c)
			aug.m[r*(n+p)+c+n] = rhs.m[r*p+c];
	}
	
	// Gather some statistical information about the augmented matrix:
	bool numeric_flag = true;
        for (const auto & elem : aug.m)
		if (!elem.info(info_flags::numeric)) {
			numeric_flag = false;
                        break;
                }
	
	// Here is the heuristics in case this routine has to decide:
	if (algo == solve_algo::automatic) {
		// Bareiss (fraction-free) elimination is generally a good guess:
		algo = solve_algo::bareiss;
		// For m<3, Bareiss elimination is equivalent to division free
		// elimination but has more logistic overhead
		if (mm<3)
			algo = solve_algo::divfree;
		// This overrides any prior decisions.
		if (numeric_flag)
			algo = solve_algo::gauss;
	}
	
	// Eliminate the augmented matrix:
        switch (algo) {
        case solve_algo::gauss:
                aug.gauss_elimination();
                break;
        case solve_algo::divfree:
                aug.division_free_elimination();
                break;
        case solve_algo::bareiss:
        default:
                aug.fraction_free_elimination();
        }

        // assemble the solution matrix:
	matrix sol(n,p);
	for (unsigned co=0; co<p; ++co) {
		unsigned last_assigned_sol = n+1;
		for (int r=mm-1; r>=0; --r) {
			unsigned fnz = 1;    // first non-zero in row
			while ((fnz<=n) && (aug.m[r*(n+p)+(fnz-1)].normal().is_zero()))
				++fnz;
			if (fnz>n) {
				// row consists only of zeros, corresponding rhs must be 0, too
				if (!aug.m[r*(n+p)+n+co].normal().is_zero()) {
					throw (std::runtime_error("matrix::solve(): inconsistent linear system"));
				}
			} else {
				// assign solutions for vars between fnz+1 and
				// last_assigned_sol-1: free parameters
				for (unsigned c=fnz; c<last_assigned_sol-1; ++c)
					sol(c,co) = vars.m[c*p+co];
				ex e = aug.m[r*(n+p)+n+co];
				for (unsigned c=fnz; c<n; ++c)
					e -= aug.m[r*(n+p)+c]*sol.m[c*p+co];
				sol(fnz-1,co) = (e/(aug.m[r*(n+p)+(fnz-1)])).normal();
				last_assigned_sol = fnz;
			}
		}
		// assign solutions for vars between 1 and
		// last_assigned_sol-1: free parameters
		for (unsigned ro=0; ro<last_assigned_sol-1; ++ro)
			sol(ro,co) = vars(ro,co);
	}
	
	return sol;
}


/** Compute the rank of this matrix. */
unsigned matrix::rank() const
{
	// Method:
	// Transform this matrix into upper echelon form and then count the
	// number of non-zero rows.

	GINAC_ASSERT(row*col==m.capacity());

	// Actually, any elimination scheme will do since we are only
	// interested in the echelon matrix' zeros.
	matrix to_eliminate = *this;
	to_eliminate.fraction_free_elimination();

	unsigned r = row*col;  // index of last non-zero element
	while ((r--) != 0u) {
		if (!to_eliminate.m[r].is_zero())
			return 1+r/col;
	}
	return 0;
}


// protected

/** Recursive determinant for small matrices having at least one symbolic
 *  entry.  The basic algorithm, known as Laplace-expansion, is enhanced by
 *  some bookkeeping to avoid calculation of the same submatrices ("minors")
 *  more than once.  According to W.M.Gentleman and S.C.Johnson this algorithm
 *  is better than elimination schemes for matrices of sparse multivariate
 *  polynomials and also for matrices of dense univariate polynomials if the
 *  matrix' dimesion is larger than 7.
 *
 *  @return the determinant as a new expression (in expanded form)
 *  @see matrix::determinant() */
ex matrix::determinant_minor() const
{
	// for small matrices the algorithm does not make any sense:
	const unsigned n = this->cols();
	if (n==1)
		return m[0].expand();
	if (n==2)
		return (m[0]*m[3]-m[2]*m[1]).expand();
	if (n==3)
		return (m[0]*m[4]*m[8]-m[0]*m[5]*m[7]-
		        m[1]*m[3]*m[8]+m[2]*m[3]*m[7]+
		        m[1]*m[5]*m[6]-m[2]*m[4]*m[6]).expand();
	
	// This algorithm can best be understood by looking at a naive
	// implementation of Laplace-expansion, like this one:
	// ex det;
	// matrix minorM(this->rows()-1,this->cols()-1);
	// for (unsigned r1=0; r1<this->rows(); ++r1) {
	//     // shortcut if element(r1,0) vanishes
	//     if (m[r1*col].is_zero())
	//         continue;
	//     // assemble the minor matrix
	//     for (unsigned r=0; r<minorM.rows(); ++r) {
	//         for (unsigned c=0; c<minorM.cols(); ++c) {
	//             if (r<r1)
	//                 minorM(r,c) = m[r*col+c+1];
	//             else
	//                 minorM(r,c) = m[(r+1)*col+c+1];
	//         }
	//     }
	//     // recurse down and care for sign:
	//     if (r1%2)
	//         det -= m[r1*col] * minorM.determinant_minor();
	//     else
	//         det += m[r1*col] * minorM.determinant_minor();
	// }
	// return det.expand();
	// What happens is that while proceeding down many of the minors are
	// computed more than once.  In particular, there are binomial(n,k)
	// kxk minors and each one is computed factorial(n-k) times.  Therefore
	// it is reasonable to store the results of the minors.  We proceed from
	// right to left.  At each column c we only need to retrieve the minors
	// calculated in step c-1.  We therefore only have to store at most 
	// 2*binomial(n,n/2) minors.
	
	// Unique flipper counter for partitioning into minors
	std::vector<unsigned> Pkey;
	Pkey.reserve(n);
	// key for minor determinant (a subpartition of Pkey)
	std::vector<unsigned> Mkey;
	Mkey.reserve(n-1);
	// we store our subminors in maps, keys being the rows they arise from
	typedef std::map<std::vector<unsigned>,class ex> Rmap;
	typedef std::map<std::vector<unsigned>,class ex>::value_type Rmap_value;
	Rmap A;
	Rmap B;
	ex det;
	// initialize A with last column:
	for (unsigned r=0; r<n; ++r) {
		Pkey.erase(Pkey.begin(),Pkey.end());
		Pkey.push_back(r);
		A.insert(Rmap_value(Pkey,m[n*(r+1)-1]));
	}
	// proceed from right to left through matrix
	for (int c=n-2; c>=0; --c) {
		Pkey.erase(Pkey.begin(),Pkey.end());  // don't change capacity
		Mkey.erase(Mkey.begin(),Mkey.end());
		for (unsigned i=0; i<n-c; ++i)
			Pkey.push_back(i);
		unsigned fc = 0;  // controls logic for our strange flipper counter
		do {
			det = _ex0;
			for (unsigned r=0; r<n-c; ++r) {
				// maybe there is nothing to do?
				if (m[Pkey[r]*n+c].is_zero())
					continue;
				// create the sorted key for all possible minors
				Mkey.erase(Mkey.begin(),Mkey.end());
				for (unsigned i=0; i<n-c; ++i)
					if (i!=r)
						Mkey.push_back(Pkey[i]);
				// Fetch the minors and compute the new determinant
				if ((r%2) != 0u)
					det -= m[Pkey[r]*n+c]*A[Mkey];
				else
					det += m[Pkey[r]*n+c]*A[Mkey];
			}
			// prevent build-up of deep nesting of expressions saves time:
			det = det.expand();
			// store the new determinant at its place in B:
			if (!det.is_zero())
				B.insert(Rmap_value(Pkey,det));
			// increment our strange flipper counter
			for (fc=n-c; fc>0; --fc) {
				++Pkey[fc-1];
				if (Pkey[fc-1]<fc+c)
					break;
			}
			if (fc<n-c && fc>0)
				for (unsigned j=fc; j<n-c; ++j)
					Pkey[j] = Pkey[j-1]+1;
		} while(fc != 0u);
		// next column, so change the role of A and B:
		A.swap(B);
		B.clear();
	}
	
	return det;
}


/** Perform the steps of an ordinary Gaussian elimination to bring the m x n
 *  matrix into an upper echelon form.  The algorithm is ok for matrices
 *  with numeric coefficients but quite unsuited for symbolic matrices.
 *
 *  @param det may be set to true to save a lot of space if one is only
 *  interested in the diagonal elements (i.e. for calculating determinants).
 *  The others are set to zero in this case.
 *  @return sign is 1 if an even number of rows was swapped, -1 if an odd
 *  number of rows was swapped and 0 if the matrix is singular. */
int matrix::gauss_elimination(const bool det)
{
	ensure_if_modifiable();
	const unsigned mm = this->rows();
	const unsigned n = this->cols();
	GINAC_ASSERT(!det || n==mm);
	int sign = 1;
	
	unsigned r0 = 0;
	for (unsigned c0=0; c0<n && r0<mm-1; ++c0) {
		int indx = pivot(r0, c0, true);
		if (indx == -1) {
			sign = 0;
			if (det)
				return 0;  // leaves *this in a messy state
		}
		if (indx>=0) {
			if (indx > 0)
				sign = -sign;
			for (unsigned r2=r0+1; r2<mm; ++r2) {
				if (!this->m[r2*n+c0].is_zero()) {
					// yes, there is something to do in this row
					ex piv = this->m[r2*n+c0] / this->m[r0*n+c0];
					for (unsigned c=c0+1; c<n; ++c) {
						this->m[r2*n+c] -= piv * this->m[r0*n+c];
						if (!this->m[r2*n+c].info(info_flags::numeric))
							this->m[r2*n+c] = this->m[r2*n+c].normal();
					}
				}
				// fill up left hand side with zeros
				for (unsigned c=r0; c<=c0; ++c)
					this->m[r2*n+c] = _ex0;
			}
			if (det) {
				// save space by deleting no longer needed elements
				for (unsigned c=r0+1; c<n; ++c)
					this->m[r0*n+c] = _ex0;
			}
			++r0;
		}
	}
	// clear remaining rows
	for (unsigned r=r0+1; r<mm; ++r) {
		for (unsigned c=0; c<n; ++c)
			this->m[r*n+c] = _ex0;
	}

	return sign;
}


/** Perform the steps of division free elimination to bring the m x n matrix
 *  into an upper echelon form.
 *
 *  @param det may be set to true to save a lot of space if one is only
 *  interested in the diagonal elements (i.e. for calculating determinants).
 *  The others are set to zero in this case.
 *  @return sign is 1 if an even number of rows was swapped, -1 if an odd
 *  number of rows was swapped and 0 if the matrix is singular. */
int matrix::division_free_elimination(const bool det)
{
	ensure_if_modifiable();
	const unsigned mm = this->rows();
	const unsigned n = this->cols();
	GINAC_ASSERT(!det || n==mm);
	int sign = 1;
	
	unsigned r0 = 0;
	for (unsigned c0=0; c0<n && r0<mm-1; ++c0) {
		int indx = pivot(r0, c0, true);
		if (indx==-1) {
			sign = 0;
			if (det)
				return 0;  // leaves *this in a messy state
		}
		if (indx>=0) {
			if (indx>0)
				sign = -sign;
			for (unsigned r2=r0+1; r2<mm; ++r2) {
				for (unsigned c=c0+1; c<n; ++c)
					this->m[r2*n+c] = (this->m[r0*n+c0]*this->m[r2*n+c] - this->m[r2*n+c0]*this->m[r0*n+c]).expand();
				// fill up left hand side with zeros
				for (unsigned c=r0; c<=c0; ++c)
					this->m[r2*n+c] = _ex0;
			}
			if (det) {
				// save space by deleting no longer needed elements
				for (unsigned c=r0+1; c<n; ++c)
					this->m[r0*n+c] = _ex0;
			}
			++r0;
		}
	}
	// clear remaining rows
	for (unsigned r=r0+1; r<mm; ++r) {
		for (unsigned c=0; c<n; ++c)
			this->m[r*n+c] = _ex0;
	}

	return sign;
}


/** Perform the steps of Bareiss' one-step fraction free elimination to bring
 *  the matrix into an upper echelon form.  Fraction free elimination means
 *  that divide is used straightforwardly, without computing GCDs first.  This
 *  is possible, since we know the divisor at each step.
 *  
 *  @param det may be set to true to save a lot of space if one is only
 *  interested in the last element (i.e. for calculating determinants). The
 *  others are set to zero in this case.
 *  @return sign is 1 if an even number of rows was swapped, -1 if an odd
 *  number of rows was swapped and 0 if the matrix is singular. */
int matrix::fraction_free_elimination(const bool det)
{
	// Method:
	// (single-step fraction free elimination scheme, already known to Jordan)
	//
	// Usual division-free elimination sets m[0](r,c) = m(r,c) and then sets
	//     m[k+1](r,c) = m[k](k,k) * m[k](r,c) - m[k](r,k) * m[k](k,c).
	//
	// Bareiss (fraction-free) elimination in addition divides that element
	// by m[k-1](k-1,k-1) for k>1, where it can be shown by means of the
	// Sylvester identity that this really divides m[k+1](r,c).
	//
	// We also allow rational functions where the original prove still holds.
	// However, we must care for numerator and denominator separately and
	// "manually" work in the integral domains because of subtle cancellations
	// (see below).  This blows up the bookkeeping a bit and the formula has
	// to be modified to expand like this (N{x} stands for numerator of x,
	// D{x} for denominator of x):
	//     N{m[k+1](r,c)} = N{m[k](k,k)}*N{m[k](r,c)}*D{m[k](r,k)}*D{m[k](k,c)}
	//                     -N{m[k](r,k)}*N{m[k](k,c)}*D{m[k](k,k)}*D{m[k](r,c)}
	//     D{m[k+1](r,c)} = D{m[k](k,k)}*D{m[k](r,c)}*D{m[k](r,k)}*D{m[k](k,c)}
	// where for k>1 we now divide N{m[k+1](r,c)} by
	//     N{m[k-1](k-1,k-1)}
	// and D{m[k+1](r,c)} by
	//     D{m[k-1](k-1,k-1)}.
	
	ensure_if_modifiable();
	const unsigned mm = this->rows();
	const unsigned n = this->cols();
	GINAC_ASSERT(!det || n==mm);
	int sign = 1;
	if (mm==1)
		return 1;
	ex divisor_n = 1;
	ex divisor_d = 1;
	ex dividend_n;
	ex dividend_d;
	
	// We populate temporary matrices to subsequently operate on.  There is
	// one holding numerators and another holding denominators of entries.
	// This is a must since the evaluator (or even earlier mul's constructor)
	// might cancel some trivial element which causes divide() to fail.  The
	// elements are normalized first (yes, even though this algorithm doesn't
	// need GCDs) since the elements of *this might be unnormalized, which
	// makes things more complicated than they need to be.
	matrix tmp_n(*this);
	matrix tmp_d(mm,n);  // for denominators, if needed
	exmap srl;  // symbol replacement list
	exvector::const_iterator cit = this->m.begin(), citend = this->m.end();
	auto tmp_n_it = tmp_n.m.begin(), tmp_d_it = tmp_d.m.begin();
	while (cit != citend) {
		ex nd = cit->normal().to_rational(srl).numer_denom();
		++cit;
		*tmp_n_it++ = nd.op(0);
		*tmp_d_it++ = nd.op(1);
	}
	
	unsigned r0 = 0;
	for (unsigned c0=0; c0<n && r0<mm-1; ++c0) {
		// When trying to find a pivot, we should try a bit harder than expand().
		// Searching the first non-zero element in-place here instead of calling
		// pivot() allows us to do no more substitutions and back-substitutions
		// than are actually necessary.
		int indx = r0;
		while ((static_cast<unsigned>(indx)<mm) &&
		       (tmp_n[indx*n+c0].subs(srl, subs_options::no_pattern).expand().is_zero()))
			++indx;
		if (static_cast<unsigned> (indx)==mm) {
			// all elements in column c0 below row r0 vanish
			sign = 0;
			if (det)
				return 0;
		} else {
			if (static_cast<unsigned>(indx)>r0) {
				// Matrix needs pivoting, swap rows r0 and indx of tmp_n and tmp_d.
				sign = -sign;
				for (unsigned c=c0; c<n; ++c) {
					tmp_n.m[n*indx+c].swap(tmp_n.m[n*r0+c]);
					tmp_d.m[n*indx+c].swap(tmp_d.m[n*r0+c]);
				}
			}
			for (unsigned r2=r0+1; r2<mm; ++r2) {
				for (unsigned c=c0+1; c<n; ++c) {
					dividend_n = (tmp_n.m[r0*n+c0]*tmp_n.m[r2*n+c]*
					              tmp_d.m[r2*n+c0]*tmp_d.m[r0*n+c]
					             -tmp_n.m[r2*n+c0]*tmp_n.m[r0*n+c]*
					              tmp_d.m[r0*n+c0]*tmp_d.m[r2*n+c]).expand();
					dividend_d = (tmp_d.m[r2*n+c0]*tmp_d.m[r0*n+c]*
					              tmp_d.m[r0*n+c0]*tmp_d.m[r2*n+c]).expand();
					bool check = divide(dividend_n, divisor_n,
					                    tmp_n.m[r2*n+c], true);
					check &= static_cast<int>(divide(dividend_d, divisor_d,
					                tmp_d.m[r2*n+c], true));
#ifdef DO_GINAC_ASSERT
					GINAC_ASSERT(check);
#endif
				}
				// fill up left hand side with zeros
				for (unsigned c=r0; c<=c0; ++c)
					tmp_n.m[r2*n+c] = _ex0;
			}
			if (c0<n && r0<mm-1) {
				// compute next iteration's divisor
				divisor_n = tmp_n.m[r0*n+c0].expand();
				divisor_d = tmp_d.m[r0*n+c0].expand();
				if (det) {
					// save space by deleting no longer needed elements
					for (unsigned c=0; c<n; ++c) {
						tmp_n.m[r0*n+c] = _ex0;
						tmp_d.m[r0*n+c] = _ex1;
					}
				}
			}
			++r0;
		}
	}
	// clear remaining rows
	for (unsigned r=r0+1; r<mm; ++r) {
		for (unsigned c=0; c<n; ++c)
			tmp_n.m[r*n+c] = _ex0;
	}

	// repopulate *this matrix:
	auto it = this->m.begin(), itend = this->m.end();
	tmp_n_it = tmp_n.m.begin();
	tmp_d_it = tmp_d.m.begin();
	while (it != itend)
		*it++ = ((*tmp_n_it++)/(*tmp_d_it++)).subs(srl, subs_options::no_pattern);
	
	return sign;
}


/** Partial pivoting method for matrix elimination schemes.
 *  Usual pivoting (symbolic==false) returns the index to the element with the
 *  largest absolute value in column ro and swaps the current row with the one
 *  where the element was found.  With (symbolic==true) it does the same thing
 *  with the first non-zero element.
 *
 *  @param ro is the row from where to begin
 *  @param co is the column to be inspected
 *  @param symbolic signal if we want the first non-zero element to be pivoted
 *  (true) or the one with the largest absolute value (false).
 *  @return 0 if no interchange occurred, -1 if all are zero (usually signaling
 *  a degeneracy) and positive integer k means that rows ro and k were swapped.
 */
int matrix::pivot(unsigned ro, unsigned co, bool symbolic)
{
	unsigned k = ro;
	if (symbolic) {
		// search first non-zero element in column co beginning at row ro
		while ((k<row) && (this->m[k*col+co].expand().is_zero()))
			++k;
	} else {
		// search largest element in column co beginning at row ro
		GINAC_ASSERT(is_exactly_a<numeric>(this->m[k*col+co]));
		unsigned kmax = k+1;
		numeric mmax = abs(ex_to<numeric>(m[kmax*col+co]));
		while (kmax<row) {
			GINAC_ASSERT(is_exactly_a<numeric>(this->m[kmax*col+co]));
			const numeric& tmp = ex_to<numeric>(this->m[kmax*col+co]);
			if (abs(tmp) > mmax) {
				mmax = tmp;
				k = kmax;
			}
			++kmax;
		}
		if (!mmax.is_zero())
			k = kmax;
	}
	if (k==row)
		// all elements in column co below row ro vanish
		return -1;
	if (k==ro)
		// matrix needs no pivoting
		return 0;
	// matrix needs pivoting, so swap rows k and ro
	ensure_if_modifiable();
	for (unsigned c=0; c<col; ++c)
		this->m[k*col+c].swap(this->m[ro*col+c]);
	
	return k;
}

/** Function to check that all elements of the matrix are zero.
 */
bool matrix::is_zero_matrix() const
{
	for (const auto & elem : m) 
		if(!(elem.is_zero()))
			return false;
	return true;
}

ex lst_to_matrix(const lst & l)
{
	lst::const_iterator itr, itc;

	// Find number of rows and columns
	size_t rows = l.nops(), cols = 0;
	for (itr = l.begin(); itr != l.end(); ++itr) {
		if (!is_a<lst>(*itr))
			throw (std::invalid_argument("lst_to_matrix: argument must be a list of lists"));
		if (itr->nops() > cols)
			cols = itr->nops();
	}

	// Allocate and fill matrix
	matrix &M = *new matrix(rows, cols);
	M.setflag(status_flags::dynallocated);

	unsigned i;
	for (itr = l.begin(), i = 0; itr != l.end(); ++itr, ++i) {
		unsigned j;
		for (itc = ex_to<lst>(*itr).begin(), j = 0; itc != ex_to<lst>(*itr).end(); ++itc, ++j)
			M(i, j) = *itc;
	}

	return M;
}

ex diag_matrix(const lst & l)
{
	lst::const_iterator it;
	size_t dim = l.nops();

	// Allocate and fill matrix
	matrix &M = *new matrix(dim, dim);
	M.setflag(status_flags::dynallocated);

	unsigned i;
	for (it = l.begin(), i = 0; it != l.end(); ++it, ++i)
		M(i, i) = *it;

	return M;
}

ex unit_matrix(unsigned r, unsigned c)
{
	matrix &Id = *new matrix(r, c);
	Id.setflag(status_flags::dynallocated | status_flags::evaluated);
	for (unsigned i=0; i<r && i<c; i++)
		Id(i,i) = _ex1;

	return Id;
}

ex symbolic_matrix(unsigned r, unsigned c, const std::string & base_name, const std::string & tex_base_name)
{
	matrix &M = *new matrix(r, c);
	M.setflag(status_flags::dynallocated | status_flags::evaluated);

	bool long_format = (r > 10 || c > 10);
	bool single_row = (r == 1 || c == 1);

	for (unsigned i=0; i<r; i++) {
		for (unsigned j=0; j<c; j++) {
			std::ostringstream s1, s2;
			s1 << base_name;
			s2 << tex_base_name << "_{";
			if (single_row) {
				if (c == 1) {
					s1 << i;
					s2 << i << '}';
				} else {
					s1 << j;
					s2 << j << '}';
				}
			} else {
				if (long_format) {
					s1 << '_' << i << '_' << j;
					s2 << i << ';' << j << "}";
				} else {
					s1 << i << j;
					s2 << i << j << '}';
				}
			}
			M(i, j) = symbol(s1.str(), s2.str());
		}
	}

	return M;
}

ex reduced_matrix(const matrix& m, unsigned r, unsigned c)
{
	if (r+1>m.rows() || c+1>m.cols() || m.cols()<2 || m.rows()<2)
		throw std::runtime_error("minor_matrix(): index out of bounds");

	const unsigned rows = m.rows()-1;
	const unsigned cols = m.cols()-1;
	matrix &M = *new matrix(rows, cols);
	M.setflag(status_flags::dynallocated | status_flags::evaluated);

	unsigned ro = 0;
	unsigned ro2 = 0;
	while (ro2<rows) {
		if (ro==r)
			++ro;
		unsigned co = 0;
		unsigned co2 = 0;
		while (co2<cols) {
			if (co==c)
				++co;
			M(ro2,co2) = m(ro, co);
			++co;
			++co2;
		}
		++ro;
		++ro2;
	}

	return M;
}

ex sub_matrix(const matrix&m, unsigned r, unsigned nr, unsigned c, unsigned nc)
{
	if (r+nr>m.rows() || c+nc>m.cols())
		throw std::runtime_error("sub_matrix(): index out of bounds");

	matrix &M = *new matrix(nr, nc);
	M.setflag(status_flags::dynallocated | status_flags::evaluated);

	for (unsigned ro=0; ro<nr; ++ro) {
		for (unsigned co=0; co<nc; ++co) {
			M(ro,co) = m(ro+r,co+c);
		}
	}

	return M;
}

} // namespace GiNaC
