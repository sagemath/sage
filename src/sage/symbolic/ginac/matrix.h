/** @file matrix.h
 *
 *  Interface to symbolic matrices */

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

#ifndef __GINAC_MATRIX_H__
#define __GINAC_MATRIX_H__

#include "basic.h"
#include "ex.h"

#include <vector>
#include <string>

namespace GiNaC {


/** Helper template to allow initialization of matrices via an overloaded
 *  comma operator (idea stolen from Blitz++). */
template <typename T, typename It>
class matrix_init {
public:
	matrix_init(It i) : iter(std::move(i)) {}

	matrix_init<T, It> operator,(const T & x)
	{
		*iter = x;
		return matrix_init<T, It>(++iter);
	}

	// The following specializations produce much tighter code than the
	// general case above

	matrix_init<T, It> operator,(int x)
	{
		*iter = T(x);
		return matrix_init<T, It>(++iter);
	}

	matrix_init<T, It> operator,(unsigned int x)
	{
		*iter = T(x);
		return matrix_init<T, It>(++iter);
	}

	matrix_init<T, It> operator,(long x)
	{
		*iter = T(x);
		return matrix_init<T, It>(++iter);
	}

	matrix_init<T, It> operator,(unsigned long x)
	{
		*iter = T(x);
		return matrix_init<T, It>(++iter);
	}

	matrix_init<T, It> operator,(double x)
	{
		*iter = T(x);
		return matrix_init<T, It>(++iter);
	}

	matrix_init<T, It> operator,(const symbol & x)
	{
		*iter = T(x);
		return matrix_init<T, It>(++iter);
	}

private:
	matrix_init();
	It iter;
};


/** Symbolic matrices. */
class matrix : public basic
{
	GINAC_DECLARE_REGISTERED_CLASS(matrix, basic)
	
	// other constructors
public:
	matrix(unsigned r, unsigned c);
	matrix(unsigned r, unsigned c, exvector  m2);
	matrix(unsigned r, unsigned c, const lst & l);

	// First step of initialization of matrix with a comma-separated seqeuence
	// of expressions. Subsequent steps are handled by matrix_init<>::operator,().
	matrix_init<ex, exvector::iterator> operator=(const ex & x)
	{
		m[0] = x;
		return matrix_init<ex, exvector::iterator>(++m.begin());
	}
        static ex unarchive(const archive_node &n, lst &sym_lst)
        {
                return (new matrix(n, sym_lst))->
                        setflag(status_flags::dynallocated);
        }
	
	// functions overriding virtual functions from base classes
public:
	size_t nops() const override;
	const ex op(size_t i) const override;
	ex & let_op(size_t i) override;
	ex eval(int level=0) const override;
	ex evalm() const override {return *this;}
	ex subs(const exmap & m, unsigned options = 0) const override;
	ex conjugate() const override;
	ex real_part() const override;
	ex imag_part() const override;

protected:
	bool match_same_type(const basic & other) const override;
	unsigned return_type() const override { return return_types::noncommutative; };
	
	// non-virtual functions in this class
public:
	unsigned rows() const        /// Get number of rows.
		{ return row; }
	unsigned cols() const        /// Get number of columns.
		{ return col; }
	matrix add(const matrix & other) const;
	matrix sub(const matrix & other) const;
	matrix mul(const matrix & other) const;
	matrix mul(const numeric & other) const;
	matrix mul_scalar(const ex & other) const;
	matrix pow(const ex & expn) const;
	const ex & operator() (unsigned ro, unsigned co) const;
	ex & operator() (unsigned ro, unsigned co);
	matrix & set(unsigned ro, unsigned co, const ex & value) { (*this)(ro, co) = value; return *this; }
	matrix transpose() const;
	ex determinant(unsigned algo = determinant_algo::automatic) const;
	ex trace() const;
	ex charpoly(const ex & lambda) const;
	matrix inverse() const;
	matrix solve(const matrix & vars, const matrix & rhs,
	             unsigned algo = solve_algo::automatic) const;
	unsigned rank() const;
	bool is_zero_matrix() const;
protected:
	ex determinant_minor() const;
	int gauss_elimination(const bool det = false);
	int division_free_elimination(const bool det = false);
	int fraction_free_elimination(const bool det = false);
	int pivot(unsigned ro, unsigned co, bool symbolic = true);

	void print_elements(const print_context & c, const char *row_start, const char *row_end, const char *row_sep, const char *col_sep) const;
	void do_print(const print_context & c, unsigned level) const override;
	void do_print_latex(const print_latex & c, unsigned level) const;
	void do_print_python_repr(const print_python_repr & c, unsigned level) const override;
	
// member variables
protected:
	unsigned row;             ///< number of rows
	unsigned col;             ///< number of columns
	exvector m;               ///< representation (cols indexed first)
};


// wrapper functions around member functions

inline size_t nops(const matrix & m)
{ return m.nops(); }

inline ex expand(const matrix & m, unsigned options = 0)
{ return m.expand(options); }

inline ex eval(const matrix & m, int level = 0)
{ return m.eval(level); }

inline ex evalf(const matrix & m, int level = 0)
{ return m.evalf(level); }

inline unsigned rows(const matrix & m)
{ return m.rows(); }

inline unsigned cols(const matrix & m)
{ return m.cols(); }

inline matrix transpose(const matrix & m)
{ return m.transpose(); }

inline ex determinant(const matrix & m, unsigned options = determinant_algo::automatic)
{ return m.determinant(options); }

inline ex trace(const matrix & m)
{ return m.trace(); }

inline ex charpoly(const matrix & m, const ex & lambda)
{ return m.charpoly(lambda); }

inline matrix inverse(const matrix & m)
{ return m.inverse(); }

inline unsigned rank(const matrix & m)
{ return m.rank(); }

// utility functions

/** Convert list of lists to matrix. */
extern ex lst_to_matrix(const lst & l);

/** Convert list of diagonal elements to matrix. */
extern ex diag_matrix(const lst & l);

/** Create an r times c unit matrix. */
extern ex unit_matrix(unsigned r, unsigned c);

/** Create a x times x unit matrix. */
inline ex unit_matrix(unsigned x)
{ return unit_matrix(x, x); }

/** Create an r times c matrix of newly generated symbols consisting of the
 *  given base name plus the numeric row/column position of each element.
 *  The base name for LaTeX output is specified separately. */
extern ex symbolic_matrix(unsigned r, unsigned c, const std::string & base_name, const std::string & tex_base_name);

/** Return the reduced matrix that is formed by deleting the rth row and cth
 *  column of matrix m. The determinant of the result is the Minor r, c. */
extern ex reduced_matrix(const matrix& m, unsigned r, unsigned c);

/** Return the nr times nc submatrix starting at position r, c of matrix m. */
extern ex sub_matrix(const matrix&m, unsigned r, unsigned nr, unsigned c, unsigned nc);

/** Create an r times c matrix of newly generated symbols consisting of the
 *  given base name plus the numeric row/column position of each element. */
inline ex symbolic_matrix(unsigned r, unsigned c, const std::string & base_name)
{ return symbolic_matrix(r, c, base_name, base_name); }

} // namespace GiNaC

#endif // ndef __GINAC_MATRIX_H__
