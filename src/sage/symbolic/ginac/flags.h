/** @file flags.h
 *
 *  Collection of all flags used through the GiNaC framework. */

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

#ifndef __GINAC_FLAGS_H__
#define __GINAC_FLAGS_H__

namespace GiNaC {

/** Flags to control the behavior of normal(). */
class normal_options {
public:
	enum {
		no_expand_combined_numer = 0x0001,      ///< y/x+1/(x+1) --> (y*(x+1)+y)/x/(x+1)
		no_expand_fraction_numer = 0x0002,      ///< (x-1)*(x+1)/x --> id
	};
};

/** Flags to control the behavior of expand(). */
class expand_options {
public:
	enum {
		expand_indexed = 0x0001,      ///< expands (a+b).i to a.i+b.i
		expand_function_args = 0x0002, ///< expands the arguments of functions
		expand_rename_idx = 0x0004, ///< no longer used 
		expand_transcendental = 0x0008, ///< expands trancendental functions like log and exp
		expand_only_numerators = 0x0010 ///< does not expand fraction denominators
	};
};

/** Flags to control the behavior of has(). */
class has_options {
public:
	enum {
		algebraic = 0x0001              ///< enable algebraic matching
	};
};

/** Flags to control the behavior of subs(). */
class subs_options {
public:
	enum {
		no_pattern = 0x0001,             ///< disable pattern matching
		subs_no_pattern = 0x0001, // for backwards compatibility
		algebraic = 0x0002,              ///< enable algebraic substitutions
		subs_algebraic = 0x0002,  // for backwards compatibility
		pattern_is_product = 0x0004,     ///< used internally by expairseq::subschildren()
		pattern_is_not_product = 0x0008, ///< used internally by expairseq::subschildren()
		no_index_renaming = 0x0010,
		// To indicate that we want to substitute an index by something that is
		// is not an index. Without this flag the index value would be
		// substituted in that case.
		really_subs_idx = 0x0020,
                subs_inside_piecewise = 0x0040
	};
};

/** Domain of an object */
class domain {
public:
	enum {
		complex,
		real,
		positive,
                negative,
		infinity,
                rational,
		integer,
                even
	};
};

/** Flags to control series expansion. */
class series_options {
public:
	enum {
		/** Suppress branch cuts in series expansion.  Branch cuts manifest
		 *  themselves as step functions, if this option is not passed.  If
		 *  it is passed and expansion at a point on a cut is performed, then
		 *  the analytic continuation of the function is expanded. */
		suppress_branchcut = 0x0001,
                try_univariate_flint = 0x0002,
                truncate = 0x0004
	};
};

/** Switch to control algorithm for determinant computation. */
class determinant_algo {
public:
	enum {
		/** Let the system choose.  A heuristics is applied for automatic
		 *  determination of a suitable algorithm. */
		automatic,
		/** Gauss elimination.  If \f$m_{i,j}^{(0)}\f$ are the entries of the
		 *  original matrix, then the matrix is transformed into triangular
		 *  form by applying the rules
		 *  \f[
		 *      m_{i,j}^{(k+1)} = m_{i,j}^{(k)} - m_{i,k}^{(k)} m_{k,j}^{(k)} / m_{k,k}^{(k)}
		 *  \f]
		 *  The determinant is then just the product of diagonal elements.
		 *  Choose this algorithm only for purely numerical matrices. */
		gauss,
		/** Division-free elimination.  This is a modification of Gauss
		 *  elimination where the division by the pivot element is not
		 *  carried out.  If \f$m_{i,j}^{(0)}\f$ are the entries of the
		 *  original matrix, then the matrix is transformed into triangular
		 *  form by applying the rules
		 *  \f[
		 *      m_{i,j}^{(k+1)} = m_{i,j}^{(k)} m_{k,k}^{(k)} - m_{i,k}^{(k)} m_{k,j}^{(k)}
		 *  \f]
		 *  The determinant can later be computed by inspecting the diagonal
		 *  elements only.  This algorithm is only there for the purpose of
		 *  cross-checks.  It is never fast. */
		divfree,
		/** Laplace elimination.  This is plain recursive elimination along
		 *  minors although multiple minors are avoided by the algorithm.
		 *  Although the algorithm is exponential in complexity it is
		 *  frequently the fastest one when the matrix is populated by
		 *  complicated symbolic expressions. */
		laplace,
		/** Bareiss fraction-free elimination.  This is a modification of
		 *  Gauss elimination where the division by the pivot element is
		 *  <EM>delayed</EM> until it can be carried out without computing
		 *  GCDs.  If \f$m_{i,j}^{(0)}\f$ are the entries of the original
		 *  matrix, then the matrix is transformed into triangular form by
		 *  applying the rules
		 *  \f[
		 *      m_{i,j}^{(k+1)} = (m_{i,j}^{(k)} m_{k,k}^{(k)} - m_{i,k}^{(k)} m_{k,j}^{(k)}) / m_{k-1,k-1}^{(k-1)}
		 *  \f]
		 *  (We have set \f$m_{-1,-1}^{(-1)}=1\f$ in order to avoid a case
		 *  distinction in above formula.)  It can be shown that nothing more
		 *  than polynomial long division is needed for carrying out the
		 *  division.  The determinant can then be read of from the lower
		 *  right entry.  This algorithm is rarely fast for computing
		 *  determinants. */
		bareiss
	};
};

/** Switch to control algorithm for linear system solving. */
class solve_algo {
public:
	enum {
		/** Let the system choose.  A heuristics is applied for automatic
		 *  determination of a suitable algorithm. */
		automatic,
		/** Gauss elimination.  If \f$m_{i,j}^{(0)}\f$ are the entries of the
		 *  original matrix, then the matrix is transformed into triangular
		 *  form by applying the rules
		 *  \f[
		 *      m_{i,j}^{(k+1)} = m_{i,j}^{(k)} - m_{i,k}^{(k)} m_{k,j}^{(k)} / m_{k,k}^{(k)}
		 *  \f]
		 *  This algorithm is well-suited for numerical matrices but generally
		 *  suffers from the expensive division (and computation of GCDs) at
		 *  each step. */
		gauss,
		/** Division-free elimination.  This is a modification of Gauss
		 *  elimination where the division by the pivot element is not
		 *  carried out.  If \f$m_{i,j}^{(0)}\f$ are the entries of the
		 *  original matrix, then the matrix is transformed into triangular
		 *  form by applying the rules
		 *  \f[
		 *      m_{i,j}^{(k+1)} = m_{i,j}^{(k)} m_{k,k}^{(k)} - m_{i,k}^{(k)} m_{k,j}^{(k)}
		 *  \f]
		 *  This algorithm is only there for the purpose of cross-checks.
		 *  It suffers from exponential intermediate expression swell.  Use it
		 *  only for small systems. */
		divfree,
		/** Bareiss fraction-free elimination.  This is a modification of
		 *  Gauss elimination where the division by the pivot element is
		 *  <EM>delayed</EM> until it can be carried out without computing
		 *  GCDs.  If \f$m_{i,j}^{(0)}\f$ are the entries of the original
		 *  matrix, then the matrix is transformed into triangular form by
		 *  applying the rules
		 *  \f[
		 *      m_{i,j}^{(k+1)} = (m_{i,j}^{(k)} m_{k,k}^{(k)} - m_{i,k}^{(k)} m_{k,j}^{(k)}) / m_{k-1,k-1}^{(k-1)}
		 *  \f]
		 *  (We have set \f$m_{-1,-1}^{(-1)}=1\f$ in order to avoid a case
		 *  distinction in above formula.)  It can be shown that nothing more
		 *  than polynomial long division is needed for carrying out the
		 *  division.  This is generally the fastest algorithm for solving
		 *  linear systems.  In contrast to division-free elimination it only
		 *  has a linear expression swell.  For two-dimensional systems, the
		 *  two algorithms are equivalent, however. */
		bareiss
	};
};

/** Flags to store information about the state of an object.
 *  @see basic::flags */
class status_flags {
public:
	enum {
		dynallocated    = 0x0001, ///< heap-allocated (i.e. created by new if we want to be clever and bypass the stack, @see ex::construct_from_basic() )
		evaluated       = 0x0002, ///< .eval() has already done its job
		expanded        = 0x0004, ///< .expand(0) has already done its job (other expand() options ignore this flag)
		hash_calculated = 0x0008, ///< .calchash() has already done its job
		not_shareable   = 0x0010, ///< don't share instances of this object between different expressions unless explicitly asked to (used by ex::compare())
		has_indices	= 0x0020,
		has_no_indices	= 0x0040,  // ! (has_indices || has_no_indices) means "don't know"
		is_positive	= 0x0080,
		is_negative	= 0x0100,
		purely_indefinite = 0x0200,  // If set in a mul, then it does not contains any terms with determined signs, used in power::expand()
 		tdegree_calculated	= 0x0080  // .total_degree() has already
						  // done its job (for mul)
	};
};

/** Possible attributes an object can have. */
class info_flags {
public:
	enum {
		// answered by class numeric, add, mul, some functions
		// and symbols/constants in particular domains
		numeric,
		real,
		rational,
		integer,
		crational,
		cinteger,
		positive,
		negative,
		nonnegative,
		posint,
		negint,
		nonnegint,
		even,
		odd,
		prime,
		nonzero,
		infinity,
		inexact,

		// answered by class relation
		relation,
		relation_equal,
		relation_not_equal,
		relation_less,
		relation_less_or_equal,
		relation_greater,
		relation_greater_or_equal,

		// answered by class symbol
		symbol,

		// answered by class lst
		list,

		// answered by class exprseq
		exprseq,

		// answered by classes numeric, symbol, add, mul, power
		polynomial,
		integer_polynomial,
		cinteger_polynomial,
		rational_polynomial,
		crational_polynomial,
		rational_function,
		algebraic,

		// answered by classes numeric, symbol, add, mul, power
		expanded,

		// is meaningful for mul only
		indefinite
	};
};

class return_types {
public:
	enum {
		commutative,
		noncommutative,
		noncommutative_composite
	};
};

/** Strategies how to clean up the function remember cache.
 *  @see remember_table */
class remember_strategies {
public:
	enum {
		delete_never,   ///< Let table grow undefinitely
		delete_lru,     ///< Least recently used
		delete_lfu,     ///< Least frequently used
		delete_cyclic   ///< First (oldest) one in list
	};
};

} // namespace GiNaC

#endif // ndef __GINAC_FLAGS_H__
