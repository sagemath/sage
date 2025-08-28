/** @file ex.h
 *
 *  Interface to GiNaC's light-weight expression handles. */

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

#ifndef __GINAC_EX_H__
#define __GINAC_EX_H__

#include "basic.h"
#include "ptr.h"
#include "optional.hpp"

#include <iosfwd>
#include <iterator>
#include <functional>
#include <stack>
#include <unordered_set>

class CanonicalForm;

namespace GiNaC {
#ifdef _MSC_VER
  // MSVC produces a different symbol for _ex0 when it is declared inside   
  // ex::is_zero() than when it is declared at top level as follows
  extern const ex _ex0;
#endif


/** Helper class to initialize the library.  There must be one static object
 *  of this class in every object file that makes use of our flyweights in
 *  order to guarantee proper initialization.  Hence we put it into this
 *  file which is included by every relevant file anyways.  This is modeled
 *  after section 27.4.2.1.6 of the C++ standard, where cout and friends are
 *  set up.
 *
 *  @see utils.cpp */
class library_init {
public:
	library_init();
	~library_init();
private:
	static int count;
};
/** For construction of flyweights, etc. */
static library_init library_initializer;

/** Rotate bits of unsigned value by one bit to the left.
  * This can be necessary if the user wants to define its own hashes. */
inline unsigned rotate_left(unsigned n)
{
	return (n & 0x80000000U) ? (n << 1 | 0x00000001U) : (n << 1);
}

class scalar_products;
class const_iterator;
class const_preorder_iterator;
class const_postorder_iterator;
class symbol;
struct symbolhasher;
using symbolset = std::unordered_set<symbol,symbolhasher>;
class wildcard;
struct wildhasher;
using wildset = std::unordered_set<wildcard,wildhasher>;
using expairvec = std::vector<std::pair<ex,ex>>;
using ocvector = std::vector<numeric>;
using power_ocvector_map = std::map<ex, ocvector, GiNaC::ex_is_less>;
using opt_ex = nonstd::optional<ex>;

/** Lightweight wrapper for GiNaC's symbolic objects.  It holds a pointer to
 *  the other object in order to do garbage collection by the method of
 *  reference counting.  I.e., it is a smart pointer.  Also, the constructor
 *  ex::ex(const basic & other) calls the methods that do automatic
 *  evaluation.  E.g., x-x turns automatically into 0. */
class ex {
	friend class archive_node;
	friend inline bool are_ex_trivially_equal(const ex &, const ex &);
	template<class T> friend inline const T &ex_to(const ex &);
	template<class T> friend inline bool is_a(const ex &);
	template<class T> friend inline bool is_exactly_a(const ex &);
	
	friend class print_order;
	friend class print_order_mul;
	friend class print_order_pair;
	// default constructor, copy constructor and assignment operator
public:
	ex() throw();

	// other constructors
public:
	ex(const basic & other);
	ex(int i);
	ex(unsigned int i);
	ex(long i);
	ex(double const d);
	ex(PyObject* o);

	/** Construct ex from string and a list of symbols. The input grammar is
	 *  similar to the GiNaC output format. All symbols and indices to be used
	 *  in the expression must be specified in a lst in the second argument.
	 *  Undefined symbols and other parser errors will throw an exception. */
	ex(const std::string &s, const ex &l);
	
public:
	// non-virtual functions in this class
public:
	/** Efficiently swap the contents of two expressions. */
	void swap(ex & other) throw()
	{
		GINAC_ASSERT(bp->flags & status_flags::dynallocated);
		GINAC_ASSERT(other.bp->flags & status_flags::dynallocated);
		bp.swap(other.bp);
	}

	// iterators
	const_iterator begin() const throw();
	const_iterator end() const throw();
	const_preorder_iterator preorder_begin() const;
	const_preorder_iterator preorder_end() const throw();
	const_postorder_iterator postorder_begin() const;
	const_postorder_iterator postorder_end() const throw();

	// evaluation
	ex eval(int level = 0) const { return bp->eval(level); }
	ex evalf(int level = 0, PyObject* parent=nullptr) const 
	{ return bp->evalf(level, parent); }

	// printing
	void print(const print_context & c, unsigned level = 0) const;
	void dbgprint() const;
	void dbgprinttree() const;

	// info
	bool info(unsigned inf) const { return bp->info(inf); }
	bool is_integer() const { return bp->is_integer(); }
	bool is_real() const { return bp->is_real(); }
	bool is_positive() const { return bp->is_positive(); }
	bool is_negative() const { return bp->is_negative(); }

	// operand access
	size_t nops() const { return bp->nops(); }
        size_t treesize() const;
	size_t nsymbols() const;
        bool get_first_symbol(ex &x) const;
        symbolset symbols() const;
        symbolset free_symbols() const;
        std::unordered_set<unsigned> functions() const;
        wildset wilds() const;
	const ex op(size_t i) const { return bp->op(i); }
	ex sorted_op(size_t i) const;
	ex operator[](const ex & index) const { return (*bp)[index]; }
	ex operator[](size_t i) const { return (*bp)[i]; }
	ex & let_op(size_t i);
	ex & operator[](const ex & index);
	ex & operator[](size_t i);
        void set_epseq_from(size_t i, ex e) { bp->set_epseq_from(i, e); }
	ex lhs() const;
	ex rhs() const;

	// function for complex expressions
	ex conjugate() const { return bp->conjugate(); }
	ex real_part() const { return bp->real_part(); }
	ex imag_part() const { return bp->imag_part(); }

	// pattern matching
	bool has(const ex & pattern, unsigned options = 0) const
            { return bp->has(pattern, options); }
	bool find(const ex & pattern, lst & found) const;
	bool match(const ex & pattern) const;
	bool match(const ex & pattern, lst & repl_lst) const;
	bool match(const ex & pattern, exmap& map) const
            { return bp->match(pattern, map); }
	bool match(const ex & pattern, exvector& vec) const;
        std::vector<exmap> all_matches(const ex & pattern) const;

	// substitutions
	ex subs(const exmap & m, unsigned options = 0) const;
	ex subs(const lst & ls, const lst & lr, unsigned options = 0) const;
	ex subs(const ex & e, unsigned options = 0) const;

	// function mapping
	ex map(map_function & f) const { return bp->map(f); }
	ex map(ex (*f)(const ex & e)) const;

	// visitors and tree traversal
	void accept(visitor & v) const { bp->accept(v); }
	void traverse_preorder(visitor & v) const;
	void traverse_postorder(visitor & v) const;
	void traverse(visitor & v) const { traverse_preorder(v); }

	// degree/coeff
	bool is_polynomial(const ex & vars) const;
	numeric degree(const ex & s) const;
	numeric ldegree(const ex & s) const;
	ex coeff(const ex & s, const ex & n) const { return bp->coeff(s, n); }
	ex lcoeff(const ex & s) const;
	ex tcoeff(const ex & s) const;
        void coefficients(const ex & s, expairvec & vec) const;

	// expand/collect
	ex expand(unsigned options=0) const;
	ex collect(const ex & s, bool distributed = false) const { return bp->collect(s, distributed); }
        ex collect_powers() const;

	// differentiation and series expansion
	ex diff(const symbol & s, unsigned nth = 1) const;
	ex series(const ex & r, int order, unsigned options = 0) const;
        void useries(flint_series_t& fp, int order) const { return bp->useries(fp, order); }

	// rational functions
	ex normal(int level = 0, bool noexpand_combined=false,
                        bool noexpand_numer=true) const;
	ex to_rational(exmap & repl) const;
	ex to_rational(lst & repl_lst) const;
	ex to_polynomial(exmap & repl) const;
	ex to_polynomial(lst & repl_lst) const;
        const CanonicalForm to_canonical(ex_int_umap& map,
                        power_ocvector_map& pomap, exvector& revmap) const;
        void collect_powers(power_ocvector_map& pomap) const;
	ex numer() const;
	ex denom() const;
	ex numer_denom() const;
        ex combine_fractions(bool deep=true) const;

	// polynomial algorithms
	ex unit(const ex &x) const;
	ex content(const ex &x) const;
	numeric integer_content() const;
	ex primpart(const ex &x) const;
	ex primpart(const ex &x, const ex &cont) const;
	void unitcontprim(const ex &x, ex &u, ex &c, ex &p) const;
	ex smod(const numeric &xi) const { return bp->smod(xi); }
	numeric max_coefficient() const;
        bool is_linear(const symbol& x, ex& a, ex& b) const;
        bool is_quadratic(const symbol& x, ex& a, ex& b, ex& c) const;
        bool is_binomial(const symbol& x, ex& a, ex& j, ex& b, ex& n) const;

	// domains
	void set_domain(unsigned d);

	// comparison
	int compare(const ex & other) const;
	bool is_equal(const ex & other) const;
	bool is_zero() const;
        bool is_one() const;
        bool is_minus_one() const;
        bool is_negative_or_minus() const;
        bool is_num_integer() const;
        bool is_num_fraction() const;
	
	// noncommutativity
	unsigned return_type() const { return bp->return_type(); }
	tinfo_t return_type_tinfo() const { return bp->return_type_tinfo(); }

	long gethash() const { return bp->gethash(); }

	static ptr<basic> construct_from_basic(const basic & other);
	static basic & construct_from_int(int i);
	static basic & construct_from_pyobject(PyObject* o);
	static basic & construct_from_uint(unsigned int i);
	static basic & construct_from_long(long i);
	static basic & construct_from_double(double d);
	static ptr<basic> construct_from_string_and_lst(const std::string &s, const ex &l);
	void makewriteable();
	void share(const ex & other) const;
        static ex deep_combine_fractions(ex e);
        struct combine_map_function : public map_function {
                ex operator()(const ex & e) override { return deep_combine_fractions(e); }
        };


// member variables

	mutable ptr<basic> bp;  ///< pointer to basic object managed by this
};

struct ex_hash {
    long operator()(const ex& e) const { return e.gethash(); }
};


// performance-critical inlined method implementations

// This needs to be a basic* because we don't know that numeric is derived
// from basic and we need a basic& for the ex default constructor
extern const basic *_num0_bp;

inline
ex::ex() throw() : bp(*const_cast<basic *>(_num0_bp))
{
	GINAC_ASSERT(bp->flags & status_flags::dynallocated);
}

inline
ex::ex(const basic & other) : bp(construct_from_basic(other))
{
	GINAC_ASSERT(bp->flags & status_flags::dynallocated);
}

inline
ex::ex(int i) : bp(construct_from_int(i))
{
	GINAC_ASSERT(bp->flags & status_flags::dynallocated);
}

inline
ex::ex(PyObject* o) : bp(construct_from_pyobject(o))
{
	GINAC_ASSERT(bp->flags & status_flags::dynallocated);
}

inline
ex::ex(unsigned int i) : bp(construct_from_uint(i))
{
	GINAC_ASSERT(bp->flags & status_flags::dynallocated);
}

inline
ex::ex(long i) : bp(construct_from_long(i))
{
	GINAC_ASSERT(bp->flags & status_flags::dynallocated);
}

inline
ex::ex(double const d) : bp(construct_from_double(d))
{
	GINAC_ASSERT(bp->flags & status_flags::dynallocated);
}

inline
ex::ex(const std::string &s, const ex &l) : bp(construct_from_string_and_lst(s, l))
{
	GINAC_ASSERT(bp->flags & status_flags::dynallocated);
}

// Iterators

class const_iterator : public std::iterator<std::random_access_iterator_tag, ex, ptrdiff_t, const ex *, const ex &> {
	friend class ex;
	friend class const_preorder_iterator;
	friend class const_postorder_iterator;

public:
	const_iterator() throw() {}

private:
	const_iterator(ex e_, size_t i_) throw() : e(std::move(e_)), i(i_) {}

public:
	// This should return an ex&, but that would be a reference to a
	// temporary value
	ex operator*() const
	{
		return e.op(i);
	}

	// This should return an ex*, but that would be a pointer to a
	// temporary value
	std::unique_ptr<ex> operator->() const
	{
		return std::unique_ptr<ex>(new ex(operator*()));
	}

	ex operator[](difference_type n) const
	{
		return e.op(i + n);
	}

	const_iterator &operator++() throw()
	{
		++i;
		return *this;
	}

	const_iterator operator++(int) throw()
	{
		const_iterator tmp = *this;
		++i;
		return tmp;
	}

	const_iterator &operator+=(difference_type n) throw()
	{
		i += n;
		return *this;
	}

	const_iterator operator+(difference_type n) const throw()
	{
		return const_iterator(e, i + n);
	}

	inline friend const_iterator operator+(difference_type n, const const_iterator &it) throw()
	{
		return const_iterator(it.e, it.i + n);
	}

	const_iterator &operator--() throw()
	{
		--i;
		return *this;
	}

	const_iterator operator--(int) throw()
	{
		const_iterator tmp = *this;
		--i;
		return tmp;
	}

	const_iterator &operator-=(difference_type n) throw()
	{
		i -= n;
		return *this;
	}

	const_iterator operator-(difference_type n) const throw()
	{
		return const_iterator(e, i - n);
	}

	inline friend difference_type operator-(const const_iterator &lhs, const const_iterator &rhs) throw()
	{
		return lhs.i - rhs.i;
	}

	bool operator==(const const_iterator &other) const throw()
	{
		return are_ex_trivially_equal(e, other.e) && i == other.i;
	}

	bool operator!=(const const_iterator &other) const throw()
	{
		return !(*this == other);
	}

	bool operator<(const const_iterator &other) const throw()
	{
		return i < other.i;
	}

	bool operator>(const const_iterator &other) const throw()
	{
		return other < *this;
	}

	bool operator<=(const const_iterator &other) const throw()
	{
		return !(other < *this);
	}

	bool operator>=(const const_iterator &other) const throw()
	{
		return !(*this < other);
	}

protected:
	ex e; // this used to be a "const basic *", but in view of object fusion that wouldn't be safe
	size_t i;
};

namespace internal {

struct _iter_rep {
	_iter_rep(ex e_, size_t i_, size_t i_end_) : e(std::move(e_)), i(i_), i_end(i_end_) {}

	bool operator==(const _iter_rep &other) const throw()
	{
		return are_ex_trivially_equal(e, other.e) && i == other.i;
	}

	bool operator!=(const _iter_rep &other) const throw()
	{
		return !(*this == other);
	}

	ex e;
	size_t i;
	size_t i_end;
};

} // namespace internal

class const_preorder_iterator : public std::iterator<std::forward_iterator_tag, ex, ptrdiff_t, const ex *, const ex &> {
public:
	const_preorder_iterator() throw() {}

	const_preorder_iterator(const ex &e, size_t n)
	{
		s.push(internal::_iter_rep(e, 0, n));
	}

public:
	reference operator*() const
	{
		return s.top().e;
	}

	pointer operator->() const
	{
		return &(s.top().e);
	}

	const_preorder_iterator &operator++()
	{
		increment();
		return *this;
	}

	const_preorder_iterator operator++(int)
	{
		const_preorder_iterator tmp = *this;
		increment();
		return tmp;
	}

	bool operator==(const const_preorder_iterator &other) const throw()
	{
		return s == other.s;
	}

	bool operator!=(const const_preorder_iterator &other) const throw()
	{
		return !(*this == other);
	}

private:
	std::stack<internal::_iter_rep, std::vector<internal::_iter_rep> > s;

	void increment()
	{
		while (!s.empty() && s.top().i == s.top().i_end) {
			s.pop();
			if (s.empty())
				return;
			++s.top().i;
		}

		internal::_iter_rep & current = s.top();

		if (current.i != current.i_end) {
			const ex & child = current.e.op(current.i);
			s.push(internal::_iter_rep(child, 0, child.nops()));
		}
	}
};

class const_postorder_iterator : public std::iterator<std::forward_iterator_tag, ex, ptrdiff_t, const ex *, const ex &> {
public:
	const_postorder_iterator() throw() {}

	const_postorder_iterator(const ex &e, size_t n)
	{
		s.push(internal::_iter_rep(e, 0, n));
		descend();
	}

public:
	reference operator*() const
	{
		return s.top().e;
	}

	pointer operator->() const
	{
		return &(s.top().e);
	}

	const_postorder_iterator &operator++()
	{
		increment();
		return *this;
	}

	const_postorder_iterator operator++(int)
	{
		const_postorder_iterator tmp = *this;
		increment();
		return tmp;
	}

	bool operator==(const const_postorder_iterator &other) const throw()
	{
		return s == other.s;
	}

	bool operator!=(const const_postorder_iterator &other) const throw()
	{
		return !(*this == other);
	}

private:
	std::stack<internal::_iter_rep, std::vector<internal::_iter_rep> > s;

	void descend()
	{
		while (s.top().i != s.top().i_end) {
			internal::_iter_rep & current = s.top();
			const ex & child = current.e.op(current.i);
			s.push(internal::_iter_rep(child, 0, child.nops()));
		}
	}

	void increment()
	{
		if (s.top().i == s.top().i_end)
			s.pop();
		if (!s.empty()) {
			++s.top().i;
			descend();
		}
	}
};

inline const_iterator ex::begin() const throw()
{
	return const_iterator(*this, 0);
}

inline const_iterator ex::end() const throw()
{
	return const_iterator(*this, nops());
}

inline const_preorder_iterator ex::preorder_begin() const
{
	return const_preorder_iterator(*this, nops());
}

inline const_preorder_iterator ex::preorder_end() const throw()
{
	return const_preorder_iterator();
}

inline const_postorder_iterator ex::postorder_begin() const
{
	return const_postorder_iterator(*this, nops());
}

inline const_postorder_iterator ex::postorder_end() const throw()
{
	return const_postorder_iterator();
}

/** Compare two objects of class quickly without doing a deep tree traversal.
 *  @return "true" if they are equal
 *          "false" if equality cannot be established quickly (e1 and e2 may
 *          still be equal, in this case. */
inline bool are_ex_trivially_equal(const ex &e1, const ex &e2)
{
	return e1.bp == e2.bp;
}


// Make it possible to print exvectors and exmaps
std::ostream & operator<<(std::ostream & os, const exvector & e);
std::ostream & operator<<(std::ostream & os, const exset & e);
std::ostream & operator<<(std::ostream & os, const exmap & e);

/* Function objects for STL sort() etc. */
struct ex_is_less {
	bool operator() (const ex &lh, const ex &rh) const { return lh.compare(rh) < 0; }
};

struct ex_is_equal {
	bool operator() (const ex &lh, const ex &rh) const { return lh.is_equal(rh); }
};

struct op0_is_equal {
	bool operator() (const ex &lh, const ex &rh) const { return lh.op(0).is_equal(rh.op(0)); }
};

struct ex_swap {
	void operator() (ex &lh, ex &rh) const { lh.swap(rh); }
};

/* Convert function pointer to function object suitable for map(). */
class pointer_to_map_function : public map_function {
protected:
	ex (*ptr)(const ex &);
public:
	explicit pointer_to_map_function(ex x(const ex &)) : ptr(x) {}
	ex operator()(const ex & e) override { return ptr(e); }
};

template<class T1>
class pointer_to_map_function_1arg : public map_function {
protected:
	ex (*ptr)(const ex &, T1);
	T1 arg1;
public:
	explicit pointer_to_map_function_1arg(ex x(const ex &, T1), T1 a1) : ptr(x), arg1(a1) {}
	ex operator()(const ex & e) override { return ptr(e, arg1); }
};

template<class T1, class T2>
class pointer_to_map_function_2args : public map_function {
protected:
	ex (*ptr)(const ex &, T1, T2);
	T1 arg1;
	T2 arg2;
public:
	explicit pointer_to_map_function_2args(ex x(const ex &, T1, T2), T1 a1, T2 a2) : ptr(x), arg1(a1), arg2(a2) {}
	ex operator()(const ex & e) override { return ptr(e, arg1, arg2); }
};

template<class T1, class T2, class T3>
class pointer_to_map_function_3args : public map_function {
protected:
	ex (*ptr)(const ex &, T1, T2, T3);
	T1 arg1;
	T2 arg2;
	T3 arg3;
public:
	explicit pointer_to_map_function_3args(ex x(const ex &, T1, T2, T3), T1 a1, T2 a2, T3 a3) : ptr(x), arg1(a1), arg2(a2), arg3(a3) {}
	ex operator()(const ex & e) override { return ptr(e, arg1, arg2, arg3); }
};

template<class C>
class pointer_to_member_to_map_function : public map_function {
protected:
	ex (C::*ptr)(const ex &);
	C &c;
public:
	explicit pointer_to_member_to_map_function(ex (C::*member)(const ex &), C &obj) : ptr(member), c(obj) {}
	ex operator()(const ex & e) override { return (c.*ptr)(e); }
};

template<class C, class T1>
class pointer_to_member_to_map_function_1arg : public map_function {
protected:
	ex (C::*ptr)(const ex &, T1);
	C &c;
	T1 arg1;
public:
	explicit pointer_to_member_to_map_function_1arg(ex (C::*member)(const ex &, T1), C &obj, T1 a1) : ptr(member), c(obj), arg1(a1) {}
	ex operator()(const ex & e) override { return (c.*ptr)(e, arg1); }
};

template<class C, class T1, class T2>
class pointer_to_member_to_map_function_2args : public map_function {
protected:
	ex (C::*ptr)(const ex &, T1, T2);
	C &c;
	T1 arg1;
	T2 arg2;
public:
	explicit pointer_to_member_to_map_function_2args(ex (C::*member)(const ex&, T1, T2), C &obj, T1 a1, T2 a2) : ptr(member), c(obj), arg1(a1), arg2(a2) {}
	ex operator()(const ex & e) override { return (c.*ptr)(e, arg1, arg2); }
};

template<class C, class T1, class T2, class T3>
class pointer_to_member_to_map_function_3args : public map_function {
protected:
	ex (C::*ptr)(const ex &, T1, T2, T3);
	C &c;
	T1 arg1;
	T2 arg2;
	T3 arg3;
public:
	explicit pointer_to_member_to_map_function_3args(ex (C::*member)(const ex &, T1, T2, T3), C &obj, T1 a1, T2 a2, T3 a3) : ptr(member), c(obj), arg1(a1), arg2(a2), arg3(a3) {}
	ex operator()(const ex & e) override { return (c.*ptr)(e, arg1, arg2, arg3); }
};

inline ex ex::map(ex f(const ex &)) const
{
	pointer_to_map_function fcn(f);
	return bp->map(fcn);
}

// convenience type checker template functions

/** Check if ex is a handle to a T, including base classes. */
template <class T>
inline bool is_a(const ex &obj)
{
	return is_a<T>(*obj.bp);
}

/** Check if ex is a handle to a T, not including base classes. */
template <class T>
inline bool is_exactly_a(const ex &obj)
{
	return is_exactly_a<T>(*obj.bp);
}

/** Return a reference to the basic-derived class T object embedded in an
 *  expression.  This is fast but unsafe: the result is undefined if the
 *  expression does not contain a T object at its top level.  Hence, you
 *  should generally check the type of e first.  Also, you shouldn't cache
 *  the returned reference because GiNaC's garbage collector may destroy
 *  the referenced object any time it's used in another expression.
 *
 *  @param e expression
 *  @return reference to object of class T
 *  @see is_exactly_a<class T>() */
template <class T>
inline const T &ex_to(const ex &e)
{
	GINAC_ASSERT(is_a<T>(e));
	return static_cast<const T &>(*e.bp);
}


} // namespace GiNaC

#endif // ndef __GINAC_EX_H__
