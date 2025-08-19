/** @file basic.h
 *
 *  Interface to GiNaC's ABC. */

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

#ifndef __GINAC_BASIC_H__
#define __GINAC_BASIC_H__

#include <vector>
#include <set>
#include <map>
#include <unordered_map>
// CINT needs <algorithm> to work properly with <vector>
#include <algorithm>

#include "pynac-config.h"
#include "flags.h"
#include "ptr.h"
#include "assertion.h"
#include "registrar.h"
#include "context.h"

// PyObject forward declaration 
#ifndef PyObject_HEAD
struct _object;
typedef _object PyObject;
#endif

namespace GiNaC {

class ex;
struct ex_is_less;
struct ex_hash;
class symbol;
class numeric;
class relational;
class archive_node;
class print_context;
class flint_series_t;

typedef std::vector<ex> exvector;
typedef std::set<ex, ex_is_less> exset;
typedef std::map<ex, ex, ex_is_less> exmap;
using ex_int_map = std::map<GiNaC::ex, int, GiNaC::ex_is_less>;
using ex_int_umap = std::unordered_map<ex, int, ex_hash>;

// Define this to enable some statistical output for comparisons and hashing
#undef GINAC_COMPARE_STATISTICS

#ifdef GINAC_COMPARE_STATISTICS
class compare_statistics_t {
public:
	compare_statistics_t()
	 : total_compares(0), nontrivial_compares(0), total_basic_compares(0), compare_same_hashvalue(0), compare_same_type(0),
	   total_is_equals(0), nontrivial_is_equals(0), total_basic_is_equals(0), is_equal_same_hashvalue(0), is_equal_same_type(0),
	   total_gethash(0), gethash_cached(0) {}
	~compare_statistics_t();

	unsigned long total_compares;
	unsigned long nontrivial_compares;
	unsigned long total_basic_compares;
	unsigned long compare_same_hashvalue;
	unsigned long compare_same_type;

	unsigned long total_is_equals;
	unsigned long nontrivial_is_equals;
	unsigned long total_basic_is_equals;
	unsigned long is_equal_same_hashvalue;
	unsigned long is_equal_same_type;

	unsigned long total_gethash;
	unsigned long gethash_cached;
};

extern compare_statistics_t compare_statistics;
#endif


/** Function object for map(). */
struct map_function {
	virtual ~map_function() {}
	typedef const ex & argument_type;
	typedef ex result_type;
	virtual ex operator()(const ex & e) = 0;
};


/** Degenerate base class for visitors. basic and derivative classes
 *  support Robert C. Martin's Acyclic Visitor pattern (cf.
 *  http://objectmentor.com/publications/acv.pdf). */
class visitor {
protected:
	virtual ~visitor() {}
};


/** This class is the ABC (abstract base class) of GiNaC's class hierarchy. */
class basic : public refcounted
{
public:
	typedef void inherited;
        static const tinfo_static_t tinfo_static;
private:
	static registered_class_info reg_info;
public:
	static registered_class_info &get_class_info_static() { return reg_info; }
	virtual const registered_class_info &get_class_info() const { return basic::get_class_info_static(); }
	virtual GiNaC::registered_class_info &get_class_info() { return basic::get_class_info_static(); }
	virtual const char *class_name() const { return basic::get_class_info_static().options.get_name(); }

	basic(const archive_node &n, lst &sym_lst);
	virtual void archive(archive_node &n) const;

	class visitor {
	public:
		virtual void visit(const basic &) = 0;
		virtual ~visitor() {}; \
	};

	friend class ex;
	friend class print_order;
	friend class print_order_pair;
	// default constructor, destructor, copy constructor and assignment operator
protected:
	basic() : tinfo_key(&tinfo_static), flags(0) {}

public:
	/** basic destructor, virtual because class ex will delete objects of
	 *  derived classes via a basic*. */
	virtual ~basic()
	{
		GINAC_ASSERT((!(flags & status_flags::dynallocated)) || (get_refcount() == 0));
	}
	basic(const basic & other);
	const basic & operator=(const basic & other);

protected:
	/** Constructor with specified tinfo_key (used by derived classes instead
	 *  of the default constructor to avoid assigning tinfo_key twice). */
	basic(tinfo_t ti) : tinfo_key(ti), flags(0) {}

	// new virtual functions which can be overridden by derived classes
public: // only const functions please (may break reference counting)

	/** Create a clone of this object on the heap.  One can think of this as
	 *  simulating a virtual copy constructor which is needed for instance by
	 *  the refcounted construction of an ex from a basic. */
	virtual basic * duplicate() const { return new basic(*this); }

	// evaluation
	virtual ex eval(int level = 0) const;
	virtual ex evalf(int level = 0, PyObject* parent=nullptr) const;
	virtual ex evalm() const;

	// printing
	virtual void print(const print_context & c, unsigned level = 0) const;
	virtual void dbgprint() const;
	virtual void dbgprinttree() const;
	virtual unsigned precedence() const;

	// info
	virtual bool info(unsigned inf) const;
        virtual bool is_integer() const { return info(info_flags::integer); }
        virtual bool is_real() const { return info(info_flags::real); }
        virtual bool is_positive() const { return info(info_flags::positive); }
        virtual bool is_negative() const { return info(info_flags::negative); }

	// operand access
	virtual size_t nops() const;
	virtual const ex op(size_t i) const;
	virtual ex operator[](const ex & index) const;
	virtual ex operator[](size_t i) const;
	virtual ex & let_op(size_t i);
	virtual ex & operator[](const ex & index);
	virtual ex & operator[](size_t i);

	// pattern matching
	virtual bool has(const ex & other, unsigned options = 0) const;
	virtual bool match(const ex & pattern, exmap& map) const;
	virtual bool match_same_type(const basic & other) const;

	virtual void do_print(const print_context & c, unsigned level) const;
	virtual void do_print_tree(const print_tree & c, unsigned level) const;
	virtual void do_print_python_repr(const print_python_repr & c,
                        unsigned level) const;

	// substitutions
	virtual ex subs(const exmap & m, unsigned options = 0) const;

	// function mapping
	virtual ex map(map_function & f) const;

	// visitors and tree traversal
	virtual void accept(GiNaC::visitor & v) const
	{
		if (visitor *p = dynamic_cast<visitor *>(&v))
			p->visit(*this);
	}

	// degree/coeff
	virtual bool is_polynomial(const ex & var) const;
	virtual numeric degree(const ex & s) const;
	virtual numeric ldegree(const ex & s) const;
	virtual ex coeff(const ex & s, const ex & n) const;

	// expand/collect
	virtual ex expand(unsigned options = 0) const;
	virtual ex collect(const ex & s, bool distributed = false) const;

	// differentiation and series expansion
	virtual ex derivative(const symbol & s) const;
	virtual ex series(const relational & r, int order,
                        unsigned options = 0) const;
        virtual void useries(flint_series_t& fp, int order) const {}

	// rational functions
	virtual ex normal(exmap & repl, exmap & rev_lookup, int level = 0,
                        unsigned options = 0) const;
	virtual ex to_rational(exmap & repl) const;
	virtual ex to_polynomial(exmap & repl) const;

	// polynomial algorithms
	virtual numeric integer_content() const;
	virtual ex smod(const numeric &xi) const;
	virtual numeric max_coefficient() const;

	// noncommutativity
	virtual unsigned return_type() const;
	virtual tinfo_t return_type_tinfo() const;

	// functions for complex expressions
	virtual ex conjugate() const;
	virtual ex real_part() const;
	virtual ex imag_part() const;

	virtual int compare(const basic & other) const;
	// functions that should be called from class ex only
	virtual int compare_same_type(const basic & other) const;
	virtual bool is_equal_same_type(const basic & other) const;

	virtual long calchash() const;

	// non-virtual functions in this class
public:
	/** Like print(), but dispatch to the specified class. Can be used by
	 *  implementations of print methods to dispatch to the method of the
	 *  superclass.
	 *
	 *  @see basic::print */
	template <class T>
	void print_dispatch(const print_context & c, unsigned level) const
	{
		print_dispatch(T::get_class_info_static(), c, level);
	}

	void print_dispatch(const registered_class_info & ri, const print_context & c, unsigned level) const;

	ex subs_one_level(const exmap & m, unsigned options) const;
	ex diff(const symbol & s, unsigned nth = 1) const;
	bool is_equal(const basic & other) const;
        bool is_evaluated() const
            { return global_hold or ((flags & status_flags::evaluated) != 0u); }
	const basic & hold() const;
        void set_epseq_from(size_t i, ex e);
	long gethash() const
	{
#ifdef GINAC_COMPARE_STATISTICS
		compare_statistics.total_gethash++;
#endif
		if (flags & status_flags::hash_calculated) {
#ifdef GINAC_COMPARE_STATISTICS
			compare_statistics.gethash_cached++;
#endif
			return hashvalue;
		} 
                return calchash();
	}

	tinfo_t tinfo() const {return tinfo_key;}

	/** Set some status_flags. */
	const basic & setflag(unsigned f) const {flags |= f; return *this;}

	/** Clear some status_flags. */
	const basic & clearflag(unsigned f) const {flags &= ~f; return *this;}

	void ensure_if_modifiable() const;

	// member variables
	tinfo_t tinfo_key;                  ///< type info
	mutable unsigned flags;             ///< of type status_flags
	mutable long hashvalue=0;         ///< hash value
};


// global variables

extern int max_recursion_level;


// convenience type checker template functions

/** Check if obj is a T, including base classes. */
template <class T>
inline bool is_a(const basic &obj)
{
	return dynamic_cast<const T *>(&obj) != nullptr;
}

/** Check if obj is a T, not including base classes. */
template <class T>
inline bool is_exactly_a(const basic & obj)
{
	return obj.tinfo() == &T::tinfo_static;
}

/** Constructs a new (class basic or derived) B object on the heap.
 *
 *  This function picks the object's ctor based on the given argument types.
 *
 *  This helps the constructor of ex from basic (or a derived class B) because
 *  then the constructor doesn't have to duplicate the object onto the heap.
 *  See ex::construct_from_basic(const basic &) for more information.
 */
template<class B, typename... Args>
inline B & dynallocate(Args &&... args)
{
	return const_cast<B &>(static_cast<const B &>((new B(std::forward<Args>(args)...))->setflag(status_flags::dynallocated)));
}
/** Constructs a new (class basic or derived) B object on the heap.
 *
 *  This function is needed for GiNaC classes which have public ctors from
 *  initializer lists of expressions (which are not a type and not captured
 *  by the variadic template version).
 */
template<class B>
inline B & dynallocate(std::initializer_list<ex> il)
{
	return const_cast<B &>(static_cast<const B &>((new B(il))->setflag(status_flags::dynallocated)));
}

} // namespace GiNaC

#endif // ndef __GINAC_BASIC_H__
