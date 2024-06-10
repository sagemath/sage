/** @file print.h
 *
 *  Definition of helper classes for expression output. */

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

#ifndef __GINAC_PRINT_H__
#define __GINAC_PRINT_H__

#include <iosfwd>
#include <string>
#include <memory>
#include <typeinfo>

#include "class_info.h"

namespace GiNaC {


/** This class stores information about a registered print_context class. */
class print_context_options {
public:
	print_context_options(const char *n, const char *p, unsigned i)
	 : name(n), parent_name(p), id(i) {}

	const char *get_name() const { return name; }
	const char *get_parent_name() const { return parent_name; }
	unsigned get_id() const { return id; }

private:
	const char *name;         /**< Class name. */
	const char *parent_name;  /**< Name of superclass. */
	unsigned id;              /**< ID number (assigned automatically). */
};

typedef class_info<print_context_options> print_context_class_info;


/** Flags to control the behavior of a print_context. */
class print_options {
public:
	enum {
		print_index_dimensions = 0x0001 ///< print the dimensions of indices
	};
};

extern unsigned next_print_context_id;


/** Base class for print_contexts. */
class print_context
{
	friend class function_options;
	friend class registered_class_options;
public:
	typedef void inherited;
	static const print_context_class_info &get_class_info_static()
        {
	        static print_context_options o("print_context", "void", next_print_context_id++);
	        static print_context_class_info reg_info(o);
        	return reg_info;
        }
	virtual const print_context_class_info &get_class_info() const
          { return print_context::get_class_info_static(); }
	virtual const char *class_name() const
          { return print_context::get_class_info_static().options.get_name(); }

	print_context();
	virtual print_context * duplicate() const
          { return new print_context(*this); }
public:
	print_context(std::ostream &, unsigned options = 0);
	virtual ~print_context() {}

	std::ostream & s; /**< stream to output to */
	unsigned options; /**< option flags */
};

template <class C>
class print_context_base : public virtual print_context
{
public:
	typedef class print_context inherited;
	static const print_context_class_info &get_class_info_static()
        {
	        static print_context_options o(typeid(C).name(),
                                "print_context", next_print_context_id++);
	        static print_context_class_info reg_info(o);
        	return reg_info;
        }
	virtual const print_context_class_info &get_class_info() const
            { return C::get_class_info_static(); }
	virtual const char *class_name() const
            { return C::get_class_info_static().options.get_name(); } \
	virtual print_context * duplicate() const { return new C(*dynamic_cast<const C*>(this)); }
};

/** Context for default (ginsh-parsable) output. */
class print_dflt : public print_context_base<print_dflt>
{
public:
        print_dflt();
	print_dflt(std::ostream &, unsigned options = 0);
};

/** Context for latex-parsable output. */
class print_latex : public print_context_base<print_latex>
{
public:
	print_latex();
        print_latex(std::ostream &, unsigned options = 0);
};

/** Context for python pretty-print output. */
class print_python : public print_context_base<print_python>
{
public:
        print_python();
	print_python(std::ostream &, unsigned options = 0);
};

/** Context for python-parsable output. */
class print_python_repr : public print_context_base<print_python_repr>
{
public:
        print_python_repr();
	print_python_repr(std::ostream &, unsigned options = 0);
};

/** Context for tree-like output for debugging. */
class print_tree : public print_context_base<print_tree>
{
public:
        print_tree();
	print_tree(unsigned d);
	print_tree(std::ostream &, unsigned options = 0, unsigned d = 4);

	const unsigned delta_indent; /**< size of indentation step */
};

/** Check if obj is a T, including base classes. */
template <class T>
inline bool is_a(const print_context & obj)
{ return dynamic_cast<const T *>(&obj) != nullptr; }


class basic;

/** Base class for print_functor handlers */
class print_functor_impl {
public:
	virtual ~print_functor_impl() {}
	virtual print_functor_impl *duplicate() const = 0;
	virtual void operator()(const basic & obj, const print_context & c, unsigned level) const = 0;
};

/** print_functor handler for pointer-to-functions of class T, context type C */
template <class T, class C>
class print_ptrfun_handler : public print_functor_impl {
public:
	typedef void (*F)(const T &, const C &, unsigned);

	print_ptrfun_handler(F f_) : f(f_) {}
	print_ptrfun_handler *duplicate() const override { return new print_ptrfun_handler(*this); }

	void operator()(const basic & obj, const print_context & c, unsigned level) const override
	{
		// Call the supplied function
		f(dynamic_cast<const T &>(obj), dynamic_cast<const C &>(c), level);
	}

private:
	F f;
};

/** print_functor handler for member functions of class T, context type C */
template <class T, class C>
class print_memfun_handler : public print_functor_impl {
public:
	typedef void (T::*F)(const C & c, unsigned level) const;

	print_memfun_handler(F f_) : f(f_) {}
	print_memfun_handler *duplicate() const override { return new print_memfun_handler(*this); }

	void operator()(const basic & obj, const print_context & c, unsigned level) const override
	{
		// Call the supplied member function
		return (dynamic_cast<const T &>(obj).*f)(dynamic_cast<const C &>(c), level);
	}

private:
	F f;
};

/** This class represents a print method for a certain algebraic class and
 *  print_context type. Its main purpose is to hide the difference between
 *  member functions and nonmember functions behind one unified operator()
 *  interface. print_functor has value semantics and acts as a smart pointer
 *  (with deep copy) to a class derived from print_functor_impl which
 *  implements the actual function call. */
class print_functor {
public:
	print_functor() : impl(nullptr) {}
	print_functor(const print_functor & other) : impl(other.impl.get() ? other.impl->duplicate() : nullptr) {}
	print_functor(std::unique_ptr<print_functor_impl> impl_) : impl(std::move(impl_)) {}

	template <class T, class C>
	print_functor(void f(const T &, const C &, unsigned)) : impl(new print_ptrfun_handler<T, C>(f)) {}

	template <class T, class C>
	print_functor(void (T::*f)(const C &, unsigned) const) : impl(new print_memfun_handler<T, C>(f)) {}

	print_functor & operator=(const print_functor & other)
	{
		if (this != &other) {
			print_functor_impl *p = other.impl.get();
			impl.reset(p ? other.impl->duplicate() : nullptr);
		}
		return *this;
	}

	void operator()(const basic & obj, const print_context & c, unsigned level) const
	{
		(*impl)(obj, c, level);
	}

	bool is_valid() const { return impl.get(); }

private:
	std::unique_ptr<print_functor_impl> impl;
};


} // namespace GiNaC

#endif // ndef __GINAC_BASIC_H__
