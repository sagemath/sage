/** @file registrar.h
 *
 *  GiNaC's class registrar (for class basic and all classes derived from it). */

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

#ifndef __GINAC_REGISTRAR_H__
#define __GINAC_REGISTRAR_H__

#include "class_info.h"
#include "print.h"

#include <string>
#include <list>
#include <vector>

namespace GiNaC {

class archive_node;

template <template <class T, class = std::allocator<T> > class> class container;
typedef container<std::list> lst;

/** Definitions for the tinfo mechanism. */
typedef const void * tinfo_t;
struct tinfo_static_t {};


/** This class stores information about a registered GiNaC class. */
class registered_class_options {
public:
	registered_class_options(const char *n, const char *p, tinfo_t ti)
	 : name(n), parent_name(p), tinfo_key(ti) {}

	const char *get_name() const { return name; }
	const char *get_parent_name() const { return parent_name; }
	tinfo_t get_id() const { return tinfo_key; }
	const std::vector<print_functor> &get_print_dispatch_table() const { return print_dispatch_table; }

	template <class Ctx, class T, class C>
	registered_class_options & print_func(void f(const T &, const C & c, unsigned))
	{
		set_print_func(Ctx::get_class_info_static().options.get_id(), f);
		return *this;
	}

	template <class Ctx, class T, class C>
	registered_class_options & print_func(void (T::*f)(const C &, unsigned))
	{
		set_print_func(Ctx::get_class_info_static().options.get_id(), f);
		return *this;
	}

	template <class Ctx>
	registered_class_options & print_func(const print_functor & f)
	{
		set_print_func(Ctx::get_class_info_static().options.get_id(), f);
		return *this;
	}

	void set_print_func(unsigned id, const print_functor & f)
	{
		if (id >= print_dispatch_table.size())
			print_dispatch_table.resize(id + 1);
		print_dispatch_table[id] = f;
	}

private:
	const char *name;         /**< Class name. */
	const char *parent_name;  /**< Name of superclass. */
	tinfo_t tinfo_key;        /**< Type information key. */
	std::vector<print_functor> print_dispatch_table; /**< Method table for print() dispatch */
};

typedef class_info<registered_class_options> registered_class_info;


/** Primary macro for inclusion in the declaration of each registered class. */
#define GINAC_DECLARE_REGISTERED_CLASS_NO_CTORS(classname, supername) \
public: \
	typedef supername inherited; \
    static const GiNaC::tinfo_static_t tinfo_static; \
private: \
	static GiNaC::registered_class_info reg_info; \
public: \
	static GiNaC::registered_class_info &get_class_info_static() { return reg_info; } \
	virtual const GiNaC::registered_class_info &get_class_info() const override { return classname::get_class_info_static(); } \
	virtual GiNaC::registered_class_info &get_class_info() override { return classname::get_class_info_static(); } \
	virtual const char *class_name() const override { return classname::get_class_info_static().options.get_name(); } \
	\
	classname(const GiNaC::archive_node &n, GiNaC::lst &sym_lst); \
	virtual void archive(GiNaC::archive_node &n) const override; \
	\
	class visitor { \
	public: \
		virtual void visit(const classname &) = 0; \
		virtual ~visitor() {}; \
	};

/** Macro for inclusion in the declaration of each registered class.
 *  It declares some functions that are common to all classes derived
 *  from 'basic' as well as all required stuff for the GiNaC class
 *  registry (mainly needed for archiving). */
#define GINAC_DECLARE_REGISTERED_CLASS(classname, supername) \
	GINAC_DECLARE_REGISTERED_CLASS_NO_CTORS(classname, supername) \
public: \
	classname(); \
	virtual classname * duplicate() const override { return new classname(*this); } \
	\
	virtual void accept(GiNaC::visitor & vis) const override \
	{ \
		if (visitor *p = dynamic_cast<visitor *>(&vis)) \
			p->visit(*this); \
		else \
			inherited::accept(vis); \
	} \
protected: \
	virtual int compare_same_type(const GiNaC::basic & other) const override; \
private:


/** Macro for inclusion in the implementation of each registered class.
 *  Additional options can be specified. */
#define GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(classname, supername, options) \
	GiNaC::registered_class_info classname::reg_info = GiNaC::registered_class_info(GiNaC::registered_class_options(#classname, #supername, &classname::tinfo_static).options); \
	const GiNaC::tinfo_static_t classname::tinfo_static = {};

/** Find type information key by class name. */
extern tinfo_t find_tinfo_key(const std::string &class_name);

/** Add or replace a print method. */
template <class Alg, class Ctx, class T, class C>
extern void set_print_func(void f(const T &, const C & c, unsigned))
{
	Alg::get_class_info_static().options.set_print_func(Ctx::get_class_info_static().options.get_id(), f);
}

/** Add or replace a print method. */
template <class Alg, class Ctx, class T, class C>
extern void set_print_func(void (T::*f)(const C &, unsigned))
{
	Alg::get_class_info_static().options.set_print_func(Ctx::get_class_info_static().options.get_id(), f);
}


} // namespace GiNaC

#endif // ndef __GINAC_REGISTRAR_H__
