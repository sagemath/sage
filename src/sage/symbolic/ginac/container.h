/** @file container.h
 *
 *  Wrapper template for making GiNaC classes out of STL containers. */

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

#ifndef __GINAC_CONTAINER_H__
#define __GINAC_CONTAINER_H__

#include "ex.h"
#include "print.h"
#include "archive.h"
#include "assertion.h"

#include <iterator>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <list>
#include <memory>
#include <algorithm>

namespace GiNaC {


/** Helper template for encapsulating the reserve() mechanics of STL containers. */
template <template <class T, class = std::allocator<T> > class C>
class container_storage {
protected:
	typedef C<ex> STLT;

	container_storage() {}
	container_storage(size_t n, const ex & e) : seq(n, e) {}

	template <class In>
	container_storage(In b, In e) : seq(b, e) {}

	void reserve(size_t) {}
	static void reserve(STLT &, size_t) {}

	STLT seq;

	// disallow destruction of container through a container_storage*
protected:
	~container_storage() {}
};

template <>
inline void container_storage<std::vector>::reserve(size_t n) { seq.reserve(n); }

template <>
inline void container_storage<std::vector>::reserve(std::vector<ex> & v, size_t n) { v.reserve(n); }


/** Helper template to allow initialization of containers via an overloaded
 *  comma operator (idea stolen from Blitz++). */
template <typename T, typename STLT>
class container_init {
public:
	container_init(STLT & s) : stlt(s) {}

	container_init<T, STLT> operator,(const T & x)
	{
		stlt.push_back(x);
		return container_init<T, STLT>(stlt);
	}

	// The following specializations produce much tighter code than the
	// general case above

	container_init<T, STLT> operator,(int x)
	{
		stlt.push_back(x);
		return container_init<T, STLT>(stlt);
	}

	container_init<T, STLT> operator,(unsigned int x)
	{
		stlt.push_back(x);
		return container_init<T, STLT>(stlt);
	}

	container_init<T, STLT> operator,(long x)
	{
		stlt.push_back(x);
		return container_init<T, STLT>(stlt);
	}

	container_init<T, STLT> operator,(unsigned long x)
	{
		stlt.push_back(x);
		return container_init<T, STLT>(stlt);
	}

	container_init<T, STLT> operator,(double x)
	{
		stlt.push_back(x);
		return container_init<T, STLT>(stlt);
	}

	container_init<T, STLT> operator,(const symbol & x)
	{
		stlt.push_back(T(x));
		return container_init<T, STLT>(stlt);
	}

private:
	container_init();
	STLT & stlt;
};

/** Wrapper template for making GiNaC classes out of STL containers. */
template <template <class T, class = std::allocator<T> > class C>
class container : public basic, public container_storage<C> {
	GINAC_DECLARE_REGISTERED_CLASS(container, basic)

	friend class print_order;

protected:
	typedef typename container_storage<C>::STLT STLT;

public:
	typedef typename STLT::const_iterator const_iterator;
	typedef typename STLT::const_reverse_iterator const_reverse_iterator;

protected:
	// helpers
	static tinfo_t get_tinfo() { return nullptr; }
	static unsigned get_default_flags() { return 0; }
	static const char* get_open_delim() { return "("; }
	static const char* get_close_delim() { return ")"; }

	// constructors
public:
	container(STLT const & s, bool discardable = false) : inherited(get_tinfo())
	{
		setflag(get_default_flags());

		if (discardable)
			this->seq.swap(const_cast<STLT &>(s));
		else
			this->seq = s;
	}

	explicit container(std::unique_ptr<STLT> vp) : inherited(get_tinfo())
	{
		setflag(get_default_flags());
		this->seq.swap(*vp);
	}

	container(exvector::const_iterator b, exvector::const_iterator e)
	 : inherited(get_tinfo()), container_storage<C>(b, e)
	{
		setflag(get_default_flags());
	}

	explicit container(const ex & p1)
	 : inherited(get_tinfo()), container_storage<C>(1, p1)
	{
		setflag(get_default_flags());
	}

	container(const ex & p1, const ex & p2) : inherited(get_tinfo())
	{
		setflag(get_default_flags());
		this->reserve(this->seq, 2);
		this->seq.push_back(p1); this->seq.push_back(p2);
	}

	container(const ex & p1, const ex & p2, const ex & p3) : inherited(get_tinfo())
	{
		setflag(get_default_flags());
		this->reserve(this->seq, 3);
		this->seq.push_back(p1); this->seq.push_back(p2); this->seq.push_back(p3);
	}

	container(const ex & p1, const ex & p2, const ex & p3,
	          const ex & p4) : inherited(get_tinfo())
	{
		setflag(get_default_flags());
		this->reserve(this->seq, 4);
		this->seq.push_back(p1); this->seq.push_back(p2); this->seq.push_back(p3);
		this->seq.push_back(p4);
	}

	container(const ex & p1, const ex & p2, const ex & p3,
	          const ex & p4, const ex & p5) : inherited(get_tinfo())
	{
		setflag(get_default_flags());
		this->reserve(this->seq, 5);
		this->seq.push_back(p1); this->seq.push_back(p2); this->seq.push_back(p3);
		this->seq.push_back(p4); this->seq.push_back(p5);
	}

	container(const ex & p1, const ex & p2, const ex & p3,
	          const ex & p4, const ex & p5, const ex & p6) : inherited(get_tinfo())
	{
		setflag(get_default_flags());
		this->reserve(this->seq, 6);
		this->seq.push_back(p1); this->seq.push_back(p2); this->seq.push_back(p3);
		this->seq.push_back(p4); this->seq.push_back(p5); this->seq.push_back(p6);
	}

	container(const ex & p1, const ex & p2, const ex & p3,
	          const ex & p4, const ex & p5, const ex & p6,
	          const ex & p7) : inherited(get_tinfo())
	{
		setflag(get_default_flags());
		this->reserve(this->seq, 7);
		this->seq.push_back(p1); this->seq.push_back(p2); this->seq.push_back(p3);
		this->seq.push_back(p4); this->seq.push_back(p5); this->seq.push_back(p6);
		this->seq.push_back(p7);
	}

	container(const ex & p1, const ex & p2, const ex & p3,
	          const ex & p4, const ex & p5, const ex & p6,
	          const ex & p7, const ex & p8) : inherited(get_tinfo())
	{
		setflag(get_default_flags());
		this->reserve(this->seq, 8);
		this->seq.push_back(p1); this->seq.push_back(p2); this->seq.push_back(p3);
		this->seq.push_back(p4); this->seq.push_back(p5); this->seq.push_back(p6);
		this->seq.push_back(p7); this->seq.push_back(p8);
	}

	container(const ex & p1, const ex & p2, const ex & p3,
	          const ex & p4, const ex & p5, const ex & p6,
	          const ex & p7, const ex & p8, const ex & p9) : inherited(get_tinfo())
	{
		setflag(get_default_flags());
		this->reserve(this->seq, 9);
		this->seq.push_back(p1); this->seq.push_back(p2); this->seq.push_back(p3);
		this->seq.push_back(p4); this->seq.push_back(p5); this->seq.push_back(p6);
		this->seq.push_back(p7); this->seq.push_back(p8); this->seq.push_back(p9);
	}

	container(const ex & p1, const ex & p2, const ex & p3,
	          const ex & p4, const ex & p5, const ex & p6,
	          const ex & p7, const ex & p8, const ex & p9,
	          const ex & p10) : inherited(get_tinfo())
	{
		setflag(get_default_flags());
		this->reserve(this->seq, 10);
		this->seq.push_back(p1); this->seq.push_back(p2); this->seq.push_back(p3);
		this->seq.push_back(p4); this->seq.push_back(p5); this->seq.push_back(p6);
		this->seq.push_back(p7); this->seq.push_back(p8); this->seq.push_back(p9);
		this->seq.push_back(p10);
	}

	container(const ex & p1, const ex & p2, const ex & p3,
	          const ex & p4, const ex & p5, const ex & p6,
	          const ex & p7, const ex & p8, const ex & p9,
	          const ex & p10, const ex & p11) : inherited(get_tinfo())
	{
		setflag(get_default_flags());
		this->reserve(this->seq, 11);
		this->seq.push_back(p1); this->seq.push_back(p2); this->seq.push_back(p3);
		this->seq.push_back(p4); this->seq.push_back(p5); this->seq.push_back(p6);
		this->seq.push_back(p7); this->seq.push_back(p8); this->seq.push_back(p9);
		this->seq.push_back(p10); this->seq.push_back(p11);
	}

	container(const ex & p1, const ex & p2, const ex & p3,
	          const ex & p4, const ex & p5, const ex & p6,
	          const ex & p7, const ex & p8, const ex & p9,
	          const ex & p10, const ex & p11, const ex & p12) : inherited(get_tinfo())
	{
		setflag(get_default_flags());
		this->reserve(this->seq, 12);
		this->seq.push_back(p1); this->seq.push_back(p2); this->seq.push_back(p3);
		this->seq.push_back(p4); this->seq.push_back(p5); this->seq.push_back(p6);
		this->seq.push_back(p7); this->seq.push_back(p8); this->seq.push_back(p9);
		this->seq.push_back(p10); this->seq.push_back(p11); this->seq.push_back(p12);
	}

	container(const ex & p1, const ex & p2, const ex & p3,
	          const ex & p4, const ex & p5, const ex & p6,
	          const ex & p7, const ex & p8, const ex & p9,
	          const ex & p10, const ex & p11, const ex & p12,
	          const ex & p13) : inherited(get_tinfo())
	{
		setflag(get_default_flags());
		this->reserve(this->seq, 13);
		this->seq.push_back(p1); this->seq.push_back(p2); this->seq.push_back(p3);
		this->seq.push_back(p4); this->seq.push_back(p5); this->seq.push_back(p6);
		this->seq.push_back(p7); this->seq.push_back(p8); this->seq.push_back(p9);
		this->seq.push_back(p10); this->seq.push_back(p11); this->seq.push_back(p12);
		this->seq.push_back(p13);
	}

	container(const ex & p1, const ex & p2, const ex & p3,
	          const ex & p4, const ex & p5, const ex & p6,
	          const ex & p7, const ex & p8, const ex & p9,
	          const ex & p10, const ex & p11, const ex & p12,
	          const ex & p13, const ex & p14) : inherited(get_tinfo())
	{
		setflag(get_default_flags());
		this->reserve(this->seq, 14);
		this->seq.push_back(p1); this->seq.push_back(p2); this->seq.push_back(p3);
		this->seq.push_back(p4); this->seq.push_back(p5); this->seq.push_back(p6);
		this->seq.push_back(p7); this->seq.push_back(p8); this->seq.push_back(p9);
		this->seq.push_back(p10); this->seq.push_back(p11); this->seq.push_back(p12);
		this->seq.push_back(p13); this->seq.push_back(p14);
	}

	container(const ex & p1, const ex & p2, const ex & p3,
	          const ex & p4, const ex & p5, const ex & p6,
	          const ex & p7, const ex & p8, const ex & p9,
	          const ex & p10, const ex & p11, const ex & p12,
	          const ex & p13, const ex & p14, const ex & p15) : inherited(get_tinfo())
	{
		setflag(get_default_flags());
		this->reserve(this->seq, 15);
		this->seq.push_back(p1); this->seq.push_back(p2); this->seq.push_back(p3);
		this->seq.push_back(p4); this->seq.push_back(p5); this->seq.push_back(p6);
		this->seq.push_back(p7); this->seq.push_back(p8); this->seq.push_back(p9);
		this->seq.push_back(p10); this->seq.push_back(p11); this->seq.push_back(p12);
		this->seq.push_back(p13); this->seq.push_back(p14); this->seq.push_back(p15);
	}

	container(const ex & p1, const ex & p2, const ex & p3,
	          const ex & p4, const ex & p5, const ex & p6,
	          const ex & p7, const ex & p8, const ex & p9,
	          const ex & p10, const ex & p11, const ex & p12,
	          const ex & p13, const ex & p14, const ex & p15,
	          const ex & p16) : inherited(get_tinfo())
	{
		setflag(get_default_flags());
		this->reserve(this->seq, 16);
		this->seq.push_back(p1); this->seq.push_back(p2); this->seq.push_back(p3);
		this->seq.push_back(p4); this->seq.push_back(p5); this->seq.push_back(p6);
		this->seq.push_back(p7); this->seq.push_back(p8); this->seq.push_back(p9);
		this->seq.push_back(p10); this->seq.push_back(p11); this->seq.push_back(p12);
		this->seq.push_back(p13); this->seq.push_back(p14); this->seq.push_back(p15);
		this->seq.push_back(p16);
	}

	// First step of initialization of container with a comma-separated
	// sequence of expressions. Subsequent steps are handled by
	// container_init<>::operator,().
	container_init<ex, STLT> operator=(const ex & x)
	{
		this->seq.push_back(x);
		return container_init<ex, STLT>(this->seq);
	}

	// functions overriding virtual functions from base classes
public:
	bool info(unsigned inf) const override { return inherited::info(inf); }
	unsigned precedence() const override { return 10; }
	size_t nops() const override { return this->seq.size(); }
	const ex op(size_t i) const override;
	ex & let_op(size_t i) override;
	ex eval(int level = 0) const override;
	ex subs(const exmap & m, unsigned options = 0) const override;
        bool match(const ex & pattern, exmap& map) const override;

protected:
	ex conjugate() const override
	{
		STLT *newcont = nullptr;
		for (const auto & elem : this->seq) {
			const ex conj = elem.conjugate();
			if (newcont) {
				newcont->push_back(conj);
				continue;
			}
			if (are_ex_trivially_equal(conj, elem)) {
				continue;
			}
			newcont = new STLT;
			this->reserve(*newcont, this->seq.size());
			for (const auto & elem2 : this->seq) {
				if (&elem2 != &elem)
					newcont->push_back(elem2);
			}
			newcont->push_back(conj);
		}
		if (newcont) {
			ex result = thiscontainer(*newcont);
			delete newcont;
			return result;
		}
		return *this;
	}

	ex real_part() const override
	{
		STLT cont;
		this->reserve(cont, nops());
		for (const auto & elem : this->seq)
			cont.push_back(elem.real_part());
		return thiscontainer(cont);
	}

	ex imag_part() const override
	{
		STLT cont;
		this->reserve(cont, nops());
		for (const auto & elem : this->seq)
			cont.push_back(elem.imag_part());
		return thiscontainer(cont);
	}

	bool is_equal_same_type(const basic & other) const override;

	// new virtual functions which can be overridden by derived classes
protected:
	/** Similar to duplicate(), but with a preset sequence. Must be
	 *  overridden by derived classes. */
	virtual ex thiscontainer(const STLT & v) const { return container(v); }

	/** Similar to duplicate(), but with a preset sequence (which gets
	 *  deleted). Must be overridden by derived classes. */
	virtual ex thiscontainer(std::unique_ptr<STLT> vp) const { return container(std::move(vp)); }

	virtual void printseq(const print_context & c, const char* openbracket,
			char delim, const char* closebracket,
			unsigned this_precedence,
			unsigned upper_precedence = 0) const;

	// non-virtual functions in this class
private:
	void sort_()
	{
		std::sort(this->seq.begin(), this->seq.end(), ex_is_less());
	}

	void unique_()
	{
		auto p = std::unique(this->seq.begin(), this->seq.end(), ex_is_equal());
		this->seq.erase(p, this->seq.end());
	}

public:
	container & prepend(const ex & b);
	container & append(const ex & b);
	container & remove_first();
	container & remove_last();
	container & remove_all();
	container & sort();
	container & unique();
        static ex unarchive(const archive_node &n, lst &sym_lst);

	const_iterator begin() const {return this->seq.begin();}
	const_iterator end() const {return this->seq.end();}
	const_reverse_iterator rbegin() const {return this->seq.rbegin();}
	const_reverse_iterator rend() const {return this->seq.rend();}

protected:
	void do_print(const print_context & c, unsigned level) const override;
	void do_print_tree(const print_tree & c, unsigned level) const override;
	void do_print_python(const print_python & c, unsigned level) const;
	void do_print_python_repr(const print_python_repr & c, unsigned level) const override;
	STLT evalchildren(int level) const;
	std::unique_ptr<STLT> subschildren(const exmap & m, unsigned options = 0) const;
};

/** Default constructor */
template <template <class T, class = std::allocator<T> > class C>
container<C>::container() : inherited(get_tinfo())
{
	setflag(get_default_flags());
}

/** Construct object from archive_node. */
template <template <class T, class = std::allocator<T> > class C>
container<C>::container(const archive_node &n, lst &sym_lst) : inherited(n, sym_lst)
{
	setflag(get_default_flags());

	auto first = n.find_first("seq");
	auto last = n.find_last("seq");
	++last;
	this->reserve(this->seq, last - first);
	for (auto i=first; i<last; ++i) {
		ex e;
		n.find_ex_by_loc(i, e, sym_lst);
		this->seq.push_back(e);
	}
}

/** Unarchive the object. */
template <template <class T, class = std::allocator<T> > class C>
ex container<C>::unarchive(const archive_node &n, lst &sym_lst)
{
	return (new container(n, sym_lst))->setflag(status_flags::dynallocated);
}

/** Archive the object. */
template <template <class T, class = std::allocator<T> > class C>
void container<C>::archive(archive_node &n) const
{
	inherited::archive(n);
	for (const auto & elem : this->seq)
		n.add_ex("seq", elem);
}

template <template <class T, class = std::allocator<T> > class C>
void container<C>::do_print(const print_context & c, unsigned) const
{
	// always print brackets around seq, ignore upper_precedence
	printseq(c, get_open_delim(), ',', get_close_delim(), precedence(), precedence()+1);
}

template <template <class T, class = std::allocator<T> > class C>
void container<C>::do_print_tree(const print_tree & c, unsigned level) const
{
	c.s << std::string(level, ' ') << class_name() << " @" << this
	    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
	    << ", nops=" << nops()
	    << std::endl;
	for (const auto & elem : this->seq) {
		elem.print(c, level + c.delta_indent);
	}
	c.s << std::string(level + c.delta_indent,' ') << "=====" << std::endl;
}

template <template <class T, class = std::allocator<T> > class C>
void container<C>::do_print_python(const print_python & c, unsigned level) const
{
	printseq(c, "[", ',', "]", precedence(), precedence()+1);
}

template <template <class T, class = std::allocator<T> > class C>
void container<C>::do_print_python_repr(const print_python_repr & c, unsigned) const
{
	c.s << class_name();
	printseq(c, "(", ',', ")", precedence(), precedence()+1);
}

template <template <class T, class = std::allocator<T> > class C>
const ex container<C>::op(size_t i) const
{
	GINAC_ASSERT(i < nops());

	auto it = this->seq.begin();
	advance(it, i);
	return *it;
}

template <template <class T, class = std::allocator<T> > class C>
ex & container<C>::let_op(size_t i)
{
	GINAC_ASSERT(i < nops());

	ensure_if_modifiable();
	auto it = this->seq.begin();
	advance(it, i);
	return *it;
}

template <template <class T, class = std::allocator<T> > class C>
ex container<C>::eval(int level) const
{
	if (level == 1)
		return hold();

        return thiscontainer(evalchildren(level));
}

template <template <class T, class = std::allocator<T> > class C>
ex container<C>::subs(const exmap & m, unsigned options) const
{
	// After having subs'ed all children, this method subs'es one final
	// level, but only if the intermediate result is a container! This is
	// because if the intermediate result has eval'ed to a non-container a
	// last level substitution would be wrong, as this example involving a
	// function f and its inverse f^-1 shows:
	// f(x).subs(x==f^-1(x))
	//   -> f(f^-1(x))  [subschildren]
	//   -> x           [eval]   /* must not subs(x==f^-1(x))! */
	std::unique_ptr<STLT> vp = subschildren(m, options);
	if (vp) {
		ex result(thiscontainer(std::move(vp)));
		if (is_a<container<C> >(result))
			return ex_to<basic>(result).subs_one_level(m, options);
		return result;
	} else {
		if (is_a<container<C> >(*this))
			return subs_one_level(m, options);		
		return *this;
	}
}

/** Compare two containers of the same type. */
template <template <class T, class = std::allocator<T> > class C>
int container<C>::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<container>(other));
	const container & o = static_cast<const container &>(other);

	auto it1 = this->seq.begin(), it1end = this->seq.end(),
	               it2 = o.seq.begin(), it2end = o.seq.end();

	while (it1 != it1end && it2 != it2end) {
		int cmpval = it1->compare(*it2);
		if (cmpval)
			return cmpval;
		++it1; ++it2;
	}

	return (it1 == it1end) ? (it2 == it2end ? 0 : -1) : 1;
}

template <template <class T, class = std::allocator<T> > class C>
bool container<C>::is_equal_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<container>(other));
	const container & o = static_cast<const container &>(other);

	if (this->seq.size() != o.seq.size())
		return false;

        auto it1 = this->seq.begin(), it1end = this->seq.end(), it2 = o.seq.begin();
        while (it1 != it1end) {
                if (!it1->is_equal(*it2))
                        return false;
                ++it1; ++it2;
       }

       return true;
}

/** Add element at front. */
template <template <class T, class = std::allocator<T> > class C>
container<C> & container<C>::prepend(const ex & b)
{
	ensure_if_modifiable();
	this->seq.push_front(b);
	return *this;
}

/** Specialization of container::prepend() for std::vector. */
template<> inline container<std::vector> & container<std::vector>::prepend(const ex & b)
{
	ensure_if_modifiable();
	this->seq.insert(this->seq.begin(), b);
	return *this;
}

/** Add element at back. */
template <template <class T, class = std::allocator<T> > class C>
container<C> & container<C>::append(const ex & b)
{
	ensure_if_modifiable();
	this->seq.push_back(b);
	return *this;
}

/** Remove first element. */
template <template <class T, class = std::allocator<T> > class C>
container<C> & container<C>::remove_first()
{
	ensure_if_modifiable();
	this->seq.pop_front();
	return *this;
}

/** Specialization of container::remove_first() for std::vector. */
template<> inline container<std::vector> & container<std::vector>::remove_first()
{
	ensure_if_modifiable();
	this->seq.erase(this->seq.begin());
	return *this;
}

/** Remove last element. */
template <template <class T, class = std::allocator<T> > class C>
container<C> & container<C>::remove_last()
{
	ensure_if_modifiable();
	this->seq.pop_back();
	return *this;
}

/** Remove all elements. */
template <template <class T, class = std::allocator<T> > class C>
container<C> & container<C>::remove_all()
{
	ensure_if_modifiable();
	this->seq.clear();
	return *this;
}

/** Sort elements. */
template <template <class T, class = std::allocator<T> > class C>
container<C> & container<C>::sort()
{
	ensure_if_modifiable();
	sort_();
	return *this;
}

/** Specialization of container::sort_() for std::list. */
template<> inline void container<std::list>::sort_()
{
	this->seq.sort(ex_is_less());
}

/** Specialization of container::unique_() for std::list. */
template<> inline void container<std::list>::unique_()
{
	this->seq.unique(ex_is_equal());
}

/** Remove adjacent duplicate elements. */
template <template <class T, class = std::allocator<T> > class C>
container<C> & container<C>::unique()
{
	ensure_if_modifiable();
	unique_();
	return *this;
}

/** Print sequence of contained elements. */
template <template <class T, class = std::allocator<T> > class C>
void container<C>::printseq(const print_context & c, const char* openbracket,
		char delim, const char* closebracket, unsigned this_precedence,
		unsigned upper_precedence) const
{
	if (this_precedence <= upper_precedence)
		c.s << openbracket;

	if (!this->seq.empty()) {
		auto it = this->seq.begin(), itend = this->seq.end();
		--itend;
		while (it != itend) {
			it->print(c, this_precedence);
			c.s << delim << ' ';
			++it;
		}
		it->print(c, this_precedence);
	}

	if (this_precedence <= upper_precedence)
		c.s << closebracket;
}

template <template <class T, class = std::allocator<T> > class C>
typename container<C>::STLT container<C>::evalchildren(int level) const
{
	if (level == 1)
		return this->seq;
	if (level == -max_recursion_level)
		throw std::runtime_error("max recursion level reached");

	STLT s;
	this->reserve(s, this->seq.size());

	--level;
	for (const auto & elem : this->seq) {
		s.push_back(elem.eval(level));
	}

	return s;
}

template <template <class T, class = std::allocator<T> > class C>
std::unique_ptr<typename container<C>::STLT> container<C>::subschildren(const exmap & m, unsigned options) const
{
	// returns a NULL pointer if nothing had to be substituted
	// returns a pointer to a newly created STLT otherwise
	// (and relinquishes responsibility for the STLT)

	auto cit = this->seq.begin(), cend = this->seq.end();
	while (cit != cend) {
		const ex & subsed_ex = cit->subs(m, options);
		if (!are_ex_trivially_equal(*cit, subsed_ex)) {

			// copy first part of seq which hasn't changed
			std::unique_ptr<STLT> s(new STLT(this->seq.begin(), cit));
			this->reserve(*s, this->seq.size());

			// insert changed element
			s->push_back(subsed_ex);
			++cit;

			// copy rest
			while (cit != cend) {
				s->push_back(cit->subs(m, options));
				++cit;
			}

			return s;
		}

		++cit;
	}
	
	return std::unique_ptr<STLT>(nullptr); // nothing has changed
}

typedef container<std::vector> exprseq;
template<> const tinfo_static_t exprseq::tinfo_static;
template<> registered_class_info exprseq::reg_info;

typedef container<std::list> lst;
template<> const tinfo_static_t lst::tinfo_static;
template<> registered_class_info lst::reg_info;

template<> bool exprseq::match(const ex & pattern, exmap& map) const;
template<> bool lst::match(const ex & pattern, exmap& map) const;
} // namespace GiNaC

#endif // ndef __GINAC_CONTAINER_H__
