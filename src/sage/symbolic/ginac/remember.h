/** @file remember.h
 *
 *  Interface to helper classes for using the remember option
 *  in GiNaC functions */

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

#ifndef __GINAC_REMEMBER_H__
#define __GINAC_REMEMBER_H__

#include <iosfwd>
#include <vector>
#include <list>

namespace GiNaC {

class function;
class ex;
	
/** A single entry in the remember table of a function.
 *  Needs to be a friend of class function to access 'seq'.
 *  'last_access' and 'successful_hits' are updated at each successful
 *  'is_equal'. */
class remember_table_entry {
public:
	remember_table_entry(function const & f, ex  r);
	bool is_equal(function const & f) const;
	ex get_result() const { return result; }
	unsigned long get_last_access() const { return last_access; }
	unsigned long get_successful_hits() const { return successful_hits; };

protected:
	long hashvalue;
	exvector seq;
	ex result;
	mutable unsigned long last_access;
	mutable unsigned successful_hits;
	static unsigned long access_counter;
};    

/** A list of entries in the remember table having some least
 *  significant bits of the hashvalue in common. */
class remember_table_list : public std::list<remember_table_entry> {
public:
	remember_table_list(unsigned as, unsigned strat);
	void add_entry(function const & f, ex const & result);
	bool lookup_entry(function const & f, ex & result) const;
protected:
	unsigned max_assoc_size;
	unsigned remember_strategy;
};

/** The remember table is organized like an n-fold associative cache
 *  in a microprocessor.  The table has a width of 's' (which is rounded
 *  to table_size, some power of 2 near 's', internally) and a depth of 'as'
 *  (unless you choose that entries are never discarded). The place where
 *  an entry is stored depends on the hashvalue of the parameters of the
 *  function (this corresponds to the address of byte to be cached).
 *  The 'log_2(table_size)' least significant bits of this hashvalue
 *  give the slot in which the entry will be stored or looked up.
 *  Each slot can take up to 'as' entries. If a slot is full, an older
 *  entry is removed by one of the following strategies:
 *   - oldest entry (the first one in the list)
 *   - least recently used (the one with the lowest 'last_access')
 *   - least frequently used (the one with the lowest 'successful_hits')
 *  or all entries are kept which means that the table grows indefinitely. */
class remember_table : public std::vector<remember_table_list> {
public:
	remember_table();
	remember_table(unsigned s, unsigned as, unsigned strat);
	bool lookup_entry(function const & f, ex & result) const;
	void add_entry(function const & f, ex const & result);
	void clear_all_entries();
	void show_statistics(std::ostream & os, unsigned level) const;
	static std::vector<remember_table> & remember_tables();
protected:
	void init_table();
	unsigned table_size;
	unsigned max_assoc_size;
	unsigned remember_strategy;
};      

} // namespace GiNaC

#endif // ndef __GINAC_REMEMBER_H__
