/** @file archive.cpp
 *
 *  Archiving of GiNaC expressions. */

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

#include "archive.h"
#include "registrar.h"
#include "ex.h"
#include "lst.h"
#include "symbol.h"
#include "lst.h"
#include "container.h"
#include "numeric.h"
#include "constant.h"
#include "function.h"
#include "fderivative.h"
#include "matrix.h"
#include "pseries.h"
#include "add.h"
#include "mul.h"
#include "power.h"
#include "infinity.h"
#include "exprseq.h"
#include "relational.h"
#include "pynac-config.h"
#include "tostring.h"

#include <iostream>
#include <stdexcept>
#include <unordered_map>

namespace GiNaC {


void archive::archive_ex(const ex &e, const char *name)
{
	// Create root node (which recursively archives the whole expression tree)
	// and add it to the archive
	archive_node_id id = add_node(archive_node(*this, e));

	// Add root node ID to list of archived expressions
	archived_ex ae = archived_ex(atomize(name), id);
	exprs.push_back(ae);
}


/** Add archive_node to archive if the corresponding expression is
 *  not already archived.
 *  @return ID of archived node */
archive_node_id archive::add_node(const archive_node &n)
{
	// Look if expression is known to be in some node already.
	if (n.has_ex()) {
		auto i = exprtable.find(n.get_ex());
		if (i != exprtable.end())
			return i->second;
		nodes.push_back(n);
		exprtable[n.get_ex()] = nodes.size() - 1;
		return nodes.size() - 1;
	}

	// Not found, add archive_node to nodes vector
	nodes.push_back(n);
	return nodes.size()-1;
}


/** Retrieve archive_node by ID. */
archive_node &archive::get_node(archive_node_id id)
{
	if (id >= nodes.size())
		throw (std::range_error("archive::get_node(): archive node ID out of range"));

	return nodes[id];
}


ex archive::unarchive_ex(const lst &sym_lst, const char *name) const
{
	// Find root node
	std::string name_string = name;
	archive_atom id = atomize(name_string);
        for (const auto & elem : exprs) {
		if (elem.name == id) {
                        // Recursively unarchive all nodes, starting at the root node
                        lst sym_lst_copy = sym_lst;
                        return nodes[elem.root].unarchive(sym_lst_copy);
                }
	}
	throw (std::runtime_error("expression with name '" + name_string + "' not found in archive"));
}

ex archive::unarchive_ex(const lst &sym_lst, unsigned index) const
{
	if (index >= exprs.size())
		throw (std::range_error("index of archived expression out of range"));

	// Recursively unarchive all nodes, starting at the root node
	lst sym_lst_copy = sym_lst;
	return nodes[exprs[index].root].unarchive(sym_lst_copy);
}

ex archive::unarchive_ex(const lst &sym_lst, std::string &name, unsigned index) const
{
	if (index >= exprs.size())
		throw (std::range_error("index of archived expression out of range"));

	// Return expression name
	name = unatomize(exprs[index].name);

	// Recursively unarchive all nodes, starting at the root node
	lst sym_lst_copy = sym_lst;
	return nodes[exprs[index].root].unarchive(sym_lst_copy);
}

unsigned archive::num_expressions() const
{
	return exprs.size();
}

const archive_node &archive::get_top_node(unsigned index) const
{
	if (index >= exprs.size())
		throw (std::range_error("index of archived expression out of range"));

	return nodes[exprs[index].root];
}


/*
 *  Archive file format
 *
 *   - 4 bytes signature 'GARC'
 *   - unsigned version number
 *   - unsigned number of atoms
 *      - atom strings (each zero-terminated)
 *   - unsigned number of expressions
 *      - unsigned name atom
 *      - unsigned root node ID
 *   - unsigned number of nodes
 *      - unsigned number of properties
 *        - unsigned containing type (PTYPE_*) in its lower 3 bits and
 *          name atom in the upper bits
 *        - unsigned property value
 *
 *  Unsigned quantities are stored in a compressed format:
 *   - numbers in the range 0x00..0x7f are stored verbatim (1 byte)
 *   - numbers larger than 0x7f are stored in 7-bit packets (1 byte per
 *     packet), starting with the LSBs; all bytes except the last one have
 *     their upper bit set
 *
 *  Examples:
 *   0x00           =   0x00
 *    ..                 ..
 *   0x7f           =   0x7f
 *   0x80 0x01      =   0x80
 *    ..   ..            ..
 *   0xff 0x01      =   0xff
 *   0x80 0x02      =  0x100
 *    ..   ..            ..
 *   0xff 0x02      =  0x17f
 *   0x80 0x03      =  0x180
 *    ..   ..            ..
 *   0xff 0x7f      = 0x3fff
 *   0x80 0x80 0x01 = 0x4000
 *    ..   ..   ..       ..
 */

/** Write unsigned integer quantity to stream. */
static void write_unsigned(std::ostream &os, unsigned val)
{
	while (val >= 0x80) {
		os.put((val & 0x7f) | 0x80);
		val >>= 7;
	}
	os.put(val);
}

/** Read unsigned integer quantity from stream. */
static unsigned read_unsigned(std::istream &is)
{
	unsigned char b = 0;
	unsigned ret = 0;
	unsigned shift = 0;
	do {
		char b2 = 0;
		if (is.get(b2))
        		b = b2;
		ret |= (b & 0x7f) << shift;
		shift += 7;
	} while ((b & 0x80) != 0);
	return ret;
}

/** Write archive_node to binary data stream. */
std::ostream &operator<<(std::ostream &os, const archive_node &n)
{
	// Write properties
	unsigned num_props = n.props.size();
	write_unsigned(os, num_props);
	for (unsigned i=0; i<num_props; i++) {
		write_unsigned(os, n.props[i].type | (n.props[i].name << 3));
		write_unsigned(os, n.props[i].value);
	}
	return os;
}

/** Write archive to binary data stream. */
std::ostream &operator<<(std::ostream &os, const archive &ar)
{
	// Write header
	os.put('G');	// Signature
	os.put('A');
	os.put('R');
	os.put('C');
	write_unsigned(os, PYNAC_ARCHIVE_VERSION);

	// Write atoms
	unsigned num_atoms = ar.atoms.size();
	write_unsigned(os, num_atoms);
	for (unsigned i=0; i<num_atoms; i++)
		os << ar.atoms[i] << std::ends;

	// Write expressions
	unsigned num_exprs = ar.exprs.size();
	write_unsigned(os, num_exprs);
	for (unsigned i=0; i<num_exprs; i++) {
		write_unsigned(os, ar.exprs[i].name);
		write_unsigned(os, ar.exprs[i].root);
	}

	// Write nodes
	unsigned num_nodes = ar.nodes.size();
	write_unsigned(os, num_nodes);
	for (unsigned i=0; i<num_nodes; i++)
		os << ar.nodes[i];
	return os;
}

/** Read archive_node from binary data stream. */
std::istream &operator>>(std::istream &is, archive_node &n)
{
	// Read properties
	unsigned num_props = read_unsigned(is);
	n.props.resize(num_props);
	for (unsigned i=0; i<num_props; i++) {
		unsigned name_type = read_unsigned(is);
		n.props[i].type = static_cast<archive_node::property_type>(name_type & 7);
		n.props[i].name = name_type >> 3;
		n.props[i].value = read_unsigned(is);
	}
	return is;
}

/** Read archive from binary data stream. */
std::istream &operator>>(std::istream &is, archive &ar)
{
	// Read header
	char c1 = 0, c2 = 0, c3 = 0, c4 = 0;
	is.get(c1); is.get(c2); is.get(c3); is.get(c4);
	if (c1 != 'G' || c2 != 'A' || c3 != 'R' || c4 != 'C')
		throw (std::runtime_error("not a GiNaC archive (signature not found)"));
	unsigned version = read_unsigned(is);
	if (version > PYNAC_ARCHIVE_VERSION || version < PYNAC_ARCHIVE_VERSION - PYNAC_ARCHIVE_AGE)
		throw (std::runtime_error("archive version " + ToString(version) + " cannot be read by this GiNaC library (which supports versions " + ToString(PYNAC_ARCHIVE_VERSION-PYNAC_ARCHIVE_AGE) + " thru " + ToString(PYNAC_ARCHIVE_VERSION)));

	// Read atoms
	unsigned num_atoms = read_unsigned(is);
	ar.atoms.resize(num_atoms);
	for (unsigned i=0; i<num_atoms; i++) {
		getline(is, ar.atoms[i], '\0');
		ar.inverse_atoms[ar.atoms[i]] = i;
	}

	// Read expressions
	unsigned num_exprs = read_unsigned(is);
	ar.exprs.resize(num_exprs);
	for (unsigned i=0; i<num_exprs; i++) {
		archive_atom name = read_unsigned(is);
		archive_node_id root = read_unsigned(is);
		ar.exprs[i] = archive::archived_ex(name, root);
	}

	// Read nodes
	unsigned num_nodes = read_unsigned(is);
	ar.nodes.resize(num_nodes, ar);
	for (unsigned i=0; i<num_nodes; i++)
		is >> ar.nodes[i];
	return is;
}


/** Atomize a string (i.e. convert it into an ID number that uniquely
 *  represents the string). */
archive_atom archive::atomize(const std::string &s) const
{
	// Search for string in inverse_atoms map.
	inv_at_cit i = inverse_atoms.find(s);
	if (i!=inverse_atoms.end())
		return i->second;

	// Not found, add to atoms vector
	archive_atom id = atoms.size();
	atoms.push_back(s);
	inverse_atoms[s] = id;
	return id;
}

/** Unatomize a string (i.e. convert the ID number back to the string). */
const std::string &archive::unatomize(archive_atom id) const
{
	if (id >= atoms.size())
		throw (std::range_error("archive::unatomizee(): atom ID out of range"));

	return atoms[id];
}


/** Assignment operator of archive_node. */
const archive_node &archive_node::operator=(const archive_node &other)
{
	if (this != &other) {
		// archive &a member doesn't get copied
		props = other.props;
		has_expression = other.has_expression;
		e = other.e;
	}
	return *this;
}


/** Recursively construct archive node from expression. */
archive_node::archive_node(archive &ar, const ex &expr)
  : a(ar), has_expression(true), e(expr)
{
	expr.bp->archive(*this);
}


/** Check if the archive_node stores the same expression as another
 *  archive_node.
 *  @return "true" if expressions are the same */
bool archive_node::has_same_ex_as(const archive_node &other) const
{
	if (!has_expression || !other.has_expression)
		return false;
	return e.bp == other.e.bp;
}

archive_node::archive_node_cit
		archive_node::find_first(const std::string &name) const
{	
	archive_atom name_atom = a.atomize(name);
	for (auto i=props.begin(); i!=props.end(); ++i)
		if (i->name == name_atom)
			return i;
	return props.end();;
}

archive_node::archive_node_cit
		archive_node::find_last(const std::string &name) const
{
	archive_atom name_atom = a.atomize(name);
	for (auto i=props.end(); i!=props.begin();) {
		--i;
		if (i->name == name_atom)
			return i;
	}
	return props.end();
}

void archive_node::add_bool(const std::string &name, bool value)
{
	props.emplace_back(a.atomize(name), PTYPE_BOOL, static_cast<unsigned int>(value));
}

void archive_node::add_unsigned(const std::string &name, unsigned value)
{
	props.emplace_back(a.atomize(name), PTYPE_UNSIGNED, value);
}

void archive_node::add_string(const std::string &name, const std::string &value)
{
	props.emplace_back(a.atomize(name), PTYPE_STRING, a.atomize(value));
}

void archive_node::add_ex(const std::string &name, const ex &value)
{
	// Recursively create an archive_node and add its ID to the properties of this node
	archive_node_id id = a.add_node(archive_node(a, value));
	props.emplace_back(a.atomize(name), PTYPE_NODE, id);
}


bool archive_node::find_bool(const std::string &name, bool &ret, unsigned index) const
{
	archive_atom name_atom = a.atomize(name);
	unsigned found_index = 0;
	for (const auto & elem : props) {
		if (elem.type == PTYPE_BOOL && elem.name == name_atom) {
			if (found_index == index) {
				ret = (elem.value != 0u);
				return true;
			}
			found_index++;
		}
	}
	return false;
}

bool archive_node::find_unsigned(const std::string &name, unsigned &ret, unsigned index) const
{
	archive_atom name_atom = a.atomize(name);
	unsigned found_index = 0;
	for (const auto & elem : props) {
		if (elem.type == PTYPE_UNSIGNED && elem.name == name_atom) {
			if (found_index == index) {
				ret = elem.value;
				return true;
			}
			found_index++;
		}
	}
	return false;
}

bool archive_node::find_string(const std::string &name, std::string &ret, unsigned index) const
{
	archive_atom name_atom = a.atomize(name);
	unsigned found_index = 0;
        for (const auto & elem : props) {
		if (elem.type == PTYPE_STRING && elem.name == name_atom) {
			if (found_index == index) {
				ret = a.unatomize(elem.value);
				return true;
			}
			found_index++;
		}
	}
	return false;
}

void archive_node::find_ex_by_loc(archive_node_cit loc, ex &ret, lst &sym_lst)
		const
{
	ret = a.get_node(loc->value).unarchive(sym_lst);
}

bool archive_node::find_ex(const std::string &name, ex &ret, lst &sym_lst, unsigned index) const
{
	archive_atom name_atom = a.atomize(name);
	unsigned found_index = 0;
        for (const auto & elem : props) {
		if (elem.type == PTYPE_NODE && elem.name == name_atom) {
			if (found_index == index) {
				ret = a.get_node(elem.value).unarchive(sym_lst);
				return true;
			}
			found_index++;
		}
	}
	return false;
}

const archive_node &archive_node::find_ex_node(const std::string &name, unsigned index) const
{
	archive_atom name_atom = a.atomize(name);
	unsigned found_index = 0;
        for (const auto & elem : props) {
		if (elem.type == PTYPE_NODE && elem.name == name_atom) {
			if (found_index == index)
				return a.get_node(elem.value);
			found_index++;
		}
	}
	throw (std::runtime_error("property with name '" + name + "' not found in archive node"));
}


void archive_node::get_properties(propinfovector &v) const
{
	v.clear();
        for (const auto & elem : props) {
		property_type type = elem.type;
		std::string name = a.unatomize(elem.name);

		bool found = false;
                for (auto & velem : v) {
			if (velem.type == type && velem.name == name) {
				velem.count++;
				found = true;
				break;
			}
		}
		if (!found)
			v.emplace_back(type, name);
	}	
}

using unarch_func_t = ex (*)(const archive_node&, lst&);
using unarch_func_map = std::unordered_map<std::string, unarch_func_t>;

static unarch_func_map& fill_map()
{
        static unarch_func_map map;

        map["symbol"] = symbol::unarchive;
        map["lst"] = lst::unarchive;
        map["numeric"] = numeric::unarchive;
        map["constant"] = constant::unarchive;
        map["function"] = function::unarchive;
        map["fderivative"] = fderivative::unarchive;
        map["matrix"] = matrix::unarchive;
        map["pseries"] = pseries::unarchive;
        map["add"] = add::unarchive;
        map["mul"] = mul::unarchive;
        map["power"] = power::unarchive;
        map["infinity"] = infinity::unarchive;
        map["exprseq"] = exprseq::unarchive;
        map["relational"] = relational::unarchive;

        return map;
}

static const unarch_func_t& find_unarch_func(const std::string& s)
{
        static unarch_func_map& map = fill_map();
        auto it = map.find(s);
        if (it == map.end())
                throw std::runtime_error("can't happen");
        return it->second;
}

/** Convert archive node to GiNaC expression. */
ex archive_node::unarchive(lst &sym_lst) const
{
	// Already unarchived? Then return cached unarchived expression.
	if (has_expression)
		return e;

	// Find instantiation function for class specified in node
	std::string class_name;
	if (!find_string("class", class_name))
		throw (std::runtime_error("archive node contains no class name"));
	const unarch_func_t f = find_unarch_func(class_name);

	// Call instantiation function
	e = f(*this, sym_lst);
	has_expression = true;
	return ex(e);
}


void archive::clear()
{
	atoms.clear();
	inverse_atoms.clear();
	exprs.clear();
	nodes.clear();
	exprtable.clear();
}


/** Delete cached unarchived expressions in all archive_nodes (mainly for debugging). */
void archive::forget()
{
	for_each(nodes.begin(), nodes.end(), std::mem_fun_ref(&archive_node::forget));
}

/** Delete cached unarchived expressions from node (for debugging). */
void archive_node::forget()
{
	has_expression = false;
	e = 0;
}


/** Print archive to stream in ugly raw format (for debugging). */
void archive::printraw(std::ostream &os) const
{
	// Dump atoms
	os << "Atoms:\n";
	{
		archive_atom id = 0;
                for (const auto & elem : atoms) {
			os << " " << id << " " << elem << std::endl;
			id++;
		}
	}
	os << std::endl;

	// Dump expressions
	os << "Expressions:\n";
	{
		unsigned index = 0;
                for (const auto & elem : exprs) {
			os << " " << index << " \"" << unatomize(elem.name) << "\" root node " << elem.root << std::endl;
			index++;
		}
	}
	os << std::endl;

	// Dump nodes
	os << "Nodes:\n";
	{
		archive_node_id id = 0;
                for (const auto & elem : nodes) {
			os << " " << id << " ";
			elem.printraw(os);
			id++;
		}
	}
}

/** Output archive_node to stream in ugly raw format (for debugging). */
void archive_node::printraw(std::ostream &os) const
{
	// Dump cached unarchived expression
	if (has_expression)
		os << "(basic * " << e.bp << " = " << e << ")\n";
	else
		os << "\n";

	// Dump properties
        for (const auto & elem : props) {
		os << "  ";
		switch (elem.type) {
			case PTYPE_BOOL: os << "bool"; break;
			case PTYPE_UNSIGNED: os << "unsigned"; break;
			case PTYPE_STRING: os << "string"; break;
			case PTYPE_NODE: os << "node"; break;
			default: os << "<unknown>"; break;
		}
		os << " \"" << a.unatomize(elem.name) << "\" " << elem.value << std::endl;
	}
}

/** Create a dummy archive.  The intention is to fill archive_node's default
 *  ctor, which is currently a Cint-requirement. */
archive* archive_node::dummy_ar_creator()
{
	static auto  some_ar = new archive;
	return some_ar;
}


} // namespace GiNaC
