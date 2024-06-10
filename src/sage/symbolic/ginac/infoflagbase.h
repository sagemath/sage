/*
 * File:   infoflagbase.h
 * Author: Ralf Stephan <ralf@ark.in-berlin.de>
 *
 * Created on August 3, 2015, 9:31 AM
 */

#ifndef INFOFLAGBASE_H
#define	INFOFLAGBASE_H

#include <bitset>
#include <string>
#include "flags.h"

namespace GiNaC {

class infoflagbase {
public:
	infoflagbase();
        infoflagbase(const infoflagbase& other) : bits(other.bits) {}

	std::string to_string() const       { return bits.to_string(); }
	bool get(unsigned flag) const;
	void set(unsigned flag, bool value);
	void clear()                        { bits.reset(); }

private:
	static constexpr unsigned const flags[] = {
		info_flags::real, info_flags::rational,
		info_flags::integer, info_flags::crational,
		info_flags::cinteger, info_flags::positive,
		info_flags::negative, info_flags::nonnegative,
		info_flags::posint, info_flags::negint, info_flags::nonnegint,
		info_flags::even, info_flags::odd, info_flags::prime,
		info_flags::nonzero, info_flags::numeric,
		};

	static unsigned index[info_flags::relation];
	static void init_index();

	std::bitset<sizeof(flags)/sizeof(unsigned)> bits;


};

} // namespace GiNaC

#endif	/* INFOFLAGBASE_H */

