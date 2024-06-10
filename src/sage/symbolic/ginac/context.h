/* 
 * File:   context.h
 * Author: Ralf Stephan <ralf@ark.in-berlin.de>
 */

#ifndef CONTEXT_H
#define	CONTEXT_H

#include <string>

namespace GiNaC {

extern bool global_hold;

void set_state(const std::string&, bool);

} // namespace GiNaC

#endif	/* CONTEXT_H */

