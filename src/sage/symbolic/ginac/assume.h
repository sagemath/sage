/* 
 * File:   assume.h
 * Author: Ralf Stephan <ralf@ark.in-berlin.de>
 *
 * Created on August 5, 2015, 7:22 AM
 */

#ifndef ASSUME_H
#define	ASSUME_H

#include <string>

namespace GiNaC {

void assume(ex rel);
void assume(ex x, char* flag_desc);
void forget(ex rel);
void forget(ex x, char* flag_desc);

} // namespace GiNaC

#endif	/* ASSUME_H */

