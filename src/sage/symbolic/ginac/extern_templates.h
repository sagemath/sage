/** @file extern_templates.h
 *
 *  Include this header to declare GiNaC templates external. If
 *  nothing else, this saves time and space when compiling the
 *  application that uses libpynac. 
 * 
 *  Note that you cannot use pynac with templates that were
 *  instantiated outside of pynac since RTTI does not work across DSO
 *  boundaries. See templates.cpp for more details.
 *
 *  You should always include this header in your application after
 *  including <pynac/ginac.h>.
 * 
 *  */

#ifndef PYNAC_EXTERN_TEMPLATES__H
#define PYNAC_EXTERN_TEMPLATES__H

#ifndef __GINAC_H__
#error You must #include "ginac.h" first!
#endif

#define TEMPLATE(cls) extern template cls
#include "templates.h"
#undef TEMPLATE

#endif
