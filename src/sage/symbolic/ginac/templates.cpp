/** @file templates.cpp
 *
 *  This source file explicitly instantiates templates. This is
 *  important since RTTI does not work (reliably) across DSO
 *  boundaries. If libpynac does not contain a template instantiation
 *  but the application using libpynac does, then you will get
 *  mysterious std::bad_cast errors when you pass these instantiations
 *  into pynac. See https://github.com/sagemath/sage/issues/14780
 * 
 *  */

#ifndef PYNAC_TEMPLATES__H
#define PYNAC_TEMPLATES__H

#include "ginac.h"

#define TEMPLATE(cls) template cls
#include "templates.h"
#undef TEMPLATE

#endif

