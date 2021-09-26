/** @file templates.h
 *
 *  This file lists important c++ templates to use pynac. Depending on
 *  the definition of the TEMPLATE macro, this can be used to
 *  explicitly instantiate these templates (see templates.cpp) or
 *  declare the templates as external (see extern_templates.h)
 *  
 *  Note: You cannot use typedefs in external template
 *        declarations. The usual typedef'ed name in GiNaC is
 *        indicated in the comments below
 *
 *  You should only include this file if you know what you are doing.
 * 
 *  */



TEMPLATE(class std::vector<GiNaC::ex>);            // GiNaC::exvector;
TEMPLATE(class GiNaC::container<std::vector>);     // GiNaC::exprseq;
