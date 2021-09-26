/** @file function.h
 *
 *  Interface to class of symbolic functions. */

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

#ifndef __GINAC_FUNCTION_H__
#define __GINAC_FUNCTION_H__

#include "exprseq.h"
#include "infoflagbase.h"

#include <string>
#include <vector>

// CINT needs <algorithm> to work properly with <vector>
#include <algorithm>

// the following lines have been generated for max. 14 parameters
#define DECLARE_FUNCTION_1P(NAME) \
class NAME##_SERIAL { public: static unsigned serial; }; \
const unsigned NAME##_NPARAMS = 1; \
template<typename T1> const GiNaC::function NAME(const T1 & p1) { \
	return GiNaC::function(NAME##_SERIAL::serial, GiNaC::ex(p1)); \
}

#define DECLARE_FUNCTION_2P(NAME) \
class NAME##_SERIAL { public: static unsigned serial; }; \
const unsigned NAME##_NPARAMS = 2; \
template<typename T1, typename T2> const GiNaC::function NAME(const T1 & p1, const T2 & p2) { \
	return GiNaC::function(NAME##_SERIAL::serial, GiNaC::ex(p1), GiNaC::ex(p2)); \
}

#define DECLARE_FUNCTION_3P(NAME) \
class NAME##_SERIAL { public: static unsigned serial; }; \
const unsigned NAME##_NPARAMS = 3; \
template<typename T1, typename T2, typename T3> const GiNaC::function NAME(const T1 & p1, const T2 & p2, const T3 & p3) { \
	return GiNaC::function(NAME##_SERIAL::serial, GiNaC::ex(p1), GiNaC::ex(p2), GiNaC::ex(p3)); \
}

#define DECLARE_FUNCTION_6P(NAME) \
class NAME##_SERIAL { public: static unsigned serial; }; \
const unsigned NAME##_NPARAMS = 6; \
template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6> const GiNaC::function NAME(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4, const T5 & p5, const T6 & p6) { \
	return GiNaC::function(NAME##_SERIAL::serial, GiNaC::ex(p1), GiNaC::ex(p2), GiNaC::ex(p3), GiNaC::ex(p4), GiNaC::ex(p5), GiNaC::ex(p6)); \
}
// end of generated lines

#define REGISTER_FUNCTION(NAME,OPT) \
unsigned NAME##_SERIAL::serial = \
	GiNaC::function::register_new(GiNaC::function_options(#NAME, NAME##_NPARAMS).OPT);

namespace GiNaC {

class function;

typedef ex (* eval_funcp)();
typedef ex (* evalf_funcp)(PyObject* parent);
typedef ex (* conjugate_funcp)();
typedef ex (* real_part_funcp)();
typedef ex (* imag_part_funcp)();
typedef ex (* derivative_funcp)();
typedef ex (* expl_derivative_funcp)();
typedef ex (* power_funcp)();
typedef ex (* series_funcp)();
typedef ex (* subs_funcp)();
typedef void (* print_funcp)();

// the following lines have been generated for max. 14 parameters
typedef ex (* eval_funcp_1)(const ex &);
typedef ex (* eval_funcp_2)(const ex &, const ex &);
typedef ex (* eval_funcp_3)(const ex &, const ex &, const ex &);
typedef ex (* eval_funcp_6)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* evalf_funcp_1)(const ex &, PyObject* parent);
typedef ex (* evalf_funcp_2)(const ex &, const ex &, PyObject* parent);
typedef ex (* evalf_funcp_3)(const ex &, const ex &, const ex &, PyObject* parent);
typedef ex (* evalf_funcp_6)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, PyObject* parent);

typedef ex (* conjugate_funcp_1)(const ex &);
typedef ex (* conjugate_funcp_2)(const ex &, const ex &);
typedef ex (* conjugate_funcp_3)(const ex &, const ex &, const ex &);

typedef ex (* real_part_funcp_1)(const ex &);
typedef ex (* real_part_funcp_2)(const ex &, const ex &);
typedef ex (* real_part_funcp_3)(const ex &, const ex &, const ex &);

typedef ex (* imag_part_funcp_1)(const ex &);
typedef ex (* imag_part_funcp_2)(const ex &, const ex &);
typedef ex (* imag_part_funcp_3)(const ex &, const ex &, const ex &);
typedef ex (* derivative_funcp_1)(const ex &, unsigned);
typedef ex (* derivative_funcp_2)(const ex &, const ex &, unsigned);
typedef ex (* derivative_funcp_3)(const ex &, const ex &, const ex &, unsigned);
typedef ex (* derivative_funcp_6)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, unsigned);
typedef ex (* expl_derivative_funcp_1)(const ex &, const symbol &);
typedef ex (* expl_derivative_funcp_2)(const ex &, const ex &, const symbol &);
typedef ex (* expl_derivative_funcp_3)(const ex &, const ex &, const ex &, const symbol &);
typedef ex (* power_funcp_1)(const ex &, const ex &);
typedef ex (* power_funcp_2)(const ex &, const ex &, const ex &);
typedef ex (* power_funcp_3)(const ex &, const ex &, const ex &, const ex &);
typedef ex (* series_funcp_1)(const ex &, const relational &, int, unsigned);
typedef ex (* series_funcp_2)(const ex &, const ex &, const relational &, int, unsigned);
typedef ex (* series_funcp_3)(const ex &, const ex &, const ex &, const relational &, int, unsigned);
typedef ex (* subs_funcp_1)(const exmap&, const ex &);
typedef ex (* subs_funcp_2)(const exmap&, const ex &, const ex &);
typedef ex (* subs_funcp_3)(const exmap&, const ex &, const ex &, const ex &);

typedef void (* print_funcp_1)(const ex &, const print_context &);
typedef void (* print_funcp_2)(const ex &, const ex &, const print_context &);
typedef void (* print_funcp_3)(const ex &, const ex &, const ex &, const print_context &);
// end of generated lines

// Alternatively, an exvector may be passed into the static function, instead
// of individual ex objects.  Then, the number of arguments is not limited.
typedef ex (* eval_funcp_exvector)(const exvector &);
typedef ex (* evalf_funcp_exvector)(const exvector &, PyObject* parent);
typedef ex (* conjugate_funcp_exvector)(const exvector &);
typedef ex (* real_part_funcp_exvector)(const exvector &);
typedef ex (* imag_part_funcp_exvector)(const exvector &);
typedef ex (* derivative_funcp_exvector)(const exvector &, unsigned);
typedef ex (* expl_derivative_funcp_exvector)(const exvector &, const symbol &);
typedef ex (* power_funcp_exvector)(const exvector &, const ex &);
typedef ex (* series_funcp_exvector)(const exvector &, const relational &, int, unsigned);
typedef void (* print_funcp_exvector)(const exvector &, const print_context &);


typedef ex (* derivative_funcp_exvector_symbol)(const exvector &,
		const symbol &);

class function_options
{
	friend class function;
	friend class fderivative;
public:
	function_options();
	function_options(std::string const & n, std::string const & tn=std::string());
	function_options(std::string const & n, unsigned np);
	~function_options();
	void initialize();

	function_options & dummy() { return *this; }
	function_options & set_name(std::string const & n, std::string const & tn=std::string());
	function_options & latex_name(std::string const & tn);
// the following lines have been generated for max. 14 parameters
    function_options & eval_func(eval_funcp_1 e);
    function_options & eval_func(eval_funcp_2 e);
    function_options & eval_func(eval_funcp_3 e);
    function_options & eval_func(eval_funcp_6 e);

    function_options & evalf_func(evalf_funcp_1 ef);
    function_options & evalf_func(evalf_funcp_2 ef);
    function_options & evalf_func(evalf_funcp_3 ef);
    function_options & evalf_func(evalf_funcp_6 ef);

    function_options & conjugate_func(conjugate_funcp_1 d);
    function_options & conjugate_func(conjugate_funcp_2 d);
    function_options & conjugate_func(conjugate_funcp_3 d);

    function_options & real_part_func(real_part_funcp_1 d);
    function_options & real_part_func(real_part_funcp_2 d);
    function_options & real_part_func(real_part_funcp_3 d);

    function_options & imag_part_func(imag_part_funcp_1 d);
    function_options & imag_part_func(imag_part_funcp_2 d);
    function_options & imag_part_func(imag_part_funcp_3 d);

    function_options & derivative_func(derivative_funcp_1 d);
    function_options & derivative_func(derivative_funcp_2 d);
    function_options & derivative_func(derivative_funcp_3 d);
    function_options & derivative_func(derivative_funcp_6 d);

    function_options & expl_derivative_func(expl_derivative_funcp_1 d);
    function_options & expl_derivative_func(expl_derivative_funcp_2 d);
    function_options & expl_derivative_func(expl_derivative_funcp_3 d);

    function_options & power_func(power_funcp_1 d);
    function_options & power_func(power_funcp_2 d);
    function_options & power_func(power_funcp_3 d);

    function_options & series_func(series_funcp_1 s);
    function_options & series_func(series_funcp_2 s);
    function_options & series_func(series_funcp_3 s);

    function_options & subs_func(subs_funcp_1 s);
    function_options & subs_func(subs_funcp_2 s);
    function_options & subs_func(subs_funcp_3 s);

    template <class Ctx> function_options & print_func(print_funcp_1 p)
    {
    	test_and_set_nparams(1);
    	set_print_func(Ctx::get_class_info_static().options.get_id(), print_funcp(p));
    	return *this;
    }
    template <class Ctx> function_options & print_func(print_funcp_2 p)
    {
    	test_and_set_nparams(2);
    	set_print_func(Ctx::get_class_info_static().options.get_id(), print_funcp(p));
    	return *this;
    }
    template <class Ctx> function_options & print_func(print_funcp_3 p)
    {
    	test_and_set_nparams(3);
    	set_print_func(Ctx::get_class_info_static().options.get_id(), print_funcp(p));
    	return *this;
    }
// end of generated lines
	function_options & eval_func(eval_funcp_exvector e);
	function_options & evalf_func(evalf_funcp_exvector ef);
	function_options & conjugate_func(conjugate_funcp_exvector d);
	function_options & real_part_func(real_part_funcp_exvector d);
	function_options & imag_part_func(imag_part_funcp_exvector d);
	function_options & derivative_func(derivative_funcp_exvector d);
	function_options & power_func(power_funcp_exvector d);
	function_options & series_func(series_funcp_exvector s);
	function_options & subs_func(series_funcp_exvector s);
	function_options & derivative_func(derivative_funcp_exvector_symbol d);

	template <class Ctx> function_options & print_func(print_funcp_exvector p)
	{
		print_use_exvector_args = true;
		set_print_func(Ctx::get_class_info_static().options.get_id(), print_funcp(p));
		return *this;
	}

	// python function calls
	function_options & eval_func(PyObject* e);
	function_options & evalf_func(PyObject* e);
	function_options & conjugate_func(PyObject* e);
	function_options & real_part_func(PyObject* e);
	function_options & imag_part_func(PyObject* e);
	function_options & derivative_func(PyObject* e);
	function_options & power_func(PyObject* e);
	function_options & series_func(PyObject* e);
	function_options & subs_func(PyObject* e);

	function_options & set_return_type(unsigned rt, tinfo_t rtt=nullptr);
	function_options & do_not_evalf_params();
	function_options & do_not_apply_chain_rule();
	function_options & remember(unsigned size, unsigned assoc_size=0,
	                            unsigned strategy=remember_strategies::delete_never);
	function_options & overloaded(unsigned o);

	std::string get_name() const { return name; }
	unsigned get_nparams() const { return nparams; }

	void set_python_func() { python_func = true; }

	void set_print_latex_func(PyObject* f);
	void set_print_dflt_func(PyObject* f);

	enum {
		eval_python_f 		= 0x0001,
		evalf_python_f 		= 0x0002,
		conjugate_python_f 	= 0x0004,
		real_part_python_f	= 0x0008,
		imag_part_python_f	= 0x0010,
		derivative_python_f	= 0x0020,
		power_python_f		= 0x0040,
		series_python_f		= 0x0080,
		subs_python_f           = 0x0100,
	};

public:
	bool has_derivative() const { return derivative_f != nullptr; }
	bool has_power() const { return power_f != nullptr; }
	void test_and_set_nparams(unsigned n);
	void set_print_func(unsigned id, print_funcp f);

	std::string name;
	std::string TeX_name;

	unsigned nparams;

	eval_funcp pynac_eval_f;
	eval_funcp eval_f;
	evalf_funcp evalf_f;
	conjugate_funcp conjugate_f;
	real_part_funcp real_part_f;
	imag_part_funcp imag_part_f;
	derivative_funcp derivative_f;
        expl_derivative_funcp expl_derivative_f;
	power_funcp power_f;
	series_funcp series_f;
        subs_funcp subs_f;
	std::vector<print_funcp> print_dispatch_table;

	bool evalf_params_first;
	bool apply_chain_rule;

	bool use_return_type;
	unsigned return_type;
	tinfo_t return_type_tinfo;

	bool use_remember;
	unsigned remember_size;
	unsigned remember_assoc_size;
	unsigned remember_strategy;

	bool eval_use_exvector_args;
	bool evalf_use_exvector_args;
	bool conjugate_use_exvector_args;
	bool real_part_use_exvector_args;
	bool imag_part_use_exvector_args;
	bool derivative_use_exvector_args;
        bool expl_derivative_use_exvector_args;
	bool power_use_exvector_args;
	bool series_use_exvector_args;
	bool print_use_exvector_args;

	unsigned python_func;

	unsigned functions_with_same_name;

	ex symtree;
};


/** Exception class thrown by classes which provide their own series expansion
 *  to signal that ordinary Taylor expansion is safe. */
class do_taylor {};


/** The class function is used to implement builtin functions like sin, cos...
	and user defined functions */
class function : public exprseq
{
	GINAC_DECLARE_REGISTERED_CLASS(function, exprseq)

	friend class print_order;
	// CINT has a linking problem
#ifndef __MAKECINT__
	friend void ginsh_get_ginac_functions();
#endif // def __MAKECINT__

	friend class remember_table_entry;
	// friend class remember_table_list;
	// friend class remember_table;

// member functions

	// other constructors
public:
	function(unsigned ser);
	// the following lines have been generated for some parameters
    function(unsigned ser, const ex & param1);
    function(unsigned ser, const ex & param1, const ex & param2);
    function(unsigned ser, const ex & param1, const ex & param2, const ex & param3);
    function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4);
    function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5);
    function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6);
	// end of generated lines
	function(unsigned ser, exprseq  es);
	function(unsigned ser, const exvector & v, bool discardable = false);
	function(unsigned ser, std::unique_ptr<exvector> vp);

	// functions overriding virtual functions from base classes
public:
	void print(const print_context & c, unsigned level = 0) const override;
	unsigned precedence() const override {return 70;}
	ex expand(unsigned options=0) const override;
	ex eval(int level=0) const override;
	ex evalf(int level=0, PyObject* parent=nullptr) const override;
	long calchash() const override;
	ex series(const relational & r, int order, unsigned options = 0) const override;
        void useries(flint_series_t& fp, int order) const override;
        bool match(const ex& pattern, exmap& map) const override;
        ex subs(const exmap & m, unsigned options = 0) const override;
	ex normal(exmap & repl, exmap & rev_lookup, int level = 0, unsigned options = 0) const override;
	ex thiscontainer(const exvector & v) const override;
	ex thiscontainer(std::unique_ptr<exvector> vp) const override;
	ex conjugate() const override;
	ex real_part() const override;
	ex imag_part() const override;
	bool info(unsigned inf) const override;
	//int compare(const basic &other) const;
        static ex unarchive(const archive_node &n, lst &sym_lst);
protected:
	ex derivative(const symbol & s) const override;
	bool is_equal_same_type(const basic & other) const override;
	bool match_same_type(const basic & other) const override;
	unsigned return_type() const override;
	tinfo_t return_type_tinfo() const override;
	
	// new virtual functions which can be overridden by derived classes
	// none
	
	// non-virtual functions in this class
protected:
	ex pderivative(unsigned diff_param) const; // partial differentiation
        ex expl_derivative(const symbol & s) const; // partial differentiation
	bool lookup_remember_table(ex & result) const;
	void store_remember_table(ex const & result) const;
public:
	static std::vector<function_options> & registered_functions();
	ex power(const ex & exp) const;
	static unsigned register_new(function_options const & opt);
	static unsigned current_serial;
	static unsigned find_function(const std::string &name, unsigned nparams);
	unsigned get_serial() const {return serial;}
	std::string get_name() const;
	unsigned get_domain() const { return domain; }
	void set_domain(unsigned d);
	void set_info(unsigned flag, bool value=true) { iflags.set(flag, value); }

// member variables

protected:
	unsigned serial;
	unsigned domain;
	infoflagbase iflags;
};

// utility functions/macros

template <typename T>
inline bool is_the_function(const ex & x)
{
	return is_exactly_a<function>(x)
	    && ex_to<function>(x).get_serial() == T::serial;
}

// Check whether OBJ is the specified symbolic function.
#define is_ex_the_function(OBJ, FUNCNAME) (GiNaC::is_the_function<FUNCNAME##_SERIAL>(OBJ))

bool has_function(const ex & x);
bool has_function(const ex& x, const std::string& s);
bool has_function(const std::vector<std::string>& v);
bool has_function(const ex& x,
                const std::vector<std::string>& v,
                bool all);
bool has_symbol_or_function(const ex & x);

} // namespace GiNaC

#endif // ndef __GINAC_FUNCTION_H__

