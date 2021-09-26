/** @file utils.cpp
 *
 *  Implementation of several small and furry utilities needed within GiNaC
 *  but not of any interest to the user of the library. */

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

#include "ex.h"
#include "numeric.h"
#include "utils.h"

namespace GiNaC {


/** ctor for pole_error exception class. */
pole_error::pole_error(const std::string& what_arg, int a_degree)
	: domain_error(what_arg), deg(a_degree) { }

/** Return the degree of the pole_error exception class. */
numeric pole_error::degree() const
{
	return deg;
}

/** Integer binary logarithm */
unsigned log2(unsigned n)
{
	unsigned k;
	for (k = 0; n > 1; n >>= 1)
		++k;
	return k;
}


//////////
// flyweight chest of numbers is initialized here:
//////////

/** How many static objects were created?  Only the first one must create
 *  the static flyweights on the heap. */
int library_init::count = 0;

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wuninitialized"
// static numeric -120
const numeric *_num_120_p;
const ex _ex_120 = _ex_120;

// static numeric -60
const numeric *_num_60_p;
const ex _ex_60 = _ex_60;

// static numeric -48
const numeric *_num_48_p;
const ex _ex_48 = _ex_48;

// static numeric -30
const numeric *_num_30_p;
const ex _ex_30 = _ex_30;

// static numeric -25
const numeric *_num_25_p;
const ex _ex_25 = _ex_25;

// static numeric -24
const numeric *_num_24_p;
const ex _ex_24 = _ex_24;

// static numeric -20
const numeric *_num_20_p;
const ex _ex_20 = _ex_20;

// static numeric -18
const numeric *_num_18_p;
const ex _ex_18 = _ex_18;

// static numeric -15
const numeric *_num_15_p;
const ex _ex_15 = _ex_15;

// static numeric -12
const numeric *_num_12_p;
const ex _ex_12 = _ex_12;

// static numeric -11
const numeric *_num_11_p;
const ex _ex_11 = _ex_11;

// static numeric -10
const numeric *_num_10_p;
const ex _ex_10 = _ex_10;

// static numeric -9
const numeric *_num_9_p;
const ex _ex_9 = _ex_9;

// static numeric -8
const numeric *_num_8_p;
const ex _ex_8 = _ex_8;

// static numeric -7
const numeric *_num_7_p;
const ex _ex_7 = _ex_7;

// static numeric -6
const numeric *_num_6_p;
const ex _ex_6 = _ex_6;

// static numeric -5
const numeric *_num_5_p;
const ex _ex_5 = _ex_5;

// static numeric -4
const numeric *_num_4_p;
const ex _ex_4 = _ex_4;

// static numeric -3
const numeric *_num_3_p;
const ex _ex_3 = _ex_3;

// static numeric -2
const numeric *_num_2_p;
const ex _ex_2 = _ex_2;

// static numeric -1
const numeric *_num_1_p;
const ex _ex_1 = _ex_1;

// static numeric -1/2
const numeric *_num_1_2_p;
const ex _ex_1_2= _ex_1_2;

// static numeric -1/3
const numeric *_num_1_3_p;
const ex _ex_1_3= _ex_1_3;

// static numeric -1/4
const numeric *_num_1_4_p;
const ex _ex_1_4= _ex_1_4;

// static numeric 0
const numeric *_num0_p;
const basic *_num0_bp;
const ex _ex0 = _ex0;

// static numeric 1/4
const numeric *_num1_4_p;
const ex _ex1_4 = _ex1_4;

// static numeric 1/3
const numeric *_num1_3_p;
const ex _ex1_3 = _ex1_3;

// static numeric 1/2
const numeric *_num1_2_p;
const ex _ex1_2 = _ex1_2;

// static numeric 1
const numeric *_num1_p;
const ex _ex1 = _ex1;

// static numeric 2
const numeric *_num2_p;
const ex _ex2 = _ex2;

// static numeric 3
const numeric *_num3_p;
const ex _ex3 = _ex3;

// static numeric 4
const numeric *_num4_p;
const ex _ex4 = _ex4;

// static numeric 5
const numeric *_num5_p;
const ex _ex5 = _ex5;

// static numeric 6
const numeric *_num6_p;
const ex _ex6 = _ex6;

// static numeric 7
const numeric *_num7_p;
const ex _ex7 = _ex7;

// static numeric 8
const numeric *_num8_p;
const ex _ex8 = _ex8;

// static numeric 9
const numeric *_num9_p;
const ex _ex9 = _ex9;

// static numeric 10
const numeric *_num10_p;
const ex _ex10 = _ex10;

// static numeric 11
const numeric *_num11_p;
const ex _ex11 = _ex11;

// static numeric 12
const numeric *_num12_p;
const ex _ex12 = _ex12;

// static numeric 14
const numeric *_num14_p;
const ex _ex14 = _ex14;

// static numeric 15
const numeric *_num15_p;
const ex _ex15 = _ex15;

// static numeric 16
const numeric *_num16_p;
const ex _ex16 = _ex16;

// static numeric 18
const numeric *_num18_p;
const ex _ex18 = _ex18;

// static numeric 20
const numeric *_num20_p;
const ex _ex20 = _ex20;

// static numeric 21
const numeric *_num21_p;
const ex _ex21 = _ex21;

// static numeric 22
const numeric *_num22_p;
const ex _ex22 = _ex22;

// static numeric 24
const numeric *_num24_p;
const ex _ex24 = _ex24;

// static numeric 25
const numeric *_num25_p;
const ex _ex25 = _ex25;

// static numeric 26
const numeric *_num26_p;
const ex _ex26 = _ex26;

// static numeric 27
const numeric *_num27_p;
const ex _ex27 = _ex27;

// static numeric 28
const numeric *_num28_p;
const ex _ex28 = _ex28;

// static numeric 30
const numeric *_num30_p;
const ex _ex30 = _ex30;

// static numeric 36
const numeric *_num36_p;
const ex _ex36 = _ex36;

// static numeric 48
const numeric *_num48_p;
const ex _ex48 = _ex48;

// static numeric 60
const numeric *_num60_p;
const ex _ex60 = _ex60;

// static numeric 72
const numeric *_num72_p;
const ex _ex72 = _ex72;

// static numeric 120
const numeric *_num120_p;
const ex _ex120 = _ex120;

// static numeric 144
const numeric *_num144_p;
const ex _ex144 = _ex144;
#pragma clang diagnostic pop

/** Ctor of static initialization helpers.  The fist call to this is going
 *  to initialize the library, the others do nothing. */
library_init::library_init()
{
	if (count++==0) {
		(_num_120_p = new numeric(-120))->setflag(status_flags::dynallocated);
		(_num_60_p = new numeric(-60))->setflag(status_flags::dynallocated);
		(_num_48_p = new numeric(-48))->setflag(status_flags::dynallocated);
		(_num_30_p = new numeric(-30))->setflag(status_flags::dynallocated);
		(_num_25_p = new numeric(-25))->setflag(status_flags::dynallocated);
		(_num_24_p = new numeric(-24))->setflag(status_flags::dynallocated);
		(_num_20_p = new numeric(-20))->setflag(status_flags::dynallocated);
		(_num_18_p = new numeric(-18))->setflag(status_flags::dynallocated);
		(_num_15_p = new numeric(-15))->setflag(status_flags::dynallocated);
		(_num_12_p = new numeric(-12))->setflag(status_flags::dynallocated);
		(_num_11_p = new numeric(-11))->setflag(status_flags::dynallocated);
		(_num_10_p = new numeric(-10))->setflag(status_flags::dynallocated);
		(_num_9_p = new numeric(-9))->setflag(status_flags::dynallocated);
		(_num_8_p = new numeric(-8))->setflag(status_flags::dynallocated);
		(_num_7_p = new numeric(-7))->setflag(status_flags::dynallocated);
		(_num_6_p = new numeric(-6))->setflag(status_flags::dynallocated);
		(_num_5_p = new numeric(-5))->setflag(status_flags::dynallocated);
		(_num_4_p = new numeric(-4))->setflag(status_flags::dynallocated);
		(_num_3_p = new numeric(-3))->setflag(status_flags::dynallocated);
		(_num_2_p = new numeric(-2))->setflag(status_flags::dynallocated);
		(_num_1_p = new numeric(-1))->setflag(status_flags::dynallocated);
		(_num_1_2_p = new numeric(-1,2))->setflag(status_flags::dynallocated);
		(_num_1_3_p = new numeric(-1,3))->setflag(status_flags::dynallocated);
		(_num_1_4_p = new numeric(-1,4))->setflag(status_flags::dynallocated);
		(_num0_p = new numeric(0))->setflag(status_flags::dynallocated);
		_num0_bp  = _num0_p;  // Cf. class ex default ctor.
		(_num1_4_p = new numeric(1,4))->setflag(status_flags::dynallocated);
		(_num1_3_p = new numeric(1,3))->setflag(status_flags::dynallocated);
		(_num1_2_p = new numeric(1,2))->setflag(status_flags::dynallocated);
		(_num1_p = new numeric(1))->setflag(status_flags::dynallocated);
		(_num2_p = new numeric(2))->setflag(status_flags::dynallocated);
		(_num3_p = new numeric(3))->setflag(status_flags::dynallocated);
		(_num4_p = new numeric(4))->setflag(status_flags::dynallocated);
		(_num5_p = new numeric(5))->setflag(status_flags::dynallocated);
		(_num6_p = new numeric(6))->setflag(status_flags::dynallocated);
		(_num7_p = new numeric(7))->setflag(status_flags::dynallocated);
		(_num8_p = new numeric(8))->setflag(status_flags::dynallocated);
		(_num9_p = new numeric(9))->setflag(status_flags::dynallocated);
		(_num10_p = new numeric(10))->setflag(status_flags::dynallocated);
		(_num11_p = new numeric(11))->setflag(status_flags::dynallocated);
		(_num12_p = new numeric(12))->setflag(status_flags::dynallocated);
		(_num14_p = new numeric(14))->setflag(status_flags::dynallocated);
		(_num15_p = new numeric(15))->setflag(status_flags::dynallocated);
		(_num16_p = new numeric(16))->setflag(status_flags::dynallocated);
		(_num18_p = new numeric(18))->setflag(status_flags::dynallocated);
		(_num20_p = new numeric(20))->setflag(status_flags::dynallocated);
		(_num21_p = new numeric(21))->setflag(status_flags::dynallocated);
		(_num22_p = new numeric(22))->setflag(status_flags::dynallocated);
		(_num24_p = new numeric(24))->setflag(status_flags::dynallocated);
		(_num25_p = new numeric(25))->setflag(status_flags::dynallocated);
		(_num26_p = new numeric(26))->setflag(status_flags::dynallocated);
		(_num27_p = new numeric(27))->setflag(status_flags::dynallocated);
		(_num28_p = new numeric(28))->setflag(status_flags::dynallocated);
		(_num30_p = new numeric(30))->setflag(status_flags::dynallocated);
		(_num36_p = new numeric(36))->setflag(status_flags::dynallocated);
		(_num48_p = new numeric(48))->setflag(status_flags::dynallocated);
		(_num60_p = new numeric(60))->setflag(status_flags::dynallocated);
		(_num72_p = new numeric(72))->setflag(status_flags::dynallocated);
		(_num120_p = new numeric(120))->setflag(status_flags::dynallocated);
		(_num144_p = new numeric(144))->setflag(status_flags::dynallocated);

		new((void*)&_ex_120) ex(*_num_120_p);
		new((void*)&_ex_60) ex(*_num_60_p);
		new((void*)&_ex_48) ex(*_num_48_p);
		new((void*)&_ex_30) ex(*_num_30_p);
		new((void*)&_ex_25) ex(*_num_25_p);
		new((void*)&_ex_24) ex(*_num_24_p);
		new((void*)&_ex_20) ex(*_num_20_p);
		new((void*)&_ex_18) ex(*_num_18_p);
		new((void*)&_ex_15) ex(*_num_15_p);
		new((void*)&_ex_12) ex(*_num_12_p);
		new((void*)&_ex_11) ex(*_num_11_p);
		new((void*)&_ex_10) ex(*_num_10_p);
		new((void*)&_ex_9) ex(*_num_9_p);
		new((void*)&_ex_8) ex(*_num_8_p);
		new((void*)&_ex_7) ex(*_num_7_p);
		new((void*)&_ex_6) ex(*_num_6_p);
		new((void*)&_ex_5) ex(*_num_5_p);
		new((void*)&_ex_4) ex(*_num_4_p);
		new((void*)&_ex_3) ex(*_num_3_p);
		new((void*)&_ex_2) ex(*_num_2_p);
		new((void*)&_ex_1) ex(*_num_1_p);
		new((void*)&_ex_1_2) ex(*_num_1_2_p);
		new((void*)&_ex_1_3) ex(*_num_1_3_p);
		new((void*)&_ex_1_4) ex(*_num_1_4_p);
		new((void*)&_ex0) ex(*_num0_p);
		new((void*)&_ex1_4) ex(*_num1_4_p);
		new((void*)&_ex1_3) ex(*_num1_3_p);
		new((void*)&_ex1_2) ex(*_num1_2_p);
		new((void*)&_ex1) ex(*_num1_p);
		new((void*)&_ex2) ex(*_num2_p);
		new((void*)&_ex3) ex(*_num3_p);
		new((void*)&_ex4) ex(*_num4_p);
		new((void*)&_ex5) ex(*_num5_p);
		new((void*)&_ex6) ex(*_num6_p);
		new((void*)&_ex7) ex(*_num7_p);
		new((void*)&_ex8) ex(*_num8_p);
		new((void*)&_ex9) ex(*_num9_p);
		new((void*)&_ex10) ex(*_num10_p);
		new((void*)&_ex11) ex(*_num11_p);
		new((void*)&_ex12) ex(*_num12_p);
		new((void*)&_ex14) ex(*_num14_p);
		new((void*)&_ex15) ex(*_num15_p);
		new((void*)&_ex16) ex(*_num16_p);
		new((void*)&_ex18) ex(*_num18_p);
		new((void*)&_ex20) ex(*_num20_p);
		new((void*)&_ex21) ex(*_num21_p);
		new((void*)&_ex22) ex(*_num22_p);
		new((void*)&_ex24) ex(*_num24_p);
		new((void*)&_ex25) ex(*_num25_p);
		new((void*)&_ex26) ex(*_num26_p);
		new((void*)&_ex27) ex(*_num27_p);
		new((void*)&_ex28) ex(*_num28_p);
		new((void*)&_ex30) ex(*_num30_p);
		new((void*)&_ex36) ex(*_num36_p);
		new((void*)&_ex48) ex(*_num48_p);
		new((void*)&_ex60) ex(*_num60_p);
		new((void*)&_ex72) ex(*_num72_p);
		new((void*)&_ex120) ex(*_num120_p);
		new((void*)&_ex144) ex(*_num144_p);

		// Initialize print context class info (this is not strictly necessary
		// but we do it anyway to make print_context_class_info::dump_hierarchy()
		// output the whole hierarchy whether or not the classes are actually
		// used)
		print_context::get_class_info_static();
		print_dflt::get_class_info_static();
		print_latex::get_class_info_static();
		print_python::get_class_info_static();
		print_python_repr::get_class_info_static();
		print_tree::get_class_info_static();
	}
}


/** Dtor of static initialization helpers.  The last call to this is going
 *  to shut down the library, the others do nothing. */
library_init::~library_init()
{
	if (--count==0) {
		// In theory, we would have to clean up here.  But since we were
		// only initializing memory in the ctor and that memory is reclaimed
		// anyways by the OS when the program exits, we skip this.
	}
}

// comment skeleton for header files


// member functions

	// default constructor, destructor, copy constructor and assignment operator
	// none

	// other constructors
	// none

	// functions overriding virtual functions from base classes
	// none
	
	// new virtual functions which can be overridden by derived classes
	// none

	// non-virtual functions in this class
	// none

// member variables
// none
	


// comment skeleton for implementation files


//////////
// default constructor, destructor, copy constructor and assignment operator
//////////

// public
// protected

//////////
// other constructors
//////////

// public
// none

//////////
// functions overriding virtual functions from base classes
//////////

// public
// protected
// none

//////////
// new virtual functions which can be overridden by derived classes
//////////

// public
// protected
// none

//////////
// non-virtual functions in this class
//////////

// public
// protected
// none

//////////
// static member variables
//////////

// protected
// private
// none


} // namespace GiNaC
