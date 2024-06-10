/** @file ginac.h
 *
 *  This include file includes all other public GiNaC headers. */

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

#ifndef __GINAC_H__
#define __GINAC_H__

#include "basic.h"

#include "ex.h"
#include "normal.h"
#include "archive.h"
#include "print.h"

#include "constant.h"
#include "infinity.h"
#include "lst.h"
#include "matrix.h"
#include "numeric.h"
#include "power.h"
#include "relational.h"
#include "symbol.h"
#include "pseries.h"
#include "wildcard.h"

#include "expair.h"
#include "expairseq.h"
#include "add.h"
#include "mul.h"
#include "upoly.h"
#include "mpoly.h"

#include "exprseq.h"
#include "function.h"

#include "inifcns.h"
#include "fderivative.h"
#include "operators.h"

#include "assume.h"
#include "sum.h"

#ifdef __MAKECINT__
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#endif

#endif // ndef __GINAC_H__
