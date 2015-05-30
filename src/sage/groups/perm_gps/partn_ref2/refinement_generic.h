/*******************************************************************************
#       Copyright (C) 2012 Thomas Feulner <thomas.feulner@uni-bayreuth.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#ifndef BACKTRACK_WITHLATEX_DEBUG
#define BACKTRACK_WITHLATEX_DEBUG 0 // set to 1 if you have to debug
#endif

static long *global_refine_vals_array = 0;

static int my_comp_func(const void *a, const void *b)
{
	long val_a = global_refine_vals_array[*((int *) a)];
	long val_b = global_refine_vals_array[*((int *) b)];

	if(val_a == val_b)
		return 0;
	if(val_a < val_b)
		return -1;
	return 1;
}

int in_array(int *list, int length, int el){
	int i=0;
	for(; i < length; i++ )
	{
		if ( el == list[i] )
		{
			return 1;
		}
	}
	return 0;
}

