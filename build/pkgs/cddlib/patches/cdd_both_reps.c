/* cdd_both_reps.c: compute reduced H and V representation of polytope
   by Volker Braun <vbraun@stp.dias.ie>
   
   The input is taken from stdin and can be either a 
   H or V representation, not necessarily reduced.

   based on testcdd1.c, redcheck.c, and of course the cdd library
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Standard ftp site: ftp.ifor.math.ethz.ch, Directory: pub/fukuda/cdd
*/

/*  This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#include "setoper.h"
#include "cdd.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>





void compute_adjacency(dd_MatrixPtr Rep, dd_ErrorType* err_ptr)
{
  dd_SetFamilyPtr AdjacencyGraph;
  if (*err_ptr != dd_NoError) return;

  switch (Rep->representation) {
  case dd_Inequality: 
    printf("Facet graph\n");
    break;
  case dd_Generator: 
    printf("Vertex graph\n");
    break;
  case dd_Unspecified:
    printf("unknown representation type!\n");
  default:
    printf("This should be unreachable!\n");
    exit(2);
  }

  /* Output adjacency of vertices/rays/lines */
  if (Rep->rowsize > 0) {  /* workaround for bug with empty polyhedron */
    /* compute adjacent vertices/rays/lines */
    AdjacencyGraph = dd_Matrix2Adjacency(Rep, err_ptr);
    if (*err_ptr == dd_NoError) {
      dd_WriteSetFamily(stdout,AdjacencyGraph);
      dd_FreeSetFamily(AdjacencyGraph);
    }
  } else {
    printf("begin\n");
    printf("  0    0\n");
    printf("end\n");
  }

  printf("\n");
}


void minimal_Vrep_Hrep(dd_MatrixPtr M, 
		       dd_MatrixPtr* Vrep_ptr, dd_MatrixPtr* Hrep_ptr, 
		       dd_ErrorType* err_ptr)
{
  dd_PolyhedraPtr poly;
  dd_rowindex newpos;
  dd_rowset impl_linset,redset;
  dd_MatrixPtr Vrep, Hrep;

  if (*err_ptr != dd_NoError) return;

   /* compute the second representation */
  poly = dd_DDMatrix2Poly(M, err_ptr);
  if (*err_ptr != dd_NoError) return;

  if (*err_ptr == dd_NoError) {
    /* compute canonical H-representation */
    Hrep = dd_CopyInequalities(poly);
    if (Hrep->rowsize > 0) {  /* workaround for bug with empty matrix */
      dd_MatrixCanonicalize(&Hrep, &impl_linset, &redset, &newpos, err_ptr);
      if (*err_ptr == dd_NoError) {
	set_free(redset);
	set_free(impl_linset);
	free(newpos);
      }
    }
    if (*err_ptr == dd_NoError) (*Hrep_ptr) = Hrep;
  }

  if (*err_ptr == dd_NoError) {
    /* compute canonical V-representation */
    Vrep = dd_CopyGenerators(poly);
    if (Vrep->rowsize > 0) {  /* workaround for bug with empty matrix */
      dd_MatrixCanonicalize(&Vrep, &impl_linset, &redset, &newpos, err_ptr);
      if (*err_ptr == dd_NoError) {
	set_free(redset);
	set_free(impl_linset);
	free(newpos);
      }
    }
    if (*err_ptr == dd_NoError) (*Vrep_ptr) = Vrep;
  }

  dd_FreePolyhedra(poly);
}


void print_both_reps(dd_MatrixPtr Vrep, dd_MatrixPtr Hrep)
{
  /* Output V-representation */
  dd_WriteMatrix(stdout,Vrep);
  printf("\n");

  /* Output H-representation */
  dd_WriteMatrix(stdout,Hrep);
  printf("\n");
}


void compute_both_reps(dd_MatrixPtr M, dd_ErrorType* err_ptr)
{
  dd_MatrixPtr Vrep, Hrep;
  minimal_Vrep_Hrep(M, &Vrep, &Hrep, err_ptr);
  if (*err_ptr != dd_NoError) return;

  print_both_reps(Vrep, Hrep);
  dd_FreeMatrix(Hrep);
  dd_FreeMatrix(Vrep);
}


void compute_all(dd_MatrixPtr M, dd_ErrorType* err_ptr)
{
  dd_MatrixPtr Vrep, Hrep;
  minimal_Vrep_Hrep(M, &Vrep, &Hrep, err_ptr);
  if (*err_ptr != dd_NoError) return;

  print_both_reps(Vrep, Hrep);
  compute_adjacency(Vrep, err_ptr);
  compute_adjacency(Hrep, err_ptr);
  dd_FreeMatrix(Hrep);
  dd_FreeMatrix(Vrep);
}



void usage(char *name)
{
  printf("No known option specified, I don't know what to do!\n"
	 "Usage:\n"
	 "%s --option\n"
	 "where --option is precisely one of the following:\n\n"
	 "  --all: Compute everything.\n"
	 "    This will compute minimal H-,V-representation and vertex and facet graph.\n"
	 "\n"
	 "  --reps: Compute both a minimal H- and minimal V-representation.\n"
	 "\n"
	 "  --adjacency: Compute adjacency information only.\n"
	 "    The input is assumed to be a minimal representation, as, for example, computed\n"
	 "    by --reps. Warning, you will not get the correct answer if the input\n"
	 "    representation is not minimal! The output is the vertex or facet graph,\n"
	 "    depending on the input.\n"
	 "\n"
	 "The input data is a H- or V-representation in cdd's ine/ext format and\n"
	 "is in each case read from stdin.\n", 
	 name);
}


enum command_line_arguments { ALL, REPS, ADJACENCY };


int parse_arguments(char* arg, enum command_line_arguments* option)
{
  if (strcmp(arg,"--all")==0) {
    *option = ALL;
    return 0;
  }
  if (strcmp(arg,"--reps")==0) {
    *option = REPS;
    return 0;
  }
  if (strcmp(arg,"--adjacency")==0) {
    *option = ADJACENCY;
    return 0;
  }
  printf("Unknown option: %s\n", arg);
  return 1;
}


int main(int argc, char *argv[])
{
  dd_ErrorType err=dd_NoError;
  dd_MatrixPtr M;
  enum command_line_arguments option;

  if (argc!=2 || parse_arguments(argv[1],&option)) {
    usage(argv[0]);
    return 0;
  } 

  dd_set_global_constants(); 

  /* Read data from stdin */
  M = dd_PolyFile2Matrix(stdin, &err);
  if (err != dd_NoError) {
    printf("I was unable to parse the input data!\n");
    dd_WriteErrorMessages(stdout,err);
    dd_free_global_constants();
    return 1;
  }

  switch (option) {
  case ALL:
    compute_all(M,&err);
    break;
  case REPS:
    compute_both_reps(M,&err);
    break;
  case ADJACENCY:
    compute_adjacency(M,&err);
    break;
  default:
    printf("unreachable option %d\n", option);
    exit(3); /* unreachable */
  }
  
  /* cleanup */
  dd_FreeMatrix(M);
  if (err != dd_NoError) {
    dd_WriteErrorMessages(stdout,err);
  }

  dd_free_global_constants();
  return 0;
}



