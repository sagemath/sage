/*
Planarity-Related Graph Algorithms Project
Copyright (c) 1997-2010, John M. Boyer
All rights reserved. Includes a reference implementation of the following:

* John M. Boyer. "Simplified O(n) Algorithms for Planar Graph Embedding,
  Kuratowski Subgraph Isolation, and Related Problems". Ph.D. Dissertation,
  University of Victoria, 2001.

* John M. Boyer and Wendy J. Myrvold. "On the Cutting Edge: Simplified O(n)
  Planarity by Edge Addition". Journal of Graph Algorithms and Applications,
  Vol. 8, No. 3, pp. 241-273, 2004.

* John M. Boyer. "A New Method for Efficiently Generating Planar Graph
  Visibility Representations". In P. Eades and P. Healy, editors,
  Proceedings of the 13th International Conference on Graph Drawing 2005,
  Lecture Notes Comput. Sci., Volume 3843, pp. 508-511, Springer-Verlag, 2006.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* Neither the name of the Planarity-Related Graph Algorithms Project nor the names
  of its contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#define GRAPHSTRUCTURE_C

#include <stdlib.h>

#include "graphStructures.h"
#include "graph.h"

/* Imported functions for FUNCTION POINTERS */

extern int  _CreateFwdArcLists(graphP);
extern void _CreateDFSTreeEmbedding(graphP theGraph);
extern int  _SortVertices(graphP theGraph);
extern void _EmbedBackEdgeToDescendant(graphP theGraph, int RootSide, int RootVertex, int W, int WPrevLink);
extern void _WalkUp(graphP theGraph, int I, int J);
extern int  _WalkDown(graphP theGraph, int I, int RootVertex);
extern int  _MergeBicomps(graphP theGraph, int I, int RootVertex, int W, int WPrevLink);
extern int  _HandleBlockedDescendantBicomp(graphP theGraph, int I, int RootVertex, int R, int *pRout, int *pW, int *pWPrevLink);
extern int  _HandleInactiveVertex(graphP theGraph, int BicompRoot, int *pW, int *pWPrevLink);
extern int  _MarkDFSPath(graphP theGraph, int ancestor, int descendant);
extern int  _HandleBlockedEmbedIteration(graphP theGraph, int I);
extern int  _EmbedPostprocess(graphP theGraph, int I, int edgeEmbeddingResult);
extern int  _CheckEmbeddingIntegrity(graphP theGraph, graphP origGraph);
extern int  _CheckObstructionIntegrity(graphP theGraph, graphP origGraph);
extern int  _ReadPostprocess(graphP theGraph, void *extraData, long extraDataSize);
extern int  _WritePostprocess(graphP theGraph, void **pExtraData, long *pExtraDataSize);

/* Internal util functions for FUNCTION POINTERS */

int  _HideVertex(graphP theGraph, int vertex);
void _HideEdge(graphP theGraph, int arcPos);
void _RestoreEdge(graphP theGraph, int arcPos);
int  _ContractEdge(graphP theGraph, int e);
int  _IdentifyVertices(graphP theGraph, int u, int v, int eBefore);
int  _RestoreVertex(graphP theGraph);

/********************************************************************
 Private functions, except exported within library
 ********************************************************************/

void _ClearIsolatorContext(graphP theGraph);
void _FillVisitedFlags(graphP theGraph, int FillValue);
int  _FillVisitedFlagsInBicomp(graphP theGraph, int BicompRoot, int FillValue);
int  _FillVisitedFlagsInOtherBicomps(graphP theGraph, int BicompRoot, int FillValue);
void _FillVisitedFlagsInUnembeddedEdges(graphP theGraph, int FillValue);
int  _SetVertexTypeInBicomp(graphP theGraph, int BicompRoot, int theType);

int  _HideInternalEdges(graphP theGraph, int vertex);
int  _RestoreInternalEdges(graphP theGraph, int stackBottom);
int  _RestoreHiddenEdges(graphP theGraph, int stackBottom);

int  _GetBicompSize(graphP theGraph, int BicompRoot);
int  _DeleteUnmarkedEdgesInBicomp(graphP theGraph, int BicompRoot);
int  _ClearInvertedFlagsInBicomp(graphP theGraph, int BicompRoot);

void _InitFunctionTable(graphP theGraph);

/********************************************************************
 Private functions.
 ********************************************************************/

void _ClearGraph(graphP theGraph);

int  _GetRandomNumber(int NMin, int NMax);

/* Private functions for which there are FUNCTION POINTERS */

void _InitGraphNode(graphP theGraph, int I);
void _InitVertexRec(graphP theGraph, int I);

int  _InitGraph(graphP theGraph, int N);
void _ReinitializeGraph(graphP theGraph);
int  _EnsureArcCapacity(graphP theGraph, int requiredArcCapacity);

/********************************************************************
 gp_New()
 Constructor for graph object.
 Can create two graphs if restricted to no dynamic memory.
 ********************************************************************/

graphP gp_New()
{
graphP theGraph = (graphP) malloc(sizeof(baseGraphStructure));

     if (theGraph != NULL)
     {
         theGraph->G = NULL;
         theGraph->V = NULL;

         theGraph->BicompLists = NULL;
         theGraph->DFSChildLists = NULL;
         theGraph->theStack = NULL;

         theGraph->buckets = NULL;
         theGraph->bin = NULL;

         theGraph->extFace = NULL;

         theGraph->edgeHoles = NULL;

         theGraph->extensions = NULL;

         _InitFunctionTable(theGraph);

         _ClearGraph(theGraph);
     }

     return theGraph;
}

/********************************************************************
 _InitFunctionTable()

 If you add functions to the function table, then they must be
 initialized here, but you must also add the new function pointer
 to the definition of the graphFunctionTable in graphFunctionTable.h

 Function headers for the functions used to initialize the table are
 classified at the top of this file as either imported from other
 compilation units (extern) or private to this compilation unit.
 Search for FUNCTION POINTERS in this file to see where to add the
 function header.
 ********************************************************************/

void _InitFunctionTable(graphP theGraph)
{
     theGraph->functions.fpCreateFwdArcLists = _CreateFwdArcLists;
     theGraph->functions.fpCreateDFSTreeEmbedding = _CreateDFSTreeEmbedding;
     theGraph->functions.fpEmbedBackEdgeToDescendant = _EmbedBackEdgeToDescendant;
     theGraph->functions.fpWalkUp = _WalkUp;
     theGraph->functions.fpWalkDown = _WalkDown;
     theGraph->functions.fpMergeBicomps = _MergeBicomps;
     theGraph->functions.fpHandleBlockedDescendantBicomp = _HandleBlockedDescendantBicomp;
     theGraph->functions.fpHandleInactiveVertex = _HandleInactiveVertex;
     theGraph->functions.fpHandleBlockedEmbedIteration = _HandleBlockedEmbedIteration;
     theGraph->functions.fpEmbedPostprocess = _EmbedPostprocess;
     theGraph->functions.fpMarkDFSPath = _MarkDFSPath;
     theGraph->functions.fpCheckEmbeddingIntegrity = _CheckEmbeddingIntegrity;
     theGraph->functions.fpCheckObstructionIntegrity = _CheckObstructionIntegrity;

     theGraph->functions.fpInitGraphNode = _InitGraphNode;
     theGraph->functions.fpInitVertexRec = _InitVertexRec;

     theGraph->functions.fpInitGraph = _InitGraph;
     theGraph->functions.fpReinitializeGraph = _ReinitializeGraph;
     theGraph->functions.fpEnsureArcCapacity = _EnsureArcCapacity;
     theGraph->functions.fpSortVertices = _SortVertices;

     theGraph->functions.fpReadPostprocess = _ReadPostprocess;
     theGraph->functions.fpWritePostprocess = _WritePostprocess;

     theGraph->functions.fpHideVertex = _HideVertex;
     theGraph->functions.fpHideEdge = _HideEdge;
     theGraph->functions.fpRestoreEdge = _RestoreEdge;
     theGraph->functions.fpContractEdge = _ContractEdge;
     theGraph->functions.fpIdentifyVertices = _IdentifyVertices;
     theGraph->functions.fpRestoreVertex = _RestoreVertex;
}

/********************************************************************
 gp_InitGraph()
 Allocates memory for vertex and edge records now that N is known.
 The arcCapacity is set to (2 * DEFAULT_EDGE_LIMIT * N) unless it
	 has already been set by gp_EnsureArcCapacity()

 For G, we need N vertex nodes, N more vertex nodes for root copies,
         and arcCapacity edge records.

 For V, we need N vertex records.

 The BicompLists and DFSChildLists are of size N and start out empty.

 The stack, initially empty, is made big enough for a pair of integers
	 per edge record, or 2 * arcCapacity.

 The edgeHoles stack, initially empty, is set to arcCapacity / 2,
	 which is big enough to push every edge (to indicate an edge
	 you only need to indicate one of its two edge records)

 buckets and bin are both O(N) in size.  They are used by
        CreateSortedSeparatedDFSChildLists()

  Returns OK on success, NOTOK on all failures.
          On NOTOK, graph extensions are freed so that the graph is
          returned to the post-condition of gp_New().
 ********************************************************************/

int gp_InitGraph(graphP theGraph, int N)
{
	// valid params check
	if (theGraph == NULL || N <= 0)
        return NOTOK;

	// Should not call init a second time; use reinit
	if (theGraph->N > 0)
		return NOTOK;

    return theGraph->functions.fpInitGraph(theGraph, N);
}

int  _InitGraph(graphP theGraph, int N)
{
int I, edgeOffset, arcCapacity, Gsize, Vsize, stackSize;

/* Compute the vertex and edge capacities of the graph */

     Vsize = 2*N;
     edgeOffset = Vsize;
     arcCapacity = theGraph->arcCapacity > 0 ? theGraph->arcCapacity : 2*DEFAULT_EDGE_LIMIT*N;
     Gsize = edgeOffset + arcCapacity;
     stackSize = 2 * arcCapacity;
     stackSize = stackSize < 6*N ? 6*N : stackSize;

/* Allocate memory as described above */

     if ((theGraph->G = (graphNodeP) malloc(Gsize*sizeof(graphNode))) == NULL ||
         (theGraph->V = (vertexRecP) malloc(N*sizeof(vertexRec))) == NULL ||
         (theGraph->BicompLists = LCNew(N)) == NULL ||
         (theGraph->DFSChildLists = LCNew(N)) == NULL ||
         (theGraph->theStack = sp_New(stackSize)) == NULL ||
         (theGraph->buckets = (int *) malloc(N * sizeof(int))) == NULL ||
         (theGraph->bin = LCNew(N)) == NULL ||
         (theGraph->extFace = (extFaceLinkRecP) malloc(Vsize*sizeof(extFaceLinkRec))) == NULL ||
         (theGraph->edgeHoles = sp_New(arcCapacity / 2)) == NULL ||
         0)
     {
         _ClearGraph(theGraph);
         return NOTOK;
     }

/* Initialize memory */

     theGraph->N = N;
     theGraph->edgeOffset = edgeOffset;
     theGraph->arcCapacity = Gsize - edgeOffset;

     for (I = 0; I < Gsize; I++)
          theGraph->functions.fpInitGraphNode(theGraph, I);

     for (I = 0; I < N; I++)
          theGraph->functions.fpInitVertexRec(theGraph, I);

     for (I = 0; I < Vsize; I++)
     {
         theGraph->extFace[I].vertex[0] = theGraph->extFace[I].vertex[1] = NIL;
         theGraph->extFace[I].inversionFlag = 0;
     }

     _ClearIsolatorContext(theGraph);

     return OK;
}

/********************************************************************
 gp_ReinitializeGraph()
 Reinitializes a graph, restoring it to the state it was in immediately
 after gp_InitGraph() processed it.
 ********************************************************************/

void gp_ReinitializeGraph(graphP theGraph)
{
	if (theGraph == NULL || theGraph->N <= 0)
		return;

    theGraph->functions.fpReinitializeGraph(theGraph);
}

void _ReinitializeGraph(graphP theGraph)
{
int  I, N = theGraph->N, edgeOffset = theGraph->edgeOffset;
int  Vsize = 2*N, Gsize = edgeOffset + theGraph->arcCapacity;

     theGraph->M = 0;
     theGraph->internalFlags = theGraph->embedFlags = 0;

     for (I = 0; I < Gsize; I++)
          theGraph->functions.fpInitGraphNode(theGraph, I);

     for (I = 0; I < N; I++)
          theGraph->functions.fpInitVertexRec(theGraph, I);

     for (I = 0; I < Vsize; I++)
     {
         theGraph->extFace[I].vertex[0] = theGraph->extFace[I].vertex[1] = NIL;
         theGraph->extFace[I].inversionFlag = 0;
     }

     _ClearIsolatorContext(theGraph);

     LCReset(theGraph->BicompLists);
     LCReset(theGraph->DFSChildLists);
     sp_ClearStack(theGraph->theStack);
     LCReset(theGraph->bin);
     sp_ClearStack(theGraph->edgeHoles);
}

/********************************************************************
 gp_GetArcCapacity()
 Returns the arcCapacity of theGraph, which is twice the maximum
 number of edges that can be added to the theGraph.
 ********************************************************************/
int gp_GetArcCapacity(graphP theGraph)
{
	return theGraph->arcCapacity;
}

/********************************************************************
 gp_EnsureArcCapacity()
 This method ensures that theGraph is or will be capable of storing
 at least requiredArcCapacity edge records.  Two edge records are
 needed per edge.

 This method is most performant when invoked immediately after
 gp_New(), since it must only set the arcCapacity and then let
 normal initialization to occur through gp_InitGraph().

 This method is also a constant time operation if the graph already
 has at least the requiredArcCapacity, since it will return OK
 without making any structural changes.

 This method is generally more performant if it is invoked before
 attaching extensions to the graph.  Some extensions associate
 parallel data with edge records, which is a faster operation if
 the associated data is created and initialized only after the
 proper arcCapacity is specified.

 If the graph has been initialized and has a lower arc capacity,
 then the array of edge records is reallocated to satisfy the
 requiredArcCapacity.  The new array contains the old edges and
 edge holes at the same locations, and all newly created edge records
 are initialized.

 Also, if the arc capacity must be increased, then the
 arcCapacity member of theGraph is changed and both
 theStack and edgeHoles are expanded (since the sizes of both
 are based on the arc capacity).

 Extensions that add to data associated with edges must overload
 this method to ensure capacity in the parallel extension data
 structures.  An extension can return NOTOK if it does not
 support arc capacity expansion.  The extension function will
 not be called if arcCapacity is expanded before the graph is
 initialized, and it is assumed that extensions will allocate
 parallel data structures according to the arc capacity.

 If an extension supports arc capacity expansion, then higher
 performance can be obtained by using the method of unhooking
 the initializers for individual edge records before invoking
 the superclass version of fpEnsureArcCapacity().  Ideally,
 application authors should ensure the proper arc capacity before
 attaching extensions to achieve better performance.

 Returns NOTOK on failure to reallocate the edge record array to
         satisfy the requiredArcCapacity, or if the requested
         capacity is odd
         OK if reallocation is not required or if reallocation succeeds
 ********************************************************************/
int gp_EnsureArcCapacity(graphP theGraph, int requiredArcCapacity)
{
	if (theGraph == NULL || requiredArcCapacity <= 0)
		return NOTOK;

	// Train callers to only ask for an even number of arcs, since
	// two are required per edge or directed edge.
	if (requiredArcCapacity & 1)
		return NOTOK;

    if (theGraph->arcCapacity >= requiredArcCapacity)
    	return OK;

    // In the special case where gp_InitGraph() has not yet been called,
    // we can simply set the higher arcCapacity since normal initialization
    // will then allocate the correct number of edge records.
    if (theGraph->N == 0)
    {
    	theGraph->arcCapacity = requiredArcCapacity;
    	return OK;
    }

    // Try to expand the arc capacity
    return theGraph->functions.fpEnsureArcCapacity(theGraph, requiredArcCapacity);
}

int _EnsureArcCapacity(graphP theGraph, int requiredArcCapacity)
{
stackP newStack;
int J, Gsize=theGraph->edgeOffset + theGraph->arcCapacity;
int newGsize = theGraph->edgeOffset + requiredArcCapacity;

    // Expand theStack
    if (sp_GetCapacity(theGraph->theStack) < 2 * requiredArcCapacity)
    {
    	int stackSize = 2 * requiredArcCapacity;

    	if (stackSize < 6*theGraph->N)
    	{
			// NOTE: Since this routine only makes the stack bigger, this
    		//       calculation is not needed here because we already ensured
			//       we had stack capacity of the greater of 2*arcs and 6*N
    		//       But we do it for clarity and consistency (e.g. so this rule
    		//       is not forgotten whenever a "SetArcCapacity" method is added)
    		stackSize = 6*theGraph->N;
    	}

    	if ((newStack = sp_New(stackSize)) == NULL)
    		return NOTOK;

    	sp_CopyContent(newStack, theGraph->theStack);
    	sp_Free(&theGraph->theStack);
    	theGraph->theStack = newStack;
    }

	// Expand edgeHoles
    if ((newStack = sp_New(requiredArcCapacity / 2)) == NULL)
    	return NOTOK;

	sp_CopyContent(newStack, theGraph->edgeHoles);
    sp_Free(&theGraph->edgeHoles);
    theGraph->edgeHoles = newStack;

	// Reallocate the GraphNode array to the new size,
    theGraph->G = (graphNodeP) realloc(theGraph->G, newGsize*sizeof(graphNode));
    if (theGraph->G == NULL)
    	return NOTOK;

    // Initialize the new edge records
    for (J = Gsize; J < newGsize; J++)
         theGraph->functions.fpInitGraphNode(theGraph, J);

    // The new arcCapacity has been successfully achieved
	theGraph->arcCapacity = requiredArcCapacity;
	return OK;
}

/********************************************************************
 _InitGraphNode()
 Sets the fields in a single graph node structure to initial values
 ********************************************************************/

void _InitGraphNode(graphP theGraph, int J)
{
     theGraph->G[J].v = NIL;
     gp_SetPrevArc(theGraph, J, NIL);
     gp_SetNextArc(theGraph, J, NIL);
     theGraph->G[J].visited = 0;
     theGraph->G[J].type = TYPE_UNKNOWN;
     theGraph->G[J].flags = 0;
}

/********************************************************************
 _InitVertexRec()
 Sets the fields in a single vertex record to initial values
 ********************************************************************/

void _InitVertexRec(graphP theGraph, int I)
{
    gp_SetFirstArc(theGraph, I, gp_AdjacencyListEndMark(I));
    gp_SetLastArc(theGraph, I, gp_AdjacencyListEndMark(I));
    theGraph->V[I].leastAncestor =
    theGraph->V[I].Lowpoint = I;
    theGraph->V[I].DFSParent = NIL;
    theGraph->V[I].adjacentTo = NIL;
    theGraph->V[I].pertinentBicompList = NIL;
    theGraph->V[I].separatedDFSChildList = NIL;
    theGraph->V[I].fwdArcList = NIL;
}

/********************************************************************
 _ClearIsolatorContext()
 ********************************************************************/

void _ClearIsolatorContext(graphP theGraph)
{
isolatorContextP IC = &theGraph->IC;

     IC->minorType = 0;
     IC->v = IC->r = IC->x = IC->y = IC->w = IC->px = IC->py = IC->z =
     IC->ux = IC->dx = IC->uy = IC->dy = IC->dw = IC->uz = IC->dz = NIL;
}

/********************************************************************
 _FillVisitedFlags()
 ********************************************************************/

void _FillVisitedFlags(graphP theGraph, int FillValue)
{
int  i;
int  limit = theGraph->edgeOffset + 2*(theGraph->M + sp_GetCurrentSize(theGraph->edgeHoles));

     for (i=0; i < limit; i++)
          theGraph->G[i].visited = FillValue;
}

/********************************************************************
 _FillVisitedFlagsInBicomp()

 Places the FillValue into the 'visited' members of the vertices and
 arcs in the bicomp rooted by BicompRoot.

 This method uses the stack but preserves whatever may have been
 on it.  In debug mode, it will return NOTOK if the stack overflows.
 This method pushes at most one integer per vertex in the bicomp.

 Returns OK on success, NOTOK on implementation failure.
 ********************************************************************/

int  _FillVisitedFlagsInBicomp(graphP theGraph, int BicompRoot, int FillValue)
{
int  V, J;
int  stackBottom = sp_GetCurrentSize(theGraph->theStack);

     sp_Push(theGraph->theStack, BicompRoot);
     while (sp_GetCurrentSize(theGraph->theStack) > stackBottom)
     {
          sp_Pop(theGraph->theStack, V);
          theGraph->G[V].visited = FillValue;

          J = gp_GetFirstArc(theGraph, V);
          while (gp_IsArc(theGraph, J))
          {
             theGraph->G[J].visited = FillValue;

             if (theGraph->G[J].type == EDGE_DFSCHILD)
                 sp_Push(theGraph->theStack, theGraph->G[J].v);

             J = gp_GetNextArc(theGraph, J);
          }
     }
     return OK;
}

/********************************************************************
 _FillVisitedFlagsInOtherBicomps()
 Typically, we want to clear or set all visited flags in the graph
 (see _FillVisitedFlags).  However, in some algorithms this would be
 too costly, so it is necessary to clear or set the visited flags only
 in one bicomp (see _FillVisitedFlagsInBicomp), then do some processing
 that sets some of the flags then performs some tests.  If the tests
 are positive, then we can clear or set all the visited flags in the
 other bicomps (the processing may have set the visited flags in the
 one bicomp in a particular way that we want to retain, so we skip
 the given bicomp).
 ********************************************************************/

int  _FillVisitedFlagsInOtherBicomps(graphP theGraph, int BicompRoot, int FillValue)
{
int  R, edgeOffset = theGraph->edgeOffset;

     for (R = theGraph->N; R < edgeOffset; R++)
     {
          if (R != BicompRoot && gp_IsArc(theGraph, gp_GetFirstArc(theGraph, R)) )
          {
              if (_FillVisitedFlagsInBicomp(theGraph, R, FillValue) != OK)
            	  return NOTOK;
          }
     }
     return OK;
}

/********************************************************************
 _FillVisitedFlagsInUnembeddedEdges()
 Unembedded edges aren't part of any bicomp yet, but it may be
 necessary to fill their visited flags, typically with zero.
 ********************************************************************/

void _FillVisitedFlagsInUnembeddedEdges(graphP theGraph, int FillValue)
{
int I, J;

    for (I = 0; I < theGraph->N; I++)
    {
        J = theGraph->V[I].fwdArcList;
        while (gp_IsArc(theGraph, J))
        {
            theGraph->G[J].visited =
            theGraph->G[gp_GetTwinArc(theGraph, J)].visited = FillValue;

            J = gp_GetNextArc(theGraph, J);
            if (J == theGraph->V[I].fwdArcList)
                J = NIL;
        }
    }
}

/********************************************************************
 _SetVertexTypeInBicomp()

 Sets the 'type' member to theType for each vertex in the bicomp
 rooted by BicompRoot.

 This method uses the stack but preserves whatever may have been
 on it.  In debug mode, it will return NOTOK if the stack overflows.
 This method pushes at most one integer per vertex in the bicomp.

 Returns OK on success, NOTOK on implementation failure.
 ********************************************************************/

int  _SetVertexTypeInBicomp(graphP theGraph, int BicompRoot, int theType)
{
int  V, J;
int  stackBottom = sp_GetCurrentSize(theGraph->theStack);

     sp_Push(theGraph->theStack, BicompRoot);
     while (sp_GetCurrentSize(theGraph->theStack) > stackBottom)
     {
          sp_Pop(theGraph->theStack, V);
          theGraph->G[V].type = theType;

          J = gp_GetFirstArc(theGraph, V);
          while (gp_IsArc(theGraph, J))
          {
             if (theGraph->G[J].type == EDGE_DFSCHILD)
                 sp_Push(theGraph->theStack, theGraph->G[J].v);

             J = gp_GetNextArc(theGraph, J);
          }
     }
     return OK;
}

/********************************************************************
 _ClearGraph()
 Clears all memory used by the graph, restoring it to the state it
 was in immediately after gp_New() created it.
 ********************************************************************/

void _ClearGraph(graphP theGraph)
{
     if (theGraph->G != NULL)
     {
          free(theGraph->G);
          theGraph->G = NULL;
     }
     if (theGraph->V != NULL)
     {
          free(theGraph->V);
          theGraph->V = NULL;
     }

     theGraph->N = theGraph->M = theGraph->edgeOffset = theGraph->arcCapacity = 0;
     theGraph->internalFlags = theGraph->embedFlags = 0;

     _ClearIsolatorContext(theGraph);

     LCFree(&theGraph->BicompLists);
     LCFree(&theGraph->DFSChildLists);

     sp_Free(&theGraph->theStack);

     if (theGraph->buckets != NULL)
     {
         free(theGraph->buckets);
         theGraph->buckets = NULL;
     }

     LCFree(&theGraph->bin);

     if (theGraph->extFace != NULL)
     {
         free(theGraph->extFace);
         theGraph->extFace = NULL;
     }

     sp_Free(&theGraph->edgeHoles);

     gp_FreeExtensions(theGraph);
}

/********************************************************************
 gp_Free()
 Frees G and V, then the graph record.  Then sets your pointer to NULL
 (so you must pass the address of your pointer).
 ********************************************************************/

void gp_Free(graphP *pGraph)
{
     if (pGraph == NULL) return;
     if (*pGraph == NULL) return;

     _ClearGraph(*pGraph);

     free(*pGraph);
     *pGraph = NULL;
}

/********************************************************************
 gp_CopyGraph()
 Copies the content of the srcGraph into the dstGraph.  The dstGraph
 must have been previously initialized with the same number of
 vertices as the srcGraph (e.g. gp_InitGraph(dstGraph, srcGraph->N).

 Returns OK for success, NOTOK for failure.
 ********************************************************************/

int  gp_CopyGraph(graphP dstGraph, graphP srcGraph)
{
int  I, N = srcGraph->N, edgeOffset = srcGraph->edgeOffset;
int  Gsize = edgeOffset + srcGraph->arcCapacity;

     // Parameter checks
     if (dstGraph == NULL || srcGraph == NULL)
     {
         return NOTOK;
     }

     // The graphs need to be the same order and initialized
     if (dstGraph->N != srcGraph->N || dstGraph->N == 0)
     {
         return NOTOK;
     }

     // Ensure dstGraph has the required arc capacity; this expands
     // dstGraph if needed, but does not contract.  An error is only
     // returned if the expansion fails.
     if (gp_EnsureArcCapacity(dstGraph, srcGraph->arcCapacity) != OK)
     {
    	 return NOTOK;
     }

     // Copy the basic GraphNode structures.  Augmentations to
     // the graph node structure created by extensions are copied
     // below by gp_CopyExtensions()
     for (I = 0; I < Gsize; I++)
          dstGraph->G[I] = srcGraph->G[I];

     // Copy the basic VertexRec structures.  Augmentations to
     // the vertex structure created by extensions are copied
     // below by gp_CopyExtensions()
     for (I = 0; I < N; I++)
          dstGraph->V[I] = srcGraph->V[I];

     // Copy the external face array
     for (I = 0; I < edgeOffset; I++)
     {
         dstGraph->extFace[I].vertex[0] = srcGraph->extFace[I].vertex[0];
         dstGraph->extFace[I].vertex[1] = srcGraph->extFace[I].vertex[1];
         dstGraph->extFace[I].inversionFlag = srcGraph->extFace[I].inversionFlag;
     }

     // Give the dstGraph the same size and intrinsic properties
     dstGraph->N = srcGraph->N;
     dstGraph->M = srcGraph->M;
     dstGraph->edgeOffset = srcGraph->edgeOffset;
     dstGraph->internalFlags = srcGraph->internalFlags;
     dstGraph->embedFlags = srcGraph->embedFlags;

     dstGraph->IC = srcGraph->IC;

     LCCopy(dstGraph->BicompLists, srcGraph->BicompLists);
     LCCopy(dstGraph->DFSChildLists, srcGraph->DFSChildLists);
     sp_Copy(dstGraph->theStack, srcGraph->theStack);
     sp_Copy(dstGraph->edgeHoles, srcGraph->edgeHoles);

     // Copy the set of extensions, which includes copying the
     // extension data as well as the function overload tables
     if (gp_CopyExtensions(dstGraph, srcGraph) != OK)
    	 return NOTOK;

     // Copy the graph's function table, which has the pointers to
     // the most recent extension overloads of each function (or
     // the original function pointer if a particular function has
     // not been overloaded).
     // This must be done after copying the extension because the
     // first step of copying the extensions is to delete the
     // dstGraph extensions, which clears its function table.
     // Therefore, no good to assign the srcGraph functions *before*
     // copying the extensions because the assignment would be wiped out
     // This, in turn, means that the DupContext function of an extension
     // *cannot* depend on any extension function overloads; the extension
     // must directly invoke extension functions only.
     dstGraph->functions = srcGraph->functions;

     return OK;
}

/********************************************************************
 gp_DupGraph()
 ********************************************************************/

graphP gp_DupGraph(graphP theGraph)
{
graphP result;

     if ((result = gp_New()) == NULL) return NULL;

     if (gp_InitGraph(result, theGraph->N) != OK ||
         gp_CopyGraph(result, theGraph) != OK)
     {
         gp_Free(&result);
         return NULL;
     }

     return result;
}

/********************************************************************
 gp_CreateRandomGraph()

 Creates a randomly generated graph.  First a tree is created by
 connecting each vertex to some successor.  Then a random number of
 additional random edges are added.  If an edge already exists, then
 we retry until a non-existent edge is picked.

 This function assumes the caller has already called srand().

 Returns OK on success, NOTOK on failure
 ********************************************************************/

int  gp_CreateRandomGraph(graphP theGraph)
{
int N, I, M, u, v;

     N = theGraph->N;

/* Generate a random tree; note that this method virtually guarantees
        that the graph will be renumbered, but it is linear time.
        Also, we are not generating the DFS tree but rather a tree
        that simply ensures the resulting random graph is connected. */

     for (I=1; I < N; I++)
          if (gp_AddEdge(theGraph, _GetRandomNumber(0, I-1), 0, I, 0) != OK)
              return NOTOK;

/* Generate a random number of additional edges
        (actually, leave open a small chance that no
        additional edges will be added). */

     M = _GetRandomNumber(7*N/8, theGraph->arcCapacity/2);

     if (M > N*(N-1)/2) M = N*(N-1)/2;

     for (I=N-1; I<M; I++)
     {
          u = _GetRandomNumber(0, N-2);
          v = _GetRandomNumber(u+1, N-1);

          if (gp_IsNeighbor(theGraph, u, v))
              I--;
          else
          {
              if (gp_AddEdge(theGraph, u, 0, v, 0) != OK)
                  return NOTOK;
          }
     }

     return OK;
}

/********************************************************************
 _GetRandomNumber()
 This function generates a random number between NMin and NMax
 inclusive.  It assumes that the caller has called srand().
 It calls rand(), but before truncating to the proper range,
 it adds the high bits of the rand() result into the low bits.
 The result of this is that the randomness appearing in the
 truncated bits also has an affect on the non-truncated bits.
 I used the same technique to improve the spread of hashing functions
 in my Jan.98 Dr. Dobb's Journal article  "Resizable Data Structures".
 ********************************************************************/

int  _GetRandomNumber(int NMin, int NMax)
{
int  N = rand();

     if (NMax < NMin) return NMin;

     N += ((N&0xFFFF0000)>>16);
     N += ((N&0x0000FF00)>>8);
     N &= 0x7FFFFFF;
     N %= (NMax-NMin+1);
     return N+NMin;
}

/********************************************************************
 _getUnprocessedChild()
 Support routine for gp_Create RandomGraphEx(), this function
 obtains a child of the given vertex in the randomly generated
 tree that has not yet been processed.  NIL is returned if the
 given vertex has no unprocessed children

 ********************************************************************/

int _getUnprocessedChild(graphP theGraph, int parent)
{
int J = gp_GetFirstArc(theGraph, parent);
int JTwin = gp_GetTwinArc(theGraph, J);
int child = theGraph->G[J].v;

    // The tree edges were added to the beginning of the adjacency list,
    // and we move processed tree edge records to the end of the list,
    // so if the immediate next arc (edge record) is not a tree edge
    // then we return NIL because the vertex has no remaining
    // unprocessed children
    if (theGraph->G[J].type == TYPE_UNKNOWN)
        return NIL;

    // If the child has already been processed, then all children
    // have been pushed to the end of the list, and we have just
    // encountered the first child we processed, so there are no
    // remaining unprocessed children */
    if (theGraph->G[J].visited)
        return NIL;

    // We have found an edge leading to an unprocessed child, so
    // we mark it as processed so that it doesn't get returned
    // again in future iterations.
    theGraph->G[J].visited = 1;
    theGraph->G[JTwin].visited = 1;

    // Now we move the edge record in the parent vertex to the end
    // of the adjacency list of that vertex.
    gp_MoveArcToLast(theGraph, parent, J);

    // Now we move the edge record in the child vertex to the
    // end of the adjacency list of the child.
    gp_MoveArcToLast(theGraph, child, JTwin);

    // Now we set the child's parent and return the child.
    theGraph->V[child].DFSParent = parent;

    return child;
}

/********************************************************************
 _hasUnprocessedChild()
 Support routine for gp_Create RandomGraphEx(), this function
 obtains a child of the given vertex in the randomly generated
 tree that has not yet been processed.  False (0) is returned
 unless the given vertex has an unprocessed child.
 ********************************************************************/

int _hasUnprocessedChild(graphP theGraph, int parent)
{
int J = gp_GetFirstArc(theGraph, parent);

    if (theGraph->G[J].type == TYPE_UNKNOWN)
        return 0;

    if (theGraph->G[J].visited)
        return 0;

    return 1;
}

/********************************************************************
 gp_CreateRandomGraphEx()
 Given a graph structure with a pre-specified number of vertices N,
 this function creates a graph with the specified number of edges.

 If numEdges <= 3N-6, then the graph generated is planar.  If
 numEdges is larger, then a maximal planar graph is generated, then
 (numEdges - 3N + 6) additional random edges are added.

 This function assumes the caller has already called srand().
 ********************************************************************/

int  gp_CreateRandomGraphEx(graphP theGraph, int numEdges)
{
#define EDGE_TREE_RANDOMGEN (TYPE_UNKNOWN+1)

int N, I, arc, M, root, v, c, p, last, u, J, e;

     N = theGraph->N;

     if (numEdges > theGraph->arcCapacity/2)
         numEdges = theGraph->arcCapacity/2;

/* Generate a random tree. */

    for (I=1; I < N; I++)
    {
        v = _GetRandomNumber(0, I-1);
        if (gp_AddEdge(theGraph, v, 0, I, 0) != OK)
            return NOTOK;

        else
	    {
            arc = theGraph->edgeOffset + 2*theGraph->M - 2;
		    theGraph->G[arc].type = EDGE_TREE_RANDOMGEN;
		    theGraph->G[gp_GetTwinArc(theGraph, arc)].type = EDGE_TREE_RANDOMGEN;
		    theGraph->G[arc].visited = 0;
		    theGraph->G[gp_GetTwinArc(theGraph, arc)].visited = 0;
	    }
    }

/* Add edges up to the limit or until the graph is maximal planar. */

    M = numEdges <= 3*N - 6 ? numEdges : 3*N - 6;

    root = 0;
    v = last = _getUnprocessedChild(theGraph, root);

    while (v != root && theGraph->M < M)
    {
	     c = _getUnprocessedChild(theGraph, v);

	     if (c != NIL)
	     {
             if (last != v)
             {
		         if (gp_AddEdge(theGraph, last, 1, c, 1) != OK)
			         return NOTOK;
             }

		     if (gp_AddEdge(theGraph, root, 1, c, 1) != OK)
			     return NOTOK;

		     v = last = c;
	     }

	     else
	     {
		     p = theGraph->V[v].DFSParent;
		     while (p != NIL && (c = _getUnprocessedChild(theGraph, p)) == NIL)
		     {
			     v = p;
			     p = theGraph->V[v].DFSParent;
			     if (p != NIL && p != root)
			     {
				     if (gp_AddEdge(theGraph, last, 1, p, 1) != OK)
					     return NOTOK;
			     }
		     }

		     if (p != NIL)
		     {
                 if (p == root)
                 {
                     if (gp_AddEdge(theGraph, v, 1, c, 1) != OK)
				         return NOTOK;

                     if (v != last)
                     {
			             if (gp_AddEdge(theGraph, last, 1, c, 1) != OK)
				             return NOTOK;
                     }
                 }
                 else
                 {
			         if (gp_AddEdge(theGraph, last, 1, c, 1) != OK)
				         return NOTOK;
                 }

                 if (p != root)
                 {
			        if (gp_AddEdge(theGraph, root, 1, c, 1) != OK)
				         return NOTOK;
                    last = c;
                 }

			     v = c;
		     }
	     }
    }

/* Add additional edges if the limit has not yet been reached. */

    while (theGraph->M < numEdges)
    {
        u = _GetRandomNumber(0, N-1);
        v = _GetRandomNumber(0, N-1);

        if (u != v && !gp_IsNeighbor(theGraph, u, v))
            if (gp_AddEdge(theGraph, u, 0, v, 0) != OK)
                return NOTOK;
    }

/* Clear the edge types back to 'unknown' */

    for (e = 0; e < numEdges; e++)
    {
        J = theGraph->edgeOffset + 2*e;
        theGraph->G[J].type = TYPE_UNKNOWN;
        theGraph->G[gp_GetTwinArc(theGraph, J)].type = TYPE_UNKNOWN;
        theGraph->G[J].visited = 0;
        theGraph->G[gp_GetTwinArc(theGraph, J)].visited = 0;
    }

/* Put all DFSParent indicators back to NIL */

    for (I = 0; I < N; I++)
        theGraph->V[I].DFSParent = NIL;

    return OK;

#undef EDGE_TREE_RANDOMGEN
}

/********************************************************************
 gp_SetDirection()
 Behavior depends on edgeFlag_Direction (EDGEFLAG_DIRECTION_INONLY,
 EDGEFLAG_DIRECTION_OUTONLY, or 0).
 A direction of 0 clears directedness. Otherwise, edge record e is set
 to edgeFlag_Direction and e's twin arc is set to the opposing setting.
 ********************************************************************/

void gp_SetDirection(graphP theGraph, int e, int edgeFlag_Direction)
{
	int eTwin = gp_GetTwinArc(theGraph, e);

	if (edgeFlag_Direction == EDGEFLAG_DIRECTION_INONLY)
	{
		theGraph->G[e].flags |= EDGEFLAG_DIRECTION_INONLY;
		theGraph->G[eTwin].flags |= EDGEFLAG_DIRECTION_OUTONLY;
	}
	else if (edgeFlag_Direction == EDGEFLAG_DIRECTION_OUTONLY)
	{
		theGraph->G[e].flags |= EDGEFLAG_DIRECTION_OUTONLY;
		theGraph->G[eTwin].flags |= EDGEFLAG_DIRECTION_INONLY;
	}
	else
	{
		theGraph->G[e].flags &= ~(EDGEFLAG_DIRECTION_INONLY|EDGEFLAG_DIRECTION_OUTONLY);
		theGraph->G[eTwin].flags &= ~(EDGEFLAG_DIRECTION_INONLY|EDGEFLAG_DIRECTION_OUTONLY);
	}
}

/********************************************************************
 gp_IsNeighbor()

 Checks whether v is already in u's adjacency list, i.e. does the arc
 u -> v exist.
 If there is an edge record for v in u's list, but it is marked INONLY,
 then it represents the arc v->u but not u->v, so it is ignored.

 Returns TRUE or FALSE.
 ********************************************************************/

int  gp_IsNeighbor(graphP theGraph, int u, int v)
{
int  J;

     J = gp_GetFirstArc(theGraph, u);
     while (gp_IsArc(theGraph, J))
     {
          if (theGraph->G[J].v == v)
          {
              if (!gp_GetDirection(theGraph, J, EDGEFLAG_DIRECTION_INONLY))
            	  return TRUE;
          }
          J = gp_GetNextArc(theGraph, J);
     }
     return FALSE;
}

/********************************************************************
 gp_GetNeighborEdgeRecord()
 Searches the adjacency list of u to obtains the edge record for v.

 NOTE: The caller should check whether the edge record is INONLY;
       This method returns any edge record representing a connection
       between vertices u and v, so this method can return an
       edge record even if gp_IsNeighbor(theGraph, u, v) is false (0).
       To filter out INONLY edge records, use gp_GetDirection() on
       the edge record returned by this method.

 Returns NIL if there is no edge record indicating v in u's adjacency
         list, or the edge record location otherwise.
 ********************************************************************/

int  gp_GetNeighborEdgeRecord(graphP theGraph, int u, int v)
{
int  J;

     if (u == NIL || v == NIL)
    	 return NIL + NOTOK;

     J = gp_GetFirstArc(theGraph, u);
     while (gp_IsArc(theGraph, J))
     {
          if (theGraph->G[J].v == v)
        	  return J;

          J = gp_GetNextArc(theGraph, J);
     }
     return NIL;
}

/********************************************************************
 gp_GetVertexDegree()

 Counts the number of edge records in the adjacency list of a given
 vertex V.  The while loop condition is 2N or higher because our
 data structure keeps records at locations 0 to N-1 for vertices
 AND N to 2N-1 for copies of vertices.  So edge records are stored
 at locations 2N and above.

 Note: For digraphs, this method returns the total degree of the
       vertex, including outward arcs (undirected and OUTONLY)
       as well as INONLY arcs.  Other functions are defined to get
       the in-degree or out-degree of the vertex.

 Note: This function determines the degree by counting.  An extension
       could cache the degree value of each vertex and update the
       cached value as edges are added and deleted.
 ********************************************************************/

int  gp_GetVertexDegree(graphP theGraph, int v)
{
int  J, degree;

     if (theGraph==NULL || v==NIL)
    	 return NOTOK;

     degree = 0;

     J = gp_GetFirstArc(theGraph, v);
     while (gp_IsArc(theGraph, J))
     {
         degree++;
         J = gp_GetNextArc(theGraph, J);
     }

     return degree;
}

/********************************************************************
 gp_GetVertexInDegree()

 Counts the number of edge records in the adjacency list of a given
 vertex V that represent arcs from another vertex into V.
 This includes undirected edges and INONLY arcs, so it only excludes
 edges records that are marked as OUTONLY arcs.

 Note: This function determines the in-degree by counting.  An extension
       could cache the in-degree value of each vertex and update the
       cached value as edges are added and deleted.
 ********************************************************************/

int  gp_GetVertexInDegree(graphP theGraph, int v)
{
int  J, degree;

     if (theGraph==NULL || v==NIL)
    	 return NOTOK;

     degree = 0;

     J = gp_GetFirstArc(theGraph, v);
     while (gp_IsArc(theGraph, J))
     {
         if (!gp_GetDirection(theGraph, J, EDGEFLAG_DIRECTION_OUTONLY))
             degree++;
         J = gp_GetNextArc(theGraph, J);
     }

     return degree;
}

/********************************************************************
 gp_GetVertexOutDegree()

 Counts the number of edge records in the adjacency list of a given
 vertex V that represent arcs from V to another vertex.
 This includes undirected edges and OUTONLY arcs, so it only excludes
 edges records that are marked as INONLY arcs.

 Note: This function determines the out-degree by counting.  An extension
       could cache the out-degree value of each vertex and update the
       cached value as edges are added and deleted.
 ********************************************************************/

int  gp_GetVertexOutDegree(graphP theGraph, int v)
{
int  J, degree;

     if (theGraph==NULL || v==NIL)
    	 return NOTOK;

     degree = 0;

     J = gp_GetFirstArc(theGraph, v);
     while (gp_IsArc(theGraph, J))
     {
         if (!gp_GetDirection(theGraph, J, EDGEFLAG_DIRECTION_INONLY))
             degree++;
         J = gp_GetNextArc(theGraph, J);
     }

     return degree;
}

/********************************************************************
 gp_AttachArc()

 This routine adds newArc into v's adjacency list at a position
 adjacent to the edge record for e, either before or after e,
 depending on link.  If e is not an arc (e.g. if e is NIL),
 then link is assumed to indicate whether the new arc is to be
 placed at the beginning or end of v's adjacency list.

 NOTE: The caller can pass NIL for v if e is not NIL, since the
       vertex is implied (theGraph->G[eTwin].v)

 The arc is assumed to already exist in the data structure (i.e.
 the storage of edges), as only a whole edge (two arcs) can be
 inserted into or deleted from the data structure.  Hence there is
 no such thing as gp_InsertArc() or gp_DeleteArc().

 See also gp_DetachArc(), gp_InsertEdge() and gp_DeleteEdge()
 ********************************************************************/

void gp_AttachArc(graphP theGraph, int v, int e, int link, int newArc)
{
     if (gp_IsArc(theGraph, e))
     {
    	 int e2 = gp_GetAdjacentArc(theGraph, e, link);

         // e's link is newArc, and newArc's 1^link is e
    	 gp_SetAdjacentArc(theGraph, e, link, newArc);
    	 gp_SetAdjacentArc(theGraph, newArc, 1^link, e);

    	 // newArcs's link is e2
    	 gp_SetAdjacentArc(theGraph, newArc, link, e2);

    	 // if e2 is an arc, then e2's 1^link is newArc, else v's 1^link is newArc
    	 if (gp_IsArc(theGraph, e2))
    		 gp_SetAdjacentArc(theGraph, e2, 1^link, newArc);
    	 else
    		 gp_SetArc(theGraph, v, 1^link, newArc);
     }
     else
     {
    	 int e2 = gp_GetArc(theGraph, v, link);

    	 // v's link is newArc, and newArc's 1^link is gp_AdjacencyListEndMark(v)
    	 gp_SetArc(theGraph, v, link, newArc);
    	 gp_SetAdjacentArc(theGraph, newArc, 1^link, gp_AdjacencyListEndMark(v));

    	 // newArcs's elink is e2
    	 gp_SetAdjacentArc(theGraph, newArc, link, e2);

    	 // if e2 is an arc, then e2's 1^link is newArc, else v's 1^link is newArc
    	 if (gp_IsArc(theGraph, e2))
    		 gp_SetAdjacentArc(theGraph, e2, 1^link, newArc);
    	 else
    		 gp_SetArc(theGraph, v, 1^link, newArc);
     }
}

/****************************************************************************
 gp_DetachArc()

 This routine detaches arc from its adjacency list, but it does not delete
 it from the data structure (only a whole edge can be deleted).

 Some algorithms must temporarily detach an edge, perform some calculation,
 and eventually put the edge back. This routine supports that operation.
 The neighboring adjacency list nodes are cross-linked, but the two link
 members of the arc are retained, so the arc can be reattached later by
 invoking _RestoreArc().  A sequence of detached arcs can only be restored
 in the exact opposite order of their detachment.  Thus, algorithms do not
 directly use this method to implement the temporary detach/restore method.
 Instead, gp_HideEdge() and gp_RestoreEdge are used, and algorithms push
 edge hidden edge onto the stack.  One example of this stack usage is
 provided by detaching edges with gp_ContractEdge() or gp_IdentifyVertices(),
 and reattaching with gp_RestoreIdentifications(), which unwinds the stack
 by invoking gp_RestoreVertex().
 ****************************************************************************/

void gp_DetachArc(graphP theGraph, int arc)
{
	int nextArc = gp_GetNextArc(theGraph, arc),
	    prevArc = gp_GetPrevArc(theGraph, arc);

	    if (gp_IsArc(theGraph, nextArc))
	    	gp_SetPrevArc(theGraph, nextArc, prevArc);
	    else
	    	gp_SetLastArc(theGraph, theGraph->G[gp_GetTwinArc(theGraph, arc)].v, prevArc);

	    if (gp_IsArc(theGraph, prevArc))
	    	gp_SetNextArc(theGraph, prevArc, nextArc);
	    else
	    	gp_SetFirstArc(theGraph, theGraph->G[gp_GetTwinArc(theGraph, arc)].v, nextArc);
}

/********************************************************************
 gp_AddEdge()
 Adds the undirected edge (u,v) to the graph by placing edge records
 representing u into v's circular edge record list and v into u's
 circular edge record list.

 upos receives the location in G where the u record in v's list will be
        placed, and vpos is the location in G of the v record we placed in
 u's list.  These are used to initialize the short circuit links.

 ulink (0|1) indicates whether the edge record to v in u's list should
        become adjacent to u by its 0 or 1 link, i.e. u[ulink] == vpos.
 vlink (0|1) indicates whether the edge record to u in v's list should
        become adjacent to v by its 0 or 1 link, i.e. v[vlink] == upos.

 ********************************************************************/

int  gp_AddEdge(graphP theGraph, int u, int ulink, int v, int vlink)
{
int  upos, vpos;

     if (theGraph==NULL || u<0 || v<0 ||
         u>=2*theGraph->N || v>=2*theGraph->N)
         return NOTOK;

     /* We enforce the edge limit */

     if (theGraph->M >= theGraph->arcCapacity/2)
         return NONEMBEDDABLE;

     if (sp_NonEmpty(theGraph->edgeHoles))
     {
         sp_Pop(theGraph->edgeHoles, vpos);
     }
     else
         vpos = theGraph->edgeOffset + 2*theGraph->M;

     upos = gp_GetTwinArc(theGraph, vpos);

     theGraph->G[upos].v = v;
     gp_AttachArc(theGraph, u, NIL, ulink, upos);
     theGraph->G[vpos].v = u;
     gp_AttachArc(theGraph, v, NIL, vlink, vpos);

     theGraph->M++;
     return OK;
}

/********************************************************************
 gp_InsertEdge()

 This function adds the edge (u, v) such that the edge record added
 to the adjacency list of u is adjacent to e_u and the edge record
 added to the adjacency list of v is adjacent to e_v.
 The direction of adjacency is given by e_ulink for e_u and e_vlink
 for e_v. Specifically, the new edge will be comprised of two arcs,
 n_u and n_v.  In u's (v's) adjacency list, n_u (n_v) will be added
 so that it is indicated by e_u's (e_v's) e_ulink (e_vlink).
 If e_u (or e_v) is not an arc, then e_ulink (e_vlink) indicates
 whether to prepend or append to the adjacency list for u (v).
 ********************************************************************/

int  gp_InsertEdge(graphP theGraph, int u, int e_u, int e_ulink,
                                    int v, int e_v, int e_vlink)
{
int vertMax = 2*theGraph->N - 1,
    edgeMax = theGraph->edgeOffset + 2*theGraph->M + 2*sp_GetCurrentSize(theGraph->edgeHoles) - 1,
    upos, vpos;

     if (theGraph==NULL || u<0 || v<0 || u>vertMax || v>vertMax ||
         e_u>edgeMax || (e_u<theGraph->edgeOffset && e_u != gp_AdjacencyListEndMark(u)) ||
         e_v>edgeMax || (e_v<theGraph->edgeOffset && e_v != gp_AdjacencyListEndMark(v)) ||
         e_ulink<0 || e_ulink>1 || e_vlink<0 || e_vlink>1)
         return NOTOK;

     if (theGraph->M >= theGraph->arcCapacity/2)
         return NONEMBEDDABLE;

     if (sp_NonEmpty(theGraph->edgeHoles))
     {
         sp_Pop(theGraph->edgeHoles, vpos);
     }
     else
         vpos = theGraph->edgeOffset + 2*theGraph->M;

     upos = gp_GetTwinArc(theGraph, vpos);

     theGraph->G[upos].v = v;
     gp_AttachArc(theGraph, u, e_u, e_ulink, upos);

     theGraph->G[vpos].v = u;
     gp_AttachArc(theGraph, v, e_v, e_vlink, vpos);

     theGraph->M++;

     return OK;
}

/****************************************************************************
 gp_DeleteEdge()

 This function deletes the given edge record e and its twin, reducing the
 number of edges M in the graph.
 Before the e^th record is deleted, its 'nextLink' adjacency list neighbor
 is collected as the return result.  This is useful when iterating through
 an edge list and making deletions because the nextLink arc is the 'next'
 arc in the iteration, but it is hard to obtain *after* deleting e.
 ****************************************************************************/

int  gp_DeleteEdge(graphP theGraph, int e, int nextLink)
{
int  eTwin = gp_GetTwinArc(theGraph, e);
int  M = theGraph->M;
int  nextArc, JPos, MPos;

/* Calculate the nextArc after e so that, when e is deleted, the return result
        informs a calling loop of the next edge to be processed. */

     nextArc = gp_GetAdjacentArc(theGraph, e, nextLink);

/* Delete the edge records J and JTwin from their adjacency lists. */

     gp_DetachArc(theGraph, e);
     gp_DetachArc(theGraph, eTwin);

/* Clear the edge record contents */

    theGraph->functions.fpInitGraphNode(theGraph, e);
    theGraph->functions.fpInitGraphNode(theGraph, eTwin);

/* If records e and eTwin are not the last in the edge record array, then
     we want to record a new hole in the edge array. */

     JPos = (e < eTwin ? e : eTwin);
     MPos = theGraph->edgeOffset + 2*(M-1) + 2*sp_GetCurrentSize(theGraph->edgeHoles);

     if (JPos < MPos)
     {
         sp_Push(theGraph->edgeHoles, JPos);
     }

/* Now we reduce the number of edges in the data structure, and then
        return the previously calculated successor of J. */

     theGraph->M--;
     return nextArc;
}

/********************************************************************
 _RestoreArc()
 This routine reinserts an arc into the edge list from which it
 was previously removed by gp_DetachArc().

 The assumed processing model is that arcs will be restored in reverse
 of the order in which they were hidden, i.e. it is assumed that the
 hidden arcs will be pushed on a stack and the arcs will be popped
 from the stack for restoration.
 ********************************************************************/

void _RestoreArc(graphP theGraph, int arc)
{
int nextArc = gp_GetNextArc(theGraph, arc),
	prevArc = gp_GetPrevArc(theGraph, arc);

	if (gp_IsArc(theGraph, nextArc))
		gp_SetPrevArc(theGraph, nextArc, arc);
	else
		gp_SetLastArc(theGraph, theGraph->G[gp_GetTwinArc(theGraph, arc)].v, arc);

    if (gp_IsArc(theGraph, prevArc))
    	gp_SetNextArc(theGraph, prevArc, arc);
    else
    	gp_SetFirstArc(theGraph, theGraph->G[gp_GetTwinArc(theGraph, arc)].v, arc);
}

/********************************************************************
 gp_HideEdge()
 This routine removes the two arcs of an edge from the adjacency lists
 of its endpoint vertices, but does not delete them from the storage
 data structure.

 Many algorithms must temporarily remove an edge, perform some
 calculation, and eventually put the edge back. This routine supports
 that operation.

 For each arc, the neighboring adjacency list nodes are cross-linked,
 but the links in the arc are retained because they indicate the
 neighbor arcs to which the arc can be reattached by gp_RestoreEdge().
 ********************************************************************/

void gp_HideEdge(graphP theGraph, int e)
{
	theGraph->functions.fpHideEdge(theGraph, e);
}

void _HideEdge(graphP theGraph, int e)
{
	gp_DetachArc(theGraph, e);
	gp_DetachArc(theGraph, gp_GetTwinArc(theGraph, e));
}

/********************************************************************
 gp_RestoreEdge()
 This routine reinserts two two arcs of an edge into the adjacency
 lists of the edge's endpoints, the arcs having been previously
 removed by gp_HideEdge().

 The assumed processing model is that edges will be restored in
 reverse of the order in which they were hidden, i.e. it is assumed
 that the hidden edges will be pushed on a stack and the edges will
 be popped from the stack for restoration.

 Note: Since both arcs of an edge are restored, only one arc need
        be pushed on the stack for restoration.  This routine
        restores the two arcs in the opposite order from the order
        in which they are hidden by gp_HideEdge().
 ********************************************************************/

void gp_RestoreEdge(graphP theGraph, int e)
{
	theGraph->functions.fpRestoreEdge(theGraph, e);
}

void _RestoreEdge(graphP theGraph, int e)
{
     _RestoreArc(theGraph, gp_GetTwinArc(theGraph, e));
     _RestoreArc(theGraph, e);
}

/********************************************************************
 _HideInternalEdges()
 Pushes onto the graph's stack and hides all arc nodes of the vertex
 except the first and last arcs in the adjacency list of the vertex.
 This method is typically called on a vertex that is on the external
 face of a biconnected component, because the first and last arcs are
 the ones that attach the vertex to the external face cycle, and any
 other arcs in the adjacency list are inside that cycle.

 This method uses the stack. The caller is expected to clear the stack
 or save the stack size before invocation, since the stack size is
 needed to _RestoreInternalEdges().
 ********************************************************************/

int  _HideInternalEdges(graphP theGraph, int vertex)
{
int J = gp_GetFirstArc(theGraph, vertex);

    // If the vertex adjacency list is empty or if it contains
    // only one edge, then there are no *internal* edges to hide
    if (J == gp_GetLastArc(theGraph, vertex))
    	return OK;

    // Start with the first internal edge
    J = gp_GetNextArc(theGraph, J);

    // Cycle through all the edges, pushing each except stop
    // before pushing the last edge, which is not internal
    while (J != gp_GetLastArc(theGraph, vertex))
    {
        sp_Push(theGraph->theStack, J);
        gp_HideEdge(theGraph, J);
        J = gp_GetNextArc(theGraph, J);
    }

    return OK;
}

/********************************************************************
 _RestoreInternalEdges()
 Reverses the effects of _HideInternalEdges()
 ********************************************************************/

int  _RestoreInternalEdges(graphP theGraph, int stackBottom)
{
	return _RestoreHiddenEdges(theGraph, stackBottom);
}

/********************************************************************
 _RestoreHiddenEdges()

 Each entry on the stack, down to stackBottom, is assumed to be an
 edge record (arc) pushed in concert with invoking gp_HideEdge().
 Each edge is restored using gp_RestoreEdge() in exact reverse of the
 hiding order.  The stack is reduced in size to stackBottom.

 Returns OK on success, NOTOK on internal failure.
 ********************************************************************/

int  _RestoreHiddenEdges(graphP theGraph, int stackBottom)
{
	int  e;

	 while (sp_GetCurrentSize(theGraph->theStack) > stackBottom)
	 {
		  sp_Pop(theGraph->theStack, e);
		  if (!gp_IsArc(theGraph, e))
			  return NOTOK;
		  gp_RestoreEdge(theGraph, e);
	 }

	 return OK;
}

/********************************************************************
 gp_HideVertex()

 Pushes onto the graph's stack and hides all arc nodes of the vertex.
 Additional integers are then pushed so that the result is reversible
 by gp_RestoreVertex().  See that method for details on the expected
 stack segment.

 Returns OK for success, NOTOK for internal failure.
 ********************************************************************/

int  gp_HideVertex(graphP theGraph, int vertex)
{
	if (vertex == NIL)
		return NOTOK;

	return theGraph->functions.fpHideVertex(theGraph, vertex);
}

int  _HideVertex(graphP theGraph, int vertex)
{
	int hiddenEdgeStackBottom = sp_GetCurrentSize(theGraph->theStack);
	int J = gp_GetFirstArc(theGraph, vertex);

    // Cycle through all the edges, pushing and hiding each
    while (gp_IsArc(theGraph, J))
    {
        sp_Push(theGraph->theStack, J);
        gp_HideEdge(theGraph, J);
        J = gp_GetNextArc(theGraph, J);
    }

    // Push the additional integers needed by gp_RestoreVertex()
	sp_Push(theGraph->theStack, hiddenEdgeStackBottom);
	sp_Push(theGraph->theStack, NIL);
	sp_Push(theGraph->theStack, NIL);
	sp_Push(theGraph->theStack, NIL);
	sp_Push(theGraph->theStack, NIL);
	sp_Push(theGraph->theStack, NIL);
	sp_Push(theGraph->theStack, vertex);

    return OK;
}

/********************************************************************
 gp_ContractEdge()

 Contracts the edge e=(u,v).  This hides the edge (both e and its
 twin arc), and it also identifies vertex v with u.
 See gp_IdentifyVertices() for further details.

 Returns OK for success, NOTOK for internal failure.
 ********************************************************************/

int gp_ContractEdge(graphP theGraph, int e)
{
	if (!gp_IsArc(theGraph, e))
		return NOTOK;

	return theGraph->functions.fpContractEdge(theGraph, e);
}

int _ContractEdge(graphP theGraph, int e)
{
	int eBefore, u, v;

	if (!gp_IsArc(theGraph, e))
		return NOTOK;

	u = theGraph->G[gp_GetTwinArc(theGraph, e)].v;
	v = theGraph->G[e].v;

	eBefore = gp_GetNextArc(theGraph, e);
	sp_Push(theGraph->theStack, e);
	gp_HideEdge(theGraph, e);

	return gp_IdentifyVertices(theGraph, u, v, eBefore);
}

/********************************************************************
 gp_IdentifyVertices()

 Identifies vertex v with vertex u by transferring all adjacencies
 of v to u.  Any duplicate edges are removed as described below.
 The non-duplicate edges of v are added to the adjacency list of u
 without disturbing their relative order, and they are added before
 the edge record eBefore in u's list. If eBefore is NIL, then the
 edges are simply appended to u's list.

 If u and v are adjacent, then gp_HideEdge() is invoked to remove
 the edge e=(u,v). Then, the edges of v that indicate neighbors of
 u are also hidden.  This is done by setting the visited flags of
 u's neighbors, then traversing the adjacency list of v.  For each
 visited neighbor of v, the edge is hidden because it would duplicate
 an adjacency already expressed in u's list. Finally, the remaining
 edges of v are moved to u's list, and each twin arc is adjusted
 to indicate u as a neighbor rather than v.

 This routine assumes that the visited flags are clear beforehand,
 and visited flag settings made herein are cleared before returning.

 The following are pushed, in order, onto the graph's built-in stack:
 1) an integer for each hidden edge
 2) the stack size before any hidden edges were pushed
 3) six integers that indicate u, v and the edges moved from v to u

 An algorithm that identifies a series of vertices, either through
 directly calling this method or via gp_ContractEdge(), can unwind
 the identifications using gp_RestoreIdentifications(), which
 invokes gp_RestoreVertex() repeatedly.

 Returns OK on success, NOTOK on internal failure
 ********************************************************************/

int gp_IdentifyVertices(graphP theGraph, int u, int v, int eBefore)
{
	return theGraph->functions.fpIdentifyVertices(theGraph, u, v, eBefore);
}

int _IdentifyVertices(graphP theGraph, int u, int v, int eBefore)
{
	int e = gp_GetNeighborEdgeRecord(theGraph, u, v);
	int hiddenEdgeStackBottom, eBeforePred, J;

	// If the vertices are adjacent, then the identification is
	// essentially an edge contraction with a bit of fixup.
	if (gp_IsArc(theGraph, e))
	{
		int result = gp_ContractEdge(theGraph, e);

		// The edge contraction operation pushes one hidden edge then
	    // recursively calls this method. This method then pushes K
	    // hidden edges then an integer indicating where the top of
	    // stack was before the edges were hidden. That integer
	    // indicator must be decremented, thereby incrementing the
		// number of hidden edges to K+1.
		// After pushing the K hidden edges and the stackBottom of
		// the hidden edges, the recursive call to this method pushes
		// six more integers to indicate edges that were moved from
		// v to u, so the "hidden edges stackBottom" is in the next
		// position down.
		int hiddenEdgesStackBottomIndex = sp_GetCurrentSize(theGraph->theStack)-7;
		int hiddenEdgesStackBottomValue = sp_Get(theGraph->theStack, hiddenEdgesStackBottomIndex);

		sp_Set(theGraph->theStack, hiddenEdgesStackBottomIndex,	hiddenEdgesStackBottomValue - 1);

		return result;
	}

	// Now, u and v are not adjacent. Before we do any edge hiding or
	// moving, we record the current stack size, as this is the
	// stackBottom for the edges that will be hidden next.
	hiddenEdgeStackBottom = sp_GetCurrentSize(theGraph->theStack);

	// Mark as visited all neighbors of u
    J = gp_GetFirstArc(theGraph, u);
    while (gp_IsArc(theGraph, J))
    {
    	 if (theGraph->G[theGraph->G[J].v].visited != 0)
    		 return NOTOK;

         theGraph->G[theGraph->G[J].v].visited = 1;
         J = gp_GetNextArc(theGraph, J);
    }

	// For each edge record of v, if the neighbor is visited, then
	// push and hide the edge.
    J = gp_GetFirstArc(theGraph, v);
    while (gp_IsArc(theGraph, J))
    {
         if (theGraph->G[theGraph->G[J].v].visited)
         {
             sp_Push(theGraph->theStack, J);
             gp_HideEdge(theGraph, J);
         }
         J = gp_GetNextArc(theGraph, J);
    }

	// Mark as unvisited all neighbors of u
    J = gp_GetFirstArc(theGraph, u);
    while (gp_IsArc(theGraph, J))
    {
         theGraph->G[theGraph->G[J].v].visited = 0;
         J = gp_GetNextArc(theGraph, J);
    }

	// Push the hiddenEdgeStackBottom as a record of how many hidden
	// edges were pushed (also, see above for Contract Edge adjustment)
	sp_Push(theGraph->theStack, hiddenEdgeStackBottom);

	// Moving v's adjacency list to u is aided by knowing the predecessor
	// of u's eBefore (the edge record in u's list before which the
	// edge records of v will be added).
	eBeforePred = gp_IsArc(theGraph, eBefore)
	              ? gp_GetPrevArc(theGraph, eBefore)
	              : gp_GetLastArc(theGraph, u);

	// Turns out we only need to record six integers related to the edges
	// being moved in order to easily restore them later.
	sp_Push(theGraph->theStack, eBefore);
	sp_Push(theGraph->theStack, gp_GetLastArc(theGraph, v));
	sp_Push(theGraph->theStack, gp_GetFirstArc(theGraph, v));
	sp_Push(theGraph->theStack, eBeforePred);
	sp_Push(theGraph->theStack, u);
	sp_Push(theGraph->theStack, v);

	// For the remaining edge records of v, reassign the 'v' member
	//    of each twin arc to indicate u rather than v.
    J = gp_GetFirstArc(theGraph, v);
    while (gp_IsArc(theGraph, J))
    {
         theGraph->G[gp_GetTwinArc(theGraph, J)].v = u;
         J = gp_GetNextArc(theGraph, J);
    }

    // If v has any edges left after hiding edges, indicating common neighbors with u, ...
    if (gp_IsArc(theGraph, gp_GetFirstArc(theGraph, v)))
    {
    	// Then perform the list union of v into u between eBeforePred and eBefore
        if (gp_IsArc(theGraph, eBeforePred))
        {
        	if (gp_IsArc(theGraph, gp_GetFirstArc(theGraph, v)))
        	{
            	gp_SetNextArc(theGraph, eBeforePred, gp_GetFirstArc(theGraph, v));
            	gp_SetPrevArc(theGraph, gp_GetFirstArc(theGraph, v), eBeforePred);
        	}
        }
        else
        {
        	gp_SetFirstArc(theGraph, u, gp_GetFirstArc(theGraph, v));
        }

        if (gp_IsArc(theGraph, eBefore))
        {
        	if (gp_IsArc(theGraph, gp_GetLastArc(theGraph, v)))
        	{
            	gp_SetNextArc(theGraph, gp_GetLastArc(theGraph, v), eBefore);
            	gp_SetPrevArc(theGraph, eBefore, gp_GetLastArc(theGraph, v));
        	}
        }
        else
        {
        	gp_SetLastArc(theGraph, u, gp_GetLastArc(theGraph, v));
        }

        gp_SetFirstArc(theGraph, v, NIL);
        gp_SetLastArc(theGraph, v, NIL);
    }

    return OK;
}

/********************************************************************
 gp_RestoreVertex()

 This method assumes the built-in graph stack contents are the result
 of vertex hide, vertex identify and edge contract operations.
 This content consists of segments of integers, each segment
 corresponding to the removal of a vertex during an edge contraction
 or vertex identification in which a vertex v was merged into a
 vertex u.  The segment contains two blocks of integers.
 The first block contains information about u, v, the edge records
 in v's adjacency list that were added to u, and where in u's
 adjacency list they were added.  The second block of integers
 contains a list of edges incident to v that were hidden from the
 graph because they were incident to neighbors of v that were also
 neighbors of u (so they would have produced duplicate edges had
 they been left in v's adjacency list when it was merged with u's
 adjacency list).

 This method pops the first block of the segment off the stack and
 uses the information to help remove v's adjacency list from u and
 restore it into v.  Then, the second block is removed from the
 stack, and each indicated edge is restored from the hidden state.

 It is anticipated that this method will be overloaded by extension
 algorithms to perform some processing as each vertex is restored.
 Before restoration, the topmost segment has the following structure:

 ... FHE ... LHE HESB e_u_succ e_v_last e_v_first e_u_pred u v
      ^------------|

 FHE = First hidden edge
 LHE = Last hidden edge
 HESB = Hidden edge stack bottom
 e_u_succ, e_u_pred = The edges of u between which the edges of v
                      were inserted. NIL can appear if the edges of v
                      were added to the beginning or end of u's list
 e_v_first, e_v_last = The first and last edges of v's list, once
                       the hidden edges were removed

 Returns OK for success, NOTOK for internal failure.
 ********************************************************************/

int gp_RestoreVertex(graphP theGraph)
{
	return theGraph->functions.fpRestoreVertex(theGraph);
}

int _RestoreVertex(graphP theGraph)
{
int u, v, e_u_succ, e_u_pred, e_v_first, e_v_last, HESB, J;

    if (sp_GetCurrentSize(theGraph->theStack) < 7)
    	return NOTOK;

    sp_Pop(theGraph->theStack, v);
    sp_Pop(theGraph->theStack, u);
	sp_Pop(theGraph->theStack, e_u_pred);
	sp_Pop(theGraph->theStack, e_v_first);
	sp_Pop(theGraph->theStack, e_v_last);
	sp_Pop(theGraph->theStack, e_u_succ);

	// If u is not NIL, then vertex v was identified with u.  Otherwise, v was
	// simply hidden, so we skip to restoring the hidden edges.
	if (u != NIL)
	{
		// Remove v's adjacency list from u, including accounting for degree 0 case
		if (gp_IsArc(theGraph, e_u_pred))
		{
			gp_SetNextArc(theGraph, e_u_pred, e_u_succ);
			// If the successor edge exists, link it to the predecessor,
			// otherwise the predecessor is the new last arc
			if (gp_IsArc(theGraph, e_u_succ))
				gp_SetPrevArc(theGraph, e_u_succ, e_u_pred);
			else
				gp_SetLastArc(theGraph, u, e_u_pred);
		}
		else if (gp_IsArc(theGraph, e_u_succ))
		{
			// The successor arc exists, but not the predecessor,
			// so the successor is the new first arc
			gp_SetPrevArc(theGraph, e_u_succ, NIL);
			gp_SetFirstArc(theGraph, u, e_u_succ);
		}
		else
		{
			// Just in case u was degree zero
			gp_SetFirstArc(theGraph, u, NIL);
			gp_SetLastArc(theGraph, u, NIL);
		}

		// Place v's adjacency list into v, including accounting for degree 0 case
		gp_SetFirstArc(theGraph, v, e_v_first);
		gp_SetLastArc(theGraph, v, e_v_last);
		if (gp_IsArc(theGraph, e_v_first))
			gp_SetPrevArc(theGraph, e_v_first, NIL);
		if (gp_IsArc(theGraph, e_v_last))
			gp_SetPrevArc(theGraph, e_v_last, NIL);

		// For each edge record restored to v's adjacency list, reassign the 'v' member
		//    of each twin arc to indicate v rather than u.
	    J = e_v_first;
	    while (gp_IsArc(theGraph, J))
	    {
	         theGraph->G[gp_GetTwinArc(theGraph, J)].v = v;

	         if (J == e_v_last)
	        	 J = NIL;
	         else
	        	 J = gp_GetNextArc(theGraph, J);
	    }
	}

	// Restore the hidden edges of v, if any
	sp_Pop(theGraph->theStack, HESB);
	return _RestoreHiddenEdges(theGraph, HESB);
}

/********************************************************************
 gp_RestoreVertices()

 This method assumes the built-in graph stack has content consistent
 with numerous vertex identification or edge contraction operations.
 This method unwinds the stack, moving edges back to their original
 vertex owners and restoring hidden edges.
 This method is a simple iterator that invokes gp_RestoreVertex()
 until the stack is empty, so extension algorithms are more likely
 to overload gp_RestoreVertex().

 Returns OK for success, NOTOK for internal failure.
 ********************************************************************/

int gp_RestoreVertices(graphP theGraph)
{
    while (sp_NonEmpty(theGraph->theStack))
    {
    	if (gp_RestoreVertex(theGraph) != OK)
    		return NOTOK;
    }

    return OK;
}

/****************************************************************************
 _ComputeArcType()
 This is just a little helper function that automates a sequence of decisions
 that has to be made a number of times.
 An edge record is being added to the adjacency list of a; it indicates that
 b is a neighbor.  The edgeType can be either 'tree' (EDGE_DFSPARENT or
 EDGE_DFSCHILD) or 'cycle' (EDGE_BACK or EDGE_FORWARD).
 If a or b is a root copy, we translate to the non-virtual counterpart,
 then wedetermine which has the lesser DFI.  If a has the lower DFI then the
 edge record is a tree edge to a child (EDGE_DFSCHILD) if edgeType indicates
 a tree edge.  If edgeType indicates a cycle edge, then it is a forward cycle
 edge (EDGE_FORWARD) to a descendant.
 Symmetric conditions define the types for a > b.
 ****************************************************************************/

int  _ComputeArcType(graphP theGraph, int a, int b, int edgeType)
{
     if (a >= theGraph->N)
         a = theGraph->V[a - theGraph->N].DFSParent;
     if (b >= theGraph->N)
         b = theGraph->V[b - theGraph->N].DFSParent;

     if (a < b)
         return edgeType == EDGE_DFSPARENT || edgeType == EDGE_DFSCHILD ? EDGE_DFSCHILD : EDGE_FORWARD;

     return edgeType == EDGE_DFSPARENT || edgeType == EDGE_DFSCHILD ? EDGE_DFSPARENT : EDGE_BACK;
}

/****************************************************************************
 _SetEdgeType()
 When we are restoring an edge, we must restore its type (tree edge or cycle edge).
 We can deduce what the type was based on other information in the graph. Each
 arc of the edge gets the appropriate type setting (parent/child or back/forward).
 This method runs in constant time plus the degree of vertex u, or constant
 time if u is known to have a degree bound by a constant.
 ****************************************************************************/

int  _SetEdgeType(graphP theGraph, int u, int v)
{
int  e, eTwin, u_orig, v_orig, N;

     // If u or v is a virtual vertex (a root copy), then get the non-virtual counterpart.
     N = theGraph->N;
     u_orig = u < N ? u : (theGraph->V[u - N].DFSParent);
     v_orig = v < N ? v : (theGraph->V[v - N].DFSParent);

     // Get the edge for which we will set the type

     e = gp_GetNeighborEdgeRecord(theGraph, u, v);
     eTwin = gp_GetTwinArc(theGraph, e);

     // If u_orig is the parent of v_orig, or vice versa, then the edge is a tree edge

     if (theGraph->V[v_orig].DFSParent == u_orig ||
         theGraph->V[u_orig].DFSParent == v_orig)
     {
         if (u_orig > v_orig)
         {
             theGraph->G[e].type = EDGE_DFSPARENT;
             theGraph->G[eTwin].type = EDGE_DFSCHILD;
         }
         else
         {
             theGraph->G[eTwin].type = EDGE_DFSPARENT;
             theGraph->G[e].type = EDGE_DFSCHILD;
         }
     }

     // Otherwise it is a back edge

     else
     {
         if (u_orig > v_orig)
         {
             theGraph->G[e].type = EDGE_BACK;
             theGraph->G[eTwin].type = EDGE_FORWARD;
         }
         else
         {
             theGraph->G[eTwin].type = EDGE_BACK;
             theGraph->G[e].type = EDGE_FORWARD;
         }
     }

     return OK;
}

/********************************************************************
 _DeleteUnmarkedEdgesInBicomp()

 This function deletes from a given biconnected component all edges
 whose visited member is zero.

 The stack is used but preserved. In debug mode, NOTOK can result if
 there is a stack overflow. This method pushes at most one integer
 per vertex in the bicomp.

 Returns OK on success, NOTOK on implementation failure
 ********************************************************************/

int  _DeleteUnmarkedEdgesInBicomp(graphP theGraph, int BicompRoot)
{
int  V, J;
int  stackBottom = sp_GetCurrentSize(theGraph->theStack);

     sp_Push(theGraph->theStack, BicompRoot);
     while (sp_GetCurrentSize(theGraph->theStack) > stackBottom)
     {
          sp_Pop(theGraph->theStack, V);

          J = gp_GetFirstArc(theGraph, V);
          while (gp_IsArc(theGraph, J))
          {
             if (theGraph->G[J].type == EDGE_DFSCHILD)
                 sp_Push(theGraph->theStack, theGraph->G[J].v);

             if (!theGraph->G[J].visited)
                  J = gp_DeleteEdge(theGraph, J, 0);
             else J = gp_GetNextArc(theGraph, J);
          }
     }
     return OK;
}

/********************************************************************
 _ClearInvertedFlagsInBicomp()

 This function clears the inverted flag markers on any edges in a
 given biconnected component.

 The stack is used but preserved. In debug mode, NOTOK can result if
 there is a stack overflow. This method pushes at most one integer
 per vertex in the bicomp.

 Returns OK on success, NOTOK on implementation failure
 ********************************************************************/

int  _ClearInvertedFlagsInBicomp(graphP theGraph, int BicompRoot)
{
int  V, J;
int  stackBottom = sp_GetCurrentSize(theGraph->theStack);

     sp_Push(theGraph->theStack, BicompRoot);
     while (sp_GetCurrentSize(theGraph->theStack) > stackBottom)
     {
          sp_Pop(theGraph->theStack, V);

          J = gp_GetFirstArc(theGraph, V);
          while (gp_IsArc(theGraph, J))
          {
             if (theGraph->G[J].type == EDGE_DFSCHILD)
             {
                 sp_Push(theGraph->theStack, theGraph->G[J].v);
                 CLEAR_EDGEFLAG_INVERTED(theGraph, J);
             }

             J = gp_GetNextArc(theGraph, J);
          }
     }
     return OK;
}

/********************************************************************
 _GetBicompSize()

 Determine the number of vertices in the bicomp.

 The stack is used but preserved. In debug mode, NOTOK can result if
 there is a stack overflow. This method pushes at most one integer
 per vertex in the bicomp.

 Returns a positive number on success, NOTOK on implementation failure
 ********************************************************************/

int  _GetBicompSize(graphP theGraph, int BicompRoot)
{
int  V, J;
int  theSize = 0;
int  stackBottom = sp_GetCurrentSize(theGraph->theStack);

     sp_Push(theGraph->theStack, BicompRoot);
     while (sp_GetCurrentSize(theGraph->theStack) > stackBottom)
     {
          sp_Pop(theGraph->theStack, V);
          theSize++;
          J = gp_GetFirstArc(theGraph, V);
          while (gp_IsArc(theGraph, J))
          {
             if (theGraph->G[J].type == EDGE_DFSCHILD)
                 sp_Push(theGraph->theStack, theGraph->G[J].v);

             J = gp_GetNextArc(theGraph, J);
          }
     }
     return theSize;
}

/********************************************************************
 debugNOTOK()
 This function provides a non-void wrapper for exit().
 This is useful for debugging as it allows compilation of an exit
 command in places where NOTOK is returned.
 In exhaustive testing, we want to bail on the first NOTOK that occurs.
 Comment out the exit() call to get a stack trace.
 ********************************************************************/

int debugNOTOK(void)
{
	//exit(-1);
	return 0; // NOTOK is normally defined to be zero
}
