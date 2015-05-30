#ifndef GRAPHFUNCTIONTABLE_H
#define GRAPHFUNCTIONTABLE_H

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

#ifdef __cplusplus
extern "C" {
#endif

/*
 NOTE: If you add any FUNCTION POINTERS to this function table, then you must
       also initialize them in _InitFunctionTable() in graphUtils.c.
*/

typedef struct
{
        // These function pointers allow extension modules to overload some of
        // the behaviors of protected functions.  Only advanced applications
        // will overload these functions
        int  (*fpCreateFwdArcLists)(void *theGraph);
        void (*fpCreateDFSTreeEmbedding)(void *theGraph);
        void (*fpEmbedBackEdgeToDescendant)(void *theGraph, int rootSide, int rootVertex, int W, int WPrevLink);
        void (*fpWalkUp)(void *theGraph, int I, int J);
        int  (*fpWalkDown)(void *theGraph, int I, int rootVertex);
        int  (*fpMergeBicomps)(void *theGraph, int I, int rootVertex, int W, int WPrevLink);
        int  (*fpHandleInactiveVertex)(void *theGraph, int bicompRoot, int* pW, int *pWPrevLink);
        int  (*fpHandleBlockedDescendantBicomp)(void *theGraph, int I, int rootVertex, int R, int *pRout, int *pW, int *pWPrevLink);
        int  (*fpHandleBlockedEmbedIteration)(void *theGraph, int I);
        int  (*fpEmbedPostprocess)(void *theGraph, int I, int edgeEmbeddingResult);
        int  (*fpMarkDFSPath)(void *theGraph, int ancestor, int descendant);

        int  (*fpCheckEmbeddingIntegrity)(void *theGraph, void *origGraph);
        int  (*fpCheckObstructionIntegrity)(void *theGraph, void *origGraph);

        // These function pointers allow extension modules to overload
        // vertex and graphnode initialization. These are not part of the
        // public API, but many extensions are expected to overload them
        // if they equip vertices or edges with additional parameters
        void (*fpInitGraphNode)(void *theGraph, int I);
        void (*fpInitVertexRec)(void *theGraph, int I);

        // These function pointers allow extension modules to overload some
        // of the behaviors of gp_* function in the public API
        int  (*fpInitGraph)(void *theGraph, int size);
        void (*fpReinitializeGraph)(void *theGraph);
        int  (*fpEnsureArcCapacity)(void *theGraph, int requiredArcCapacity);
        int  (*fpSortVertices)(void *theGraph);

        int  (*fpReadPostprocess)(void *theGraph, void *extraData, long extraDataSize);
        int  (*fpWritePostprocess)(void *theGraph, void **pExtraData, long *pExtraDataSize);

        int  (*fpHideVertex)(void *theGraph, int vertex);
        void (*fpHideEdge)(void *theGraph, int e);
        void (*fpRestoreEdge)(void *theGraph, int arcPos);
        int  (*fpContractEdge)(void *theGraph, int e);
        int  (*fpIdentifyVertices)(void *theGraph, int u, int v, int eBefore);
        int  (*fpRestoreVertex)(void *theGraph);

} graphFunctionTable;

typedef graphFunctionTable * graphFunctionTableP;

#ifdef __cplusplus
}
#endif

#endif
