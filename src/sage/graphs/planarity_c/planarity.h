#ifndef PLANARITY_H
#define PLANARITY_H

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

#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include "graph.h"
#include "platformTime.h"

#include "graphK23Search.h"
#include "graphK33Search.h"
#include "graphK4Search.h"
#include "graphDrawPlanar.h"
#include "graphColorVertices.h"

void ProjectTitle(void);
int helpMessage(char *param);

/* Functions that call the Graph Library */
int SpecificGraph(char command, char *infileName, char *outfileName, char *outfile2Name);
int RandomGraph(char command, int extraEdges, int numVertices, char *outfileName, char *outfile2Name);
int RandomGraphs(char command, int, int);

int makeg_main(char command, int argc, char *argv[]);

/* Command line, Menu, and Configuration */
int commandLine(int argc, char *argv[]);
int legacyCommandLine(int argc, char *argv[]);
int menu(void);

char Mode,
     OrigOut,
     EmbeddableOut,
     ObstructedOut,
     AdjListsForEmbeddingsOut,
     quietMode;

void Reconfigure(void);

/* Low-level Utilities */
#define MAXLINE 1024
char Line[MAXLINE];

void Message(char *message);
void ErrorMessage(char *message);
void FlushConsole(FILE *f);
void Prompt(char *message);

void SaveAsciiGraph(graphP theGraph, char *filename);

int  FilesEqual(char *file1Name, char *file2Name);

int GetEmbedFlags(char command);
char *GetAlgorithmName(char command);
void AttachAlgorithm(graphP theGraph, char command);

char *ConstructInputFilename(char *infileName);
char *ConstructPrimaryOutputFilename(char *infileName, char *outfileName, char command);
void WriteAlgorithmResults(graphP theGraph, int Result, char command, platform_time start, platform_time end, char *infileName);

#ifdef __cplusplus
}
#endif

#endif
