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

#include "planarity.h"

/****************************************************************************
 Configuration
 ****************************************************************************/

char Mode='r',
     OrigOut='n',
     EmbeddableOut='n',
     ObstructedOut='n',
     AdjListsForEmbeddingsOut='n',
     quietMode='n';

void Reconfigure(void)
{
     fflush(stdin);

     Prompt("\nDo you want to \n"
    		"  Randomly generate graphs (r),\n"
    		"  Specify a graph (s),\n"
    		"  Randomly generate a maximal planar graph (m), or\n"
    		"  Randomly generate a non-planar graph (n)?");
     scanf(" %c", &Mode);

     Mode = tolower(Mode);
     if (!strchr("rsmn", Mode))
    	 Mode = 's';

     if (Mode == 'r')
     {
        Message("\nNOTE: The directories for the graphs you want must exist.\n\n");

        Prompt("Do you want original graphs in directory 'random' (last 10 max)?");
        scanf(" %c", &OrigOut);

        Prompt("Do you want adj. matrix of embeddable graphs in directory 'embedded' (last 10 max))?");
        scanf(" %c", &EmbeddableOut);

        Prompt("Do you want adj. matrix of obstructed graphs in directory 'obstructed' (last 10 max)?");
        scanf(" %c", &ObstructedOut);

        Prompt("Do you want adjacency list format of embeddings in directory 'adjlist' (last 10 max)?");
        scanf(" %c", &AdjListsForEmbeddingsOut);
     }

     FlushConsole(stdout);
}

/****************************************************************************
 MESSAGE - prints a string, but when debugging adds \n and flushes stdout
 ****************************************************************************/

#define MAXLINE 1024
char Line[MAXLINE];

void Message(char *message)
{
	if (quietMode == 'n')
	{
	    fprintf(stdout, "%s", message);

#ifdef DEBUG
//	    fprintf(stdout, "\n");
	    fflush(stdout);
#endif
	}
}

void ErrorMessage(char *message)
{
	if (quietMode == 'n')
	{
		fprintf(stderr, "%s", message);

#ifdef DEBUG
		fprintf(stderr, "\n");
		fflush(stderr);
#endif
	}
}

void FlushConsole(FILE *f)
{
#ifdef DEBUG
	    // Certain debuggers only flush completed lines of output to the console
	    fprintf(f, "\n");
#endif
	    fflush(f);
}

void Prompt(char *message)
{
	Message(message);
	FlushConsole(stdout);
}

/****************************************************************************
 ****************************************************************************/

void SaveAsciiGraph(graphP theGraph, char *filename)
{
	int  e, limit;
	FILE *outfile = fopen(filename, "wt");
	fprintf(outfile, "%s\n", filename);

	limit = theGraph->edgeOffset + 2*(theGraph->M + sp_GetCurrentSize(theGraph->edgeHoles));

	for (e = theGraph->edgeOffset; e < limit; e+=2)
	{
		if (theGraph->G[e].v != NIL)
			fprintf(outfile, "%d %d\n", theGraph->G[e].v+1, theGraph->G[e+1].v+1);
	}

	fprintf(outfile, "0 0\n");

	fclose(outfile);
}

/****************************************************************************
 ****************************************************************************/

int  FilesEqual(char *file1Name, char *file2Name)
{
	FILE *infile1 = NULL, *infile2 = NULL;
	int Result = TRUE;

	infile1 = fopen(file1Name, "r");
	infile2 = fopen(file2Name, "r");

	if (infile1 == NULL || infile2 == NULL)
		Result = FALSE;
	else
	{
		int c1=0, c2=0;

		// Read the first file to the end
		while ((c1 = fgetc(infile1)) != EOF)
		{
			// If we got a char from the first file, but not from the second
			// then the second file is shorter, so files are not equal
			if ((c2 = fgetc(infile2)) == EOF)
			{
				Result = FALSE;
				break;
			}

			// If we got a char from second file, but not equal to char from
			// first file, then files are not equal
			if (c1 != c2)
			{
				Result = FALSE;
				break;
			}
		}

		// If we got to the end of the first file without breaking the loop...
		if (c1 == EOF)
		{
			// Then attempt to read from the second file to ensure it also ends.
			if (fgetc(infile2) != EOF)
				Result = FALSE;
		}
	}

	if (infile1 != NULL) fclose(infile1);
	if (infile2 != NULL) fclose(infile2);
	return Result;
}

/****************************************************************************
 ****************************************************************************/

int GetEmbedFlags(char command)
{
	int embedFlags = 0;

	switch (command)
	{
		case 'o' : embedFlags = EMBEDFLAGS_OUTERPLANAR; break;
		case 'p' : embedFlags = EMBEDFLAGS_PLANAR; break;
		case 'd' : embedFlags = EMBEDFLAGS_DRAWPLANAR; break;
		case '2' : embedFlags = EMBEDFLAGS_SEARCHFORK23; break;
		case '3' : embedFlags = EMBEDFLAGS_SEARCHFORK33; break;
		case '4' : embedFlags = EMBEDFLAGS_SEARCHFORK4; break;
	}

	return embedFlags;
}

/****************************************************************************
 ****************************************************************************/

char *GetAlgorithmName(char command)
{
	char *algorithmName = "UnsupportedAlgorithm";

	switch (command)
	{
		case 'p' : algorithmName = "PlanarEmbed"; break;
		case 'd' : algorithmName = DRAWPLANAR_NAME;	break;
		case 'o' : algorithmName = "OuterplanarEmbed"; break;
		case '2' : algorithmName = K23SEARCH_NAME; break;
		case '3' : algorithmName = K33SEARCH_NAME; break;
		case '4' : algorithmName = K4SEARCH_NAME; break;
		case 'c' : algorithmName = COLORVERTICES_NAME; break;
	}

	return algorithmName;
}

/****************************************************************************
 ****************************************************************************/

void AttachAlgorithm(graphP theGraph, char command)
{
	switch (command)
	{
		case 'd' : gp_AttachDrawPlanar(theGraph); break;
		case '2' : gp_AttachK23Search(theGraph); break;
		case '3' : gp_AttachK33Search(theGraph); break;
		case '4' : gp_AttachK4Search(theGraph); break;
		case 'c' : gp_AttachColorVertices(theGraph); break;
	}
}

/****************************************************************************
 A string used to construct input and output filenames.

 The SUFFIXMAXLENGTH is 32 to accommodate ".out.txt" + ".render.txt" + ".test.txt"
 ****************************************************************************/

#define FILENAMEMAXLENGTH 128
#define ALGORITHMNAMEMAXLENGTH 32
#define SUFFIXMAXLENGTH 32

char theFileName[FILENAMEMAXLENGTH+1+ALGORITHMNAMEMAXLENGTH+1+SUFFIXMAXLENGTH+1];

/****************************************************************************
 ConstructInputFilename()
 Returns a string not owned by the caller (do not free string).
 String contains infileName content if infileName is non-NULL.
 If infileName is NULL, then the user is asked to supply a name.
 Returns NULL on error, or a non-NULL string on success.
 ****************************************************************************/

char *ConstructInputFilename(char *infileName)
{
	if (infileName == NULL)
	{
		Prompt("Enter graph file name: ");
		fflush(stdin);
		scanf(" %s", theFileName);

		if (!strchr(theFileName, '.'))
			strcat(theFileName, ".txt");
	}
	else
	{
		if (strlen(infileName) > FILENAMEMAXLENGTH)
		{
			ErrorMessage("Filename is too long");
			return NULL;
		}
		strcpy(theFileName, infileName);
	}

	return theFileName;
}

/****************************************************************************
 ConstructPrimaryOutputFilename()
 Returns a string not owned by the caller (do not free string).
 Reuses the same memory space as ConstructinputFilename().
 If outfileName is non-NULL, then the result string contains its content.
 If outfileName is NULL, then the infileName and the command's algorithm name
 are used to construct a string.
 Returns non-NULL string
 ****************************************************************************/

char *ConstructPrimaryOutputFilename(char *infileName, char *outfileName, char command)
{
	char *algorithmName = GetAlgorithmName(command);

	if (outfileName == NULL)
	{
		// The output filename is based on the input filename
		if (theFileName != infileName)
		    strcpy(theFileName, infileName);

		// If the primary output filename has not been given, then we use
		// the input filename + the algorithm name + a simple suffix
		if (strlen(algorithmName) <= ALGORITHMNAMEMAXLENGTH)
		{
			strcat(theFileName, ".");
			strcat(theFileName, algorithmName);
		}
		else
			ErrorMessage("Algorithm Name is too long, so it will not be used in output filename.");

	    strcat(theFileName, ".out.txt");
	}
	else
	{
		if (strlen(outfileName) > FILENAMEMAXLENGTH)
		{
			// The output filename is based on the input filename
			if (theFileName != infileName)
			    strcpy(theFileName, infileName);

	    	if (strlen(algorithmName) <= ALGORITHMNAMEMAXLENGTH)
	    	{
	    		strcat(theFileName, ".");
	    		strcat(theFileName, algorithmName);
	    	}
	        strcat(theFileName, ".out.txt");
			sprintf(Line, "Outfile filename is too long. Result placed in %s", theFileName);
			ErrorMessage(Line);
		}
		else
		{
			if (theFileName != outfileName)
			    strcpy(theFileName, outfileName);
		}
	}

	return theFileName;
}
