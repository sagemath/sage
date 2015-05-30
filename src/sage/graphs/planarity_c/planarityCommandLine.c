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

#include <unistd.h>

int callRandomGraphs(int argc, char *argv[]);
int callSpecificGraph(int argc, char *argv[]);
int callRandomMaxPlanarGraph(int argc, char *argv[]);
int callRandomNonplanarGraph(int argc, char *argv[]);

/****************************************************************************
 Command Line Processor
 ****************************************************************************/

int commandLine(int argc, char *argv[])
{
	int Result = OK;

	if (argc >= 3 && strcmp(argv[2], "-q") == 0)
		quietMode = 'y';

	if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-help") == 0)
	{
		Result = helpMessage(argc >= 3 ? argv[2] : NULL);
	}

	else if (strcmp(argv[1], "-r") == 0)
		Result = callRandomGraphs(argc, argv);

	else if (strcmp(argv[1], "-s") == 0)
		Result = callSpecificGraph(argc, argv);

	else if (strcmp(argv[1], "-rm") == 0)
		Result = callRandomMaxPlanarGraph(argc, argv);

	else if (strcmp(argv[1], "-rn") == 0)
		Result = callRandomNonplanarGraph(argc, argv);

	else
	{
		ErrorMessage("Unsupported command line.  Here is the help for this program.\n");
		helpMessage(NULL);
		Result = NOTOK;
	}

	return Result == OK ? 0 : (Result == NONEMBEDDABLE ? 1 : -1);
}

/****************************************************************************
 Legacy Command Line Processor from version 1.x
 ****************************************************************************/

int legacyCommandLine(int argc, char *argv[])
{
graphP theGraph = gp_New();
int Result;

	Result = gp_Read(theGraph, argv[1]);
	if (Result != OK)
	{
		if (Result != NONEMBEDDABLE)
		{
			if (strlen(argv[1]) > MAXLINE - 100)
				sprintf(Line, "Failed to read graph\n");
			else
				sprintf(Line, "Failed to read graph %s\n", argv[1]);
			ErrorMessage(Line);
			return -2;
		}
	}

	Result = gp_Embed(theGraph, EMBEDFLAGS_PLANAR);

	if (Result == OK)
	{
		gp_SortVertices(theGraph);
		gp_Write(theGraph, argv[2], WRITE_ADJLIST);
	}

	else if (Result == NONEMBEDDABLE)
	{
		if (argc >= 5 && strcmp(argv[3], "-n")==0)
		{
			gp_SortVertices(theGraph);
			gp_Write(theGraph, argv[4], WRITE_ADJLIST);
		}
	}
	else
		Result = NOTOK;

	gp_Free(&theGraph);

	// In the legacy 1.x versions, OK/NONEMBEDDABLE was 0 and NOTOK was -2
	return Result==OK || Result==NONEMBEDDABLE ? 0 : -2;
}


/****************************************************************************
 Quick regression test
 ****************************************************************************/

int runSpecificGraphTests(void);
int runSpecificGraphTest(char *command, char *infileName);

int runSpecificGraphTests(void)
{
	char origDir[2049];
	int retVal = 0;

	if (!getcwd(origDir, 2048))
		return -1;

	if (chdir("samples") != 0)
	{
		if (chdir("..") != 0 || chdir("samples") != 0)
		{
			// Warn but give success result
			printf("WARNING: Unable to change to samples directory to run tests on samples.\n");
			return 0;
		}
	}

	if (runSpecificGraphTest("-p", "maxPlanar5.txt") < 0)
		retVal = -1;

	if (runSpecificGraphTest("-d", "maxPlanar5.txt") < 0)
		retVal = -1;

	if (runSpecificGraphTest("-d", "drawExample.txt") < 0)
		retVal = -1;

	if (runSpecificGraphTest("-p", "Petersen.txt") < 0)
		retVal = -1;

	if (runSpecificGraphTest("-o", "Petersen.txt") < 0)
		retVal = -1;

	if (runSpecificGraphTest("-2", "Petersen.txt") < 0)
		retVal = -1;

	if (runSpecificGraphTest("-3", "Petersen.txt") < 0)
		retVal = -1;

	if (runSpecificGraphTest("-4", "Petersen.txt") < 0)
		retVal = -1;

	if (runSpecificGraphTest("-c", "maxPlanar5.txt") < 0)
		retVal = -1;

	if (runSpecificGraphTest("-c", "Petersen.txt") < 0)
		retVal = -1;

	if (runSpecificGraphTest("-c", "drawExample.txt") < 0)
		retVal = -1;

	chdir(origDir);
    FlushConsole(stdout);
	return retVal;
}

int runSpecificGraphTest(char *command, char *infileName)
{
	char *commandLine[] = {
			"planarity", "-s", "C", "infile", "outfile", "outfile2"
	};
	char *outfileName = ConstructPrimaryOutputFilename(infileName, NULL, command[1]);
	char *outfile2Name = "";
	char *testfileName = strdup(outfileName);
	int Result = 0;

	if (testfileName == NULL)
		return -1;

	outfileName = strdup(strcat(outfileName, ".test.txt"));
	if (outfileName == NULL)
	{
		free(testfileName);
		return -1;
	}

	// 'planarity -s [-q] C I O [O2]': Specific graph
	commandLine[2] = command;
	commandLine[3] = infileName;
	commandLine[4] = outfileName;
	commandLine[5] = outfile2Name;

	Result = callSpecificGraph(6, commandLine);
	if (Result == OK || Result == NONEMBEDDABLE)
		Result = 0;
	else
	{
		ErrorMessage("Test failed (graph processor returned failure result).\n");
		Result = -1;
	}

	if (Result == 0)
	{
		if (FilesEqual(testfileName, outfileName) == TRUE)
		{
			Message("Test succeeded (result equal to exemplar).\n");
			unlink(outfileName);
		}
		else
		{
			ErrorMessage("Test failed (result not equal to exemplar).\n");
			Result = -1;
		}
	}

	// For graph drawing, secondary file is outfileName + ".render.txt"

	if (command[1] == 'd' && Result == 0)
	{
		outfile2Name = ConstructPrimaryOutputFilename(NULL, outfileName, command[1]);
		free(outfileName);
		outfileName = strdup(strcat(outfile2Name, ".render.txt"));

		free(testfileName);
		testfileName = ConstructPrimaryOutputFilename(infileName, NULL, command[1]);
		testfileName = strdup(strcat(testfileName, ".render.txt"));

		if (Result == 0)
		{
			if (FilesEqual(testfileName, outfileName) == TRUE)
			{
				Message("Test succeeded (secondary result equal to exemplar).\n");
				unlink(outfileName);
			}
			else
			{
				ErrorMessage("Test failed (secondary result not equal to exemplar).\n");
				Result = -1;
			}
		}
	}

	Message("\n");

	free(outfileName);
	free(testfileName);
	return Result;
}


/****************************************************************************
 callRandomGraphs()
 ****************************************************************************/

// 'planarity -r [-q] C K N': Random graphs
int callRandomGraphs(int argc, char *argv[])
{
	char Choice = 0;
	int offset = 0, NumGraphs, SizeOfGraphs;

	if (argc < 5)
		return -1;

	if (argv[2][0] == '-' && (Choice = argv[2][1]) == 'q')
	{
		Choice = argv[3][1];
		if (argc < 6)
			return -1;
		offset = 1;
	}

	NumGraphs = atoi(argv[3+offset]);
	SizeOfGraphs = atoi(argv[4+offset]);

    return RandomGraphs(Choice, NumGraphs, SizeOfGraphs);
}

/****************************************************************************
 callSpecificGraph()
 ****************************************************************************/

// 'planarity -s [-q] C I O [O2]': Specific graph
int callSpecificGraph(int argc, char *argv[])
{
	char Choice=0, *infileName=NULL, *outfileName=NULL, *outfile2Name=NULL;
	int offset = 0;

	if (argc < 5)
		return -1;

	if (argv[2][0] == '-' && (Choice = argv[2][1]) == 'q')
	{
		Choice = argv[3][1];
		if (argc < 6)
			return -1;
		offset = 1;
	}

	infileName = argv[3+offset];
	outfileName = argv[4+offset];
	if (argc == 6+offset)
	    outfile2Name = argv[5+offset];

	return SpecificGraph(Choice, infileName, outfileName, outfile2Name);
}

/****************************************************************************
 callRandomMaxPlanarGraph()
 ****************************************************************************/

// 'planarity -rm [-q] N O [O2]': Maximal planar random graph
int callRandomMaxPlanarGraph(int argc, char *argv[])
{
	int offset = 0, numVertices;
	char *outfileName = NULL, *outfile2Name = NULL;

	if (argc < 4)
		return -1;

	if (argv[2][0] == '-' && argv[2][1] == 'q')
	{
		if (argc < 5)
			return -1;
		offset = 1;
	}

	numVertices = atoi(argv[2+offset]);
	outfileName = argv[3+offset];
	if (argc == 5+offset)
	    outfile2Name = argv[4+offset];

	return RandomGraph('p', 0, numVertices, outfileName, outfile2Name);
}

/****************************************************************************
 callRandomNonplanarGraph()
 ****************************************************************************/

// 'planarity -rn [-q] N O [O2]': Non-planar random graph (maximal planar plus edge)
int callRandomNonplanarGraph(int argc, char *argv[])
{
	int offset = 0, numVertices;
	char *outfileName = NULL, *outfile2Name = NULL;

	if (argc < 4)
		return -1;

	if (argv[2][0] == '-' && argv[2][1] == 'q')
	{
		if (argc < 5)
			return -1;
		offset = 1;
	}

	numVertices = atoi(argv[2+offset]);
	outfileName = argv[3+offset];
	if (argc == 5+offset)
	    outfile2Name = argv[4+offset];

	return RandomGraph('p', 1, numVertices, outfileName, outfile2Name);
}
