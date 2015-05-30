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

void ProjectTitle(void)
{
    Message("\n=================================================="
            "\nPlanarity version 2.2"
            "\nCopyright (c) 2010 by John M. Boyer"
    		"\nContact info: jboyer at acm.org"
            "\n=================================================="
            "\n");
}

/****************************************************************************
 MAIN
 ****************************************************************************/

int main(int argc, char *argv[])
{
	int retVal=0;

	if (argc <= 1)
		retVal = menu();

	else if (argv[1][0] == '-')
		retVal = commandLine(argc, argv);

	else
		retVal = legacyCommandLine(argc, argv);

	// Close the log file if logging
	gp_Log(NULL);

	return retVal;
}

/****************************************************************************
 helpMessage()
 ****************************************************************************/

int helpMessage(char *param)
{
	char *commandStr =
    	"C = command from menu\n"
    	"    -p = Planar embedding and Kuratowski subgraph isolation\n"
        "    -o = Outerplanar embedding and obstruction isolation\n"
        "    -d = Planar graph drawing\n"
        "    -2 = Search for subgraph homeomorphic to K_{2,3}\n"
        "    -3 = Search for subgraph homeomorphic to K_{3,3}\n"
        "    -4 = Search for subgraph homeomorphic to K_4\n"
		"    -c = Color the vertices of the graph\n"
    	"\n";

	ProjectTitle();

	if (param == NULL)
	{
	    Message(
            "'planarity': menu-driven\n"
            "'planarity (-h|-help)': this message\n"
            "'planarity (-h|-help) -menu': more help with menu-based command line\n"
    	    "'planarity -test [-q] [C]': runs tests (optional quiet mode, single test)\n"
	    	"\n"
	    );

	    Message(
	    	"Common usages\n"
	    	"-------------\n"
            "planarity -s -q -p infile.txt embedding.out [obstruction.out]\n"
	    	"Process infile.txt in quiet mode (-q), putting planar embedding in \n"
	    	"embedding.out or (optionally) a Kuratowski subgraph in Obstruction.out\n"
	    	"Process returns 0=planar, 1=nonplanar, -1=error\n"
	    	"\n"
            "planarity -s -q -d infile.txt embedding.out [drawing.out]\n"
            "If graph in infile.txt is planar, then put embedding in embedding.out \n"
            "and (optionally) an ASCII art drawing in drawing.out\n"
            "Process returns 0=planar, 1=nonplanar, -1=error\n"
	    	"\n"
	    );
	}

	else if (strcmp(param, "-menu") == 0)
	{
	    Message(
	    	"'planarity -r [-q] C K N': Random graphs\n"
	    	"'planarity -s [-q] C I O [O2]': Specific graph\n"
	        "'planarity -rm [-q] N O [O2]': Maximal planar random graph\n"
	        "'planarity -rn [-q] N O [O2]': Nonplanar random graph (maximal planar + edge)\n"
	        "'planarity I O [-n O2]': Legacy command-line (default -s -p)\n"
	    	"\n"
	    );

	    Message("-q is for quiet mode (no messages to stdout and stderr)\n\n");

	    Message(commandStr);

	    Message(
	    	"K = # of graphs to randomly generate\n"
	    	"N = # of vertices in each randomly generated graph\n"
	        "I = Input file (for work on a specific graph)\n"
	        "O = Primary output file\n"
	        "    For example, if C=-p then O receives the planar embedding\n"
	    	"    If C=-3, then O receives a subgraph containing a K_{3,3}\n"
	        "O2= Secondary output file\n"
	    	"    For -s, if C=-p or -o, then O2 receives the embedding obstruction\n"
	       	"    For -s, if C=-d, then O2 receives a drawing of the planar graph\n"
	    	"    For -m and -n, O2 contains the original randomly generated graph\n"
	    	"\n"
	    );

	    Message(
	        "planarity process results: 0=OK, -1=NOTOK, 1=NONEMBEDDABLE\n"
	    	"    1 result only produced by specific graph mode (-s)\n"
	        "      with command -2,-3,-4: found K_{2,3}, K_{3,3} or K_4\n"
	    	"      with command -p,-d: found planarity obstruction\n"
	    	"      with command -o: found outerplanarity obstruction\n"
	    );
	}

    FlushConsole(stdout);
    return 0;
}

/****************************************************************************
 MENU-DRIVEN PROGRAM
 ****************************************************************************/

int menu(void)
{
char Choice;

     do {
    	ProjectTitle();

        Message("\n"
                "P. Planar embedding and Kuratowski subgraph isolation\n"
                "D. Planar graph drawing\n"
                "O. Outerplanar embedding and obstruction isolation\n"
                "2. Search for subgraph homeomorphic to K_{2,3}\n"
                "3. Search for subgraph homeomorphic to K_{3,3}\n"
                "4. Search for subgraph homeomorphic to K_4\n"
        		"C. Color the vertices of the graph\n"
        		"H. Help message for command line version\n"
                "R. Reconfigure options\n"
                "X. Exit\n"
        		"\n"
        );

        Prompt("Enter Choice: ");
        fflush(stdin);
        scanf(" %c", &Choice);
        Choice = tolower(Choice);

        if (Choice == 'h')
        	helpMessage(NULL);

        else if (Choice == 'r')
        	Reconfigure();

        else if (Choice != 'x')
        {
        	char *secondOutfile = NULL;
        	if (Choice == 'p'  || Choice == 'o' || Choice == 'd')
        		secondOutfile ="";

            switch (tolower(Mode))
            {
                case 's' : SpecificGraph(Choice, NULL, NULL, secondOutfile); break;
                case 'r' : RandomGraphs(Choice, 0, 0); break;
                case 'm' : RandomGraph(Choice, 0, 0, NULL, NULL); break;
                case 'n' : RandomGraph(Choice, 1, 0, NULL, NULL); break;
            }
        }

        if (Choice != 'r' && Choice != 'x')
        {
            Prompt("\nPress a key then hit ENTER to continue...");
            fflush(stdin);
            scanf(" %*c");
            fflush(stdin);
            Message("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
            FlushConsole(stdout);
        }

     }  while (Choice != 'x');

     // Certain debuggers don't terminate correctly with pending output content
     FlushConsole(stdout);
     FlushConsole(stderr);

     return 0;
}
