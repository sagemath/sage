# -*- coding: utf-8 -*-
r"""
Convert output from code cells in the notebook into ReStructuredText

This is called by sws2rst

- Pablo Angulo Ardoy (2011-02-25): initial version
"""
#**************************************************
# Copyright (C) 2011 Pablo Angulo
#
# Distributed under the terms of the GPL License
#**************************************************


import re
IMAGES_DIR = 'images/'

#We parse lines one by one but keep track of current scope
#similarly to worksheet2rst.py
#Results are split into different types. Some are discarded
class States(object):
    NORMAL = 0
    HTML = 1
    MATH = 2
    TRACEBACK = 3

class LineTypes(object):
    PLAIN = 0
    IMAGE = 1
    LATEX = 2
    HTML  = 3
    TRACE = 4

class ResultsParser(object):
    """Auxiliary class for results2rst
    """
    def __init__(self, images_dir):
        ##Order matters, place more restrictive regex's before more general ones
        ##If no regex matches, line will be discarded
        ##a self transition is needes to produce any output
        self.transitions = {
            States.NORMAL:[
                #IMAGE
                     (re.compile(r"^\<html\>\<font color='black'\>"
                                 r"\<img src='cell\://(.*?)'\>"
                                 r"\</font\>\</html\>"),
                      "\n.. image:: " + images_dir + "\\1\n    :align: center\n",
                      LineTypes.IMAGE,
                      States.NORMAL),
                #SELF-CONTAINED MATH
                     (re.compile(r"^\<html\>\<div class=\"math\"\>"
                                 r"\\newcommand\{\\Bold\}\[1\]\{\\mathbf\{\#1\}\}"
                                 r"(.*?)\</div\>\</html\>$"),
                      "\n.. MATH::\n\n    \\1\n",
                      LineTypes.LATEX,
                      States.NORMAL),
                #SELF-CONTAINED MATH - BIS
                     (re.compile(r"^\<html\>\<div class=\"math\"\>"
                                 r"(.*?)\</div\>\</html\>$"),
                      "\n.. MATH::\n\n    \\1",
                      LineTypes.LATEX,
                      States.NORMAL),
                #START Traceback
                     (re.compile(r"^(Traceback.*)"),
                      "    Traceback (most recent call last):",
                      LineTypes.TRACE,
                      States.TRACEBACK),
                #START MATH
                     (re.compile(r"^\<html\>\<div class=\"math\"\>"
                                 r"\\newcommand\{\\Bold\}\[1\]\{\\mathbf\{\#1\}\}(.*?)"),
                      "\n.. MATH::\n\n    \\1",
                      LineTypes.LATEX,
                      States.MATH),
                #SELF-CONTAINED HTML
                     (re.compile(r"^\<html\>.*</html\>$"),
                      "    <html>...</html>",
                      LineTypes.HTML,
                      States.NORMAL),        
                #START HTML
                     (re.compile(r"^\<html\>.*"),
                      "    <html>...</html>",
                      LineTypes.HTML,
                      States.HTML),        
                #CONTINUE NORMAL
                     (re.compile("(.*)"),
                      "    \\1",
                      LineTypes.PLAIN,
                      States.NORMAL),                
                ],
            States.MATH:[
                 #END MATH
                     (re.compile(r"(.*?)\</div\>\</html\>$"),
                      "    \\1",
                      LineTypes.LATEX,
                      States.NORMAL),
                 #CONTINUE MATH
                     (re.compile("(.*)"),
                      "    \\1",
                      LineTypes.LATEX,
                      States.MATH),        
                ],
            States.TRACEBACK:[
                 #END Traceback
                     (re.compile(r"^(\S.*)"),
                      "    ...\n    \\1",
                      LineTypes.TRACE,
                      States.NORMAL),
                ],
            States.HTML:[
                 #END HTML
                     (re.compile(r".*</html\>$"),
                      "",
                      LineTypes.HTML,
                      States.NORMAL),
                ],
        }
    
    def parse(self, text):
        result_plain = []
        result_show = []
        state = States.NORMAL
        for line in text.splitlines():
            for regex, replacement, line_type, new_state in self.transitions[state]:
                if regex.match(line):
                    result = result_plain if line_type in (LineTypes.PLAIN, LineTypes.HTML)\
                             else result_show
                    result.append( regex.sub(replacement, line))
                    state = new_state
                    break
        result_plain.extend(result_show)
        return '\n'.join(result_plain)

def results2rst(text, images_dir):
    r"""Converts the result of evaluation of notebook cells
    into rst compatible with Sage documentation.

    Several common patterns are identified, and treated
    accordingly. Some patterns are dropped, while others
    are not recognized.

    Currently, latex and images are recognized and converted.

    INPUT:

    - ``text`` -- string -- a chunk of HTML text

    - ``images_dir`` -- string -- folder where images are stored

    OUTPUT:

    - string -- rst text

    EXAMPLES::

        >>> from sage_sws2rst.results2rst import results2rst
        >>> s="<html><font color='black'><img src='cell://sage0.png'></font></html>"
        >>> results2rst(s,'')
        '\n.. image:: sage0.png\n    :align: center\n'
        >>> results2rst("4",'')
        '    4    '
        >>> s=r'<html><div class="math">\newcommand{\Bold}[1]{\mathbf{#1}}\frac{3}{2}</div></html>'
        >>> results2rst(s,'')
        '\n.. MATH::\n\n    \\frac{3}{2}\n'
    """
    Parser = ResultsParser(images_dir)
    return Parser.parse(text)

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
