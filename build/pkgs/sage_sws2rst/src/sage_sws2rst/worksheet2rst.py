#!/usr/bin/python
# -*- coding: utf-8 -*-
r"""
Convert worksheet.html files into ReStructuredText documents

This is called by 'sage -sws2rst'. Can also be used as a commandline script 
(if BeautifulSoup is installed):

``python worksheet2rst.py worksheet.html``

or

``cat worksheet.html | python worksheet2rst.py``

AUTHOR:

- Pablo Angulo Ardoy (2011-02-25): initial version


The content of worksheet.html is split into comments, code, and output 
(the result of evaluating the code), as follows:

comments
{{{id=..|
code
///
results
}}}

Each kind of text is dealt with separately.

"""

#**************************************************
# Copyright (C) 2011 Pablo Angulo
#
# Distributed under the terms of the GPL License
#**************************************************


import sys
import os
import re
from .comments2rst import html2rst
from .results2rst import results2rst
import codecs

#We parse lines one by one but keep track of current scope
#comments
#{{{id=..|
#code
#///
#results
#}}}
#RESULT_TO_BE_DROPPED corresponds to a results section whose
#code was empty, and will be discarded, whether it's empty or not
class States(object):
    COMMENT = 0
    CODE = 1
    RESULT = 2
    RESULT_TO_BE_DROPPED = 3

# REs for splitting comments, code and results
START_CELL_RE = re.compile(r'^\{\{\{id=(\d*)\|')
END_CODE_RE   = re.compile(r'^\/\/\/')
END_CELL_RE   = re.compile(r'^\}\}\}')

#When to switch State, and which State to
transitions = {
    States.COMMENT:(
        START_CELL_RE,
        States.CODE
        ),
    States.CODE:(
        END_CODE_RE,
        States.RESULT),
    States.RESULT:(
        END_CELL_RE,
        States.COMMENT),
    States.RESULT_TO_BE_DROPPED:(
        END_CELL_RE,
        States.COMMENT)
    }

def code_parser(text):
    """
    
    Arguments:

    INPUT:

    - ``s``:sage code, may or may not start with "sage:"

    OUTPUT:

    - string -- rst text

    EXAMPLES (not used for unit test, see 
    http://groups.google.com/group/sage-devel/browse_thread/thread/d82cb049ac102f3a)

    : from sage_sws2rst.worksheet2rst import code_parser
    : s="a=2"
    : code_parser(s)
    '::\n\n    sage: a=2'
    : s="def f(n):\n    return n+1\n"
    : code_parser(s)
    '::\n\n    sage: def f(n):\n    ....:     return n+1'
    : s="sage: def f(n):\nsage:     return n+1\n"
    : code_parser(s)
    '::\n\n    sage: def f(n):\n    ....:     return n+1'
    """
    lines = ['::', '']
    for s in text.splitlines():
        l = s[6:] if s.startswith('sage: ') else s
        if not l: continue
        prefix = '    ....: ' if l[0] == ' ' else '    sage: '
        lines.append(prefix + l)
    return '\n'.join(lines)

HEADER_RE = re.compile(r'<h\d>')
def add_title_if_there_is_none(text):
    if not HEADER_RE.search(text):
        return '<h1>Please write a title for this worksheet!</h1>\n' + text
    else:
        return text

def worksheet2rst(s, images_dir=''):
    """Parses a string, tipically the content of the file
    worksheet.html inside a sws file, and converts it into
    rst compatible with Sage documentation.

    INPUT:

    - ``s`` -- string -- text, tipically the content of
                               worksheet.html

    - ``images_dir`` -- string -- folder where images are stored

    OUTPUT:

    - string -- rst text

    EXAMPLES (not used for unit test, see 
    http://groups.google.com/group/sage-devel/browse_thread/thread/d82cb049ac102f3a)

    : from sage_sws2rst.worksheet2rst import worksheet2rst
    : worksheet2rst('<p>some text</p>\n{{{id=1|\nprint 2+2\n///\n4\n}}}')
    u'.. -*- coding: utf-8 -*-\n\nPlease write a title for this worksheet!\n========================================\n\nsome text\n\n\n::\n\n    sage: print 2+2\n    4\n\n.. end of output\n'
    : s = '{{{id=2|\nshow(f)\n///\n<html><div class="math">\\sqrt{x}</div></html>\n}}}\n'
    : worksheet2rst(s)
    u'.. -*- coding: utf-8 -*-\n\nPlease write a title for this worksheet!\n========================================\n::\n\n    sage: show(f)\n\n.. MATH::\n\n    \\sqrt{x}\n\n.. end of output\n'
    """
    s = add_title_if_there_is_none(s)
    state = States.COMMENT
    result = ['.. -*- coding: utf-8 -*-\n']
    ls = []
    for line in s.splitlines():
        regex, next_state= transitions[state]
        m = regex.match(line) 
        if m:
            if state == States.COMMENT:
                last_cell_id = m.group(1)
                img_path = images_dir + os.path.sep
                result.append(html2rst('\n'.join(ls), img_path))
            elif state == States.RESULT:
                img_path = os.path.join(images_dir, 'cell_%s_' % last_cell_id)
                result.append(results2rst('\n'.join(ls),
                                             img_path))
                result.append('')
                result.append('.. end of output')
            elif state == States.CODE:
                if ls and any(ls):
                    result.append(code_parser('\n'.join(ls)))
                else:
                    next_state = States.RESULT_TO_BE_DROPPED
            ls = []
            state = next_state
        else:
            ls.append(line)
    if state == States.COMMENT:
        img_path = images_dir + os.path.sep
        result.append(html2rst('\n'.join(ls), img_path))
    elif state == States.RESULT:
        img_path = os.path.join(images_dir, 'cell_%s_' % last_cell_id)
        result.append(result_parser('\n'.join(ls),
                                     img_path))
        result.append('')
        result.append('.. end of output')
    elif state == States.CODE:
        result.append(code_parser('\n'.join(ls)))

    return '\n'.join(result)

if __name__=='__main__':
    if len(sys.argv)>1:        
        fichero = codecs.open(sys.argv[1], mode='r', encoding='utf-8')
        text = fichero.read()
        fichero.close()
    else:
        text = sys.stdin.read()
    images_dir = sys.argv[2] if len(sys.argv)>2 else ''

    print((worksheet2rst(text, images_dir).encode('utf-8')))

