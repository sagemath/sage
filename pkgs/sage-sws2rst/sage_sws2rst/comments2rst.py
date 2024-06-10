# -*- coding: utf-8 -*-
r"""
Convert html from text cells in the notebook into ReStructuredText

This is called by sws2rst

- Pablo Angulo Ardoy (2011-02-25): initial version
"""
#**************************************************
# Copyright (C) 2011 Pablo Angulo
#
# Distributed under the terms of the GPL License
#**************************************************

import re
import os
try:
    from bs4 import (BeautifulSoup, Tag,
                     CData, Comment, Declaration, ProcessingInstruction)
except ImportError:
    raise ImportError("""BeautifulSoup must be installed.

Please either install using

    sage -pip install beautifulsoup4
""")

#negative lookbehind: http://www.regular-expressions.info/lookaround.html
double_dollar = re.compile(r'(?<!\\)\$\$')


def preprocess_display_latex(text):
    r"""replace $$some display latex$$ with <display>some display latex</display>
    before the soup is built.

    Deals with the situation where <p></p> tags are mixed
    with $$, like $$<p>display_latex$$</p>, unless the mess is huge

    EXAMPLES::

        >>> from sage_sws2rst.comments2rst import preprocess_display_latex
        >>> s="$$a=2$$"
        >>> preprocess_display_latex(s)
        '<display>a=2</display>' 
        >>> s="<p>$$a=2$$</p>"
        >>> preprocess_display_latex(s)
        '<p><display>a=2</display></p>'
        >>> s="<p>$$a=2</p>$$"
        >>> preprocess_display_latex(s)
        '<p><display>a=2</display></p>'
        >>> s="$$<p>a=2</p>$$"
        >>> preprocess_display_latex(s)
        '<display>a=2</display>'
    """
    ls = []
    start_tag = True
    parts = double_dollar.split(text)
    for c in parts[:-1]:
        if start_tag:
            ls.append(c)
            ls.append('<display>')
        else:
            c0, count = prune_tags(c)
            ls.append(c0)
            ls.append('</display>')
            if count == 1:
                ls.append('<p>')
            elif count == -1:
                ls.append('</p>')
            elif abs(count)>1:
                raise Exception('display latex was messed up with html code')
        start_tag = not start_tag
    ls.append(parts[-1])
    return ''.join(ls)


def prune_tags(text):
    count = text.count('<p>') - text.count('</p>')
    return text.replace('<br/>','').replace('<br />','').replace('<p>','').replace('</p>',''), count

escapable_chars = { '+' :r'\+',
                    '*' :r'\*',
                    '|' :r'\|',
                    '-' :r'\-'}


def escape_chars(text):
    for c, r in escapable_chars.items():
        text = text.replace(c, r)
    return text

#This is supposed to be handled by BeautifulSoup, but doesn't work
xml_entities = {'&lt;':'<',
            '&gt;':'>',
            '&amp;':'&',
            '&quot;':'"',
            '&apos;':"'",
}


def replace_xml_entities(text):
    for c, r in xml_entities.items():
        text = text.replace(c, r)
    return text
 

def replace_courier(soup):
    """Lacking a better option, I use courier font to mark <code>
    within tinyMCE. And I want to turn that into real code tags.

    Most users won't be needing this(?), so this code is not called anywhere
    but kept for reference
    """
    for t in soup.findAll(lambda s: ('style' in s) and 'courier' in s['style']):
        tag = Tag(soup, 'code')
        while t.contents:
            tag.append(t.contents[0])
        t.replaceWith(tag)

#negative lookbehind: http://www.regular-expressions.info/lookaround.html
single_dollar = re.compile(r'(?<!\\)\$')
def replace_latex(soup):
    r"""Replaces inline latex by :math:`code` and escapes
    some rst special chars like +, -, * and | outside of inline latex

    does not escape chars inside display or pre tags

    EXAMPLES::

        >>> from sage_sws2rst.comments2rst import replace_latex
        >>> from bs4 import BeautifulSoup
        >>> soup = r"<p>Some <strong>latex: $e^\pi i=-1$</strong></p>"
        >>> s = BeautifulSoup(soup, features='html.parser')
        >>> replace_latex(s)
        >>> s
        <p>Some <strong>latex: :math:`e^\pi i=-1`</strong></p>

    ::

        >>> soup = "<p><strong>2+2 | 1+3</strong></p>"
        >>> s = BeautifulSoup(soup, features='html.parser')
        >>> replace_latex(s)
        >>> s
        <p><strong>2\+2 \| 1\+3</strong></p>
    """
    for t in soup.findAll(text=re.compile('.+')):
        if (t.fetchParents(name = 'display') or
            t.fetchParents(name = 'pre')        ):
            continue
        parts = single_dollar.split(t)
        even  = [escape_chars(parts[i]) for i in range(0,len(parts),2)]
        odd   = [' :math:`%s`'%parts[i] for i in range(1,len(parts),2)]
        odd.append('')
        t.replaceWith(''.join(''.join(p) for p in zip(even,odd) ))


class Soup2Rst(object):
    """builds the rst text from the Soup Tree
    """
    tags = {'h1':'header',
            'h2':'header',
            'h3':'header',
            'h4':'header',
            'h5':'header',
            'h6':'header',
            'p': 'p',
            '[document]': 'document',
            'address': 'em',
            'br': 'br',
            'b':'strong',
            'strong':'strong',
            'em':'em',
            'pre':'pre',
            'code':'code',
            'display':'display',
            'span':'inline_no_tag',
            'ul':'ul',
            'ol':'ol',
            'li':'li',
            'a':'a',
            'table':'table',
#            'tr':'tr',
            'td':'inline_no_tag',
            'th':'inline_no_tag',
            'tt':'inline_no_tag',
            'div':'block_no_tag',
            'img':'img',
#            '':'',
            }

    headers = {'h1':'=',
               'h2':'-',
               'h3':'^',
               'h4':'"',
               'h5':'~',
               'h6':'*',
               }
    
    def __init__(self, images_dir):
        self.images_dir = images_dir
        self._nested_list = -1
        self._inside_ol_or_ul = []
        self._inside_code_tag = False

    def visit(self, node):
        if isinstance(node, (CData, Comment, Declaration, ProcessingInstruction)):
            return ''
        elif hasattr(node, 'name') and node.name in self.tags:
                method = 'visit_' + self.tags[node.name]
                visitor = getattr(self, method)
                return visitor(node)
        else:
            #Assume plain string
            return str(node).replace('\n','')

    def visit_document(self, node):
        return '\n'.join(self.visit(tag) for tag in node.contents)    

    def get_plain_text(self, node):
        """Gets all text, removing all tags"""
        if hasattr(node, 'contents'):
            t = ' '.join(self.get_plain_text(tag) for tag in node.contents)
        else:
            t = str(node)
        return t.replace('\n','')
        
    def visit_header(self, node):
        s = ''.join(self.visit(tag) for tag in node.contents)
        spacer = self.headers[node.name]*len(s)
        return s.replace( '\n', '') +  '\n' + spacer

    def visit_pre(self, node):
        return '::\n\n    '+str(node)[5:-6].replace('<br />','\n').replace('<br></br>','\n').replace('\n','\n    ')

    def visit_ul(self, node):
        self._nested_list += 1
        self._inside_ol_or_ul.append(False)
        result = '\n\n'+''.join(self.visit(tag) for tag in node.contents)+'\n'
        self._inside_ol_or_ul.pop()
        self._nested_list -= 1
        return result

    def visit_ol(self, node):
        self._nested_list += 1
        self._inside_ol_or_ul.append(True)
        result = '\n\n'+''.join(self.visit(tag) for tag in node.contents)+'\n'
        self._inside_ol_or_ul.pop()
        self._nested_list -= 1
        return result

    def visit_li(self, node):
        return (' '*self._nested_list
                + ('#. ' if self._inside_ol_or_ul[-1] else '- ')
                +' '.join(self.visit(tag) for tag in node.contents)
                + '\n')

    def visit_display(self, node):
        return ('\n\n.. MATH::\n\n    ' +
                str(node)[9:-10].replace('<br></br>','\n').replace('\n','\n    ') +
                '\n\n.. end of math\n\n')

    def visit_img(self, node):
        return '.. image:: ' + os.path.join(self.images_dir, node['src'].replace(' ','_')) + '\n    :align: center\n'

    def visit_table(self,node):
        rows = []
        for elt in node.contents:
            if not hasattr(elt,'name'):
                pass
            elif elt.name == 'thead':
                rows.extend(self.prepare_tr(row)
                            for row in elt
                            if hasattr(row,'name') and
                            row.name=='tr')
                rows.append([]) #this row represents a separator
            elif (elt.name == 'tbody') or (elt.name == 'tfoot'):
                rows.extend(self.prepare_tr(row)
                            for row in elt
                            if hasattr(row,'name') and
                            row.name=='tr')
            elif elt.name == 'tr':
                rows.append(self.prepare_tr(elt))

        ncols = max(len(row) for row in rows)
        for row in rows:
            if len(row) < ncols:
                row.extend( ['']*(ncols - len(row)))
        cols_sizes = [max(len(td) for td in tds_in_col)
                      for tds_in_col in zip(*rows)]
        result = [' '.join('='*c for c in cols_sizes)]
        
        for row in rows:
            if any(td for td in row):
                result.append(' '.join(td+' '*(l - len(td))
                                       for l,td in zip(cols_sizes,row)))
            else:
                result.append(' '.join('-'*c for c in cols_sizes))
        result.append(' '.join('='*c for c in cols_sizes))
        return '\n'.join(result)

    def prepare_tr(self, node):
        return [self.visit(tag) for tag in node.contents if tag!='\n']
        
    def visit_br(self, node):
        return '\n\n'

    def visit_strong(self, node):
        if node.contents:
            content = ' '.join(self.visit(tag) for tag in node.contents).strip()
            if not content:
                return ''
            elif '``' in content:
                return content
            else:
                return ' **' + content + '** '
        else:
            return ''

    def visit_em(self,node):
        if node.contents:
            content = ' '.join(self.visit(tag) for tag in node.contents).strip()
            if not content:
                return ''
            elif '``' in content:
                return content
            else:
                return ' *' + content + '* '
        else:
            return ''

    def visit_code(self, node):
        if node.contents:
            content = self.get_plain_text(node).strip()
            return '``' + content + '``'
        else:
            return ''

    def visit_inline_no_tag(self, node):
        return (' '.join(self.visit(tag)
                         for tag in node.contents)).strip()

    def visit_block_no_tag(self, node):
        return '\n'.join(self.visit(tag) for tag in node.contents) + '\n'

    def visit_p(self, node):
        return ''.join(self.visit(tag) for tag in node.contents) + '\n\n'

    def visit_a(self, node):
        c = ' '.join(self.visit(tag) for tag in node.contents)
        try:
            link = node['href']
            if link[0]=='#':
                return ':ref:`%s <%s>`'%(c, link[1:])
            else:                    
                return '`%s <%s>`_'%(c, link)
        except KeyError:
            return '.. _%s:\n\n'%node['name']


def html2rst(text, images_dir):
    r"""
    Convert html, typically generated by tinyMCE, into rst
    compatible with Sage documentation.

    The main job is done by BeautifulSoup, which is much more
    robust than conventional parsers like HTMLParser, but also
    several details specific of this context are taken into
    account, so this code differs from generic approaches like
    those found on the web.

    INPUT:

    - ``text`` -- string -- a chunk of HTML text

    - ``images_dir`` -- string -- folder where images are stored

    OUTPUT:

    - string -- rst text

    EXAMPLES::

        >>> from sage_sws2rst.comments2rst import html2rst
        >>> text = r'<p>Some text with <em>math</em>: $e^{\pi i}=-1$</p>'
        >>> html2rst(text, '')
        'Some text with  *math* :  :math:`e^{\\pi i}=-1`\n\n'

    ::

        >>> text = '<p>Text with <em>incorrect</p> nesting</em>.'
        >>> html2rst(text, '')
        'Text with  *incorrect* \n\n nesting\n.'

    ::

        >>> text = '<pre>Preformatted: \n    a+2\n</pre><p> Not preformatted: \n    a+2\n</p>'
        >>> html2rst(text, '')
        '::\n\n    Preformatted: \n        a+2\n    \n Not preformatted:     a\\+2\n\n'

    ::

        >>> text = '&aacute;ñ&nbsp;&ntildeá'
        >>> html2rst(text, '')
        '\xe1\xf1 \xf1\xe1'

    ::

        >>> text = r'<p>some text</p><p>$$</p><p>3.183098861 \cdot 10^{-1}</p><p>$$</p>'
        >>> html2rst(text, '')
        'some text\n\n.. MATH::\n\n    3.183098861 \\cdot 10^{-1}\n\n.. end of math\n\n'

    When the content is empty::

        >>> html2rst("<strong></strong> ", '')
        '\n '
        >>> html2rst("<strong> </strong> ", '')
        '\n '
        >>> html2rst("<em></em> ", '')
        '\n '
        >>> html2rst("<em> </em> ", '')
        '\n '

    Spaces are added around *italic* or **bold** text (otherwise, it
    may be invalid ReStructuredText syntax)::

        >>> text = '<p><strong>Exercice.</strong>Let x be ...</p>'
        >>> html2rst(text, '')
        ' **Exercice.** Let x be ...\n\n'
        >>> text = '<p><em>Exercice.</em>Let x be ...</p>'
        >>> html2rst(text, '')
        ' *Exercice.* Let x be ...\n\n'

    Below is an example showing the translation from html to rst is not
    always perfect.

    Here the strong emphasis is on more than one line and is not properly
    translated::

        >>> text='<p>You will find a <em>while loop</em> helpful here. Below is a simple example:</p><p style="padding-left: 30px;"><strong>x = 0<br />while x &lt; 7:<br />&nbsp;&nbsp;&nbsp; x = x + 2<br />&nbsp;&nbsp;&nbsp; print x</strong></p>'
        >>> html2rst(text, '')
        'You will find a  *while loop*  helpful here. Below is a simple
        example:\n\n **x = 0 \n\n while x < 7: \n\n     x = x \\+ 2 \n\n
        print x** \n\n'

    """
    
    #replace $$some display latex$$ with
    #<display>some display latex</display>
    text = preprocess_display_latex(text)

    #eliminate nasty &nbsp;
    text = text.replace('&nbsp;',' ')
            
    #BeautifulSoup is better than BeautifulSoup
    #for html that wasn't generated by humans (like tinyMCE)
    soup = BeautifulSoup(text,
                         features='html.parser',
                         # https://stackoverflow.com/questions/11856011/beautifulsoup-has-no-attribute-html-entities
                         ##convertEntities=BeautifulSoup.ALL_ENTITIES
                         )        

    #remove all comments
    comments = soup.findAll(text=lambda text:isinstance(text, Comment))
    for comment in comments:
        comment.extract()

#    replace_courier(soup)
    replace_latex(soup)
    v = Soup2Rst(images_dir)

#    return v.visit(soup)
    text = v.visit(soup)
    more_than_2_blank_lines = re.compile(r'\n\n+', re.MULTILINE)
    text = more_than_2_blank_lines.sub('\n\n', text)
    text = replace_xml_entities(text)
    return text
    
if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
