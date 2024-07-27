"Interface to the free online MAGMA calculator"

# ****************************************************************************
# Copyright (C) 2007 William Stein <wstein@gmail.com>
#
# Distributed under the terms of the GNU General Public License (GPL)
#
# This code is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# The full text of the GPL is available at:
#
# https://www.gnu.org/licenses/
# ****************************************************************************


class MagmaExpr(str):
    def __repr__(self):
        return str(self)

def magma_free_eval(code, strip=True, columns=0):
    """
    Use the free online MAGMA calculator to evaluate the given input code 
    and return the answer as a string.

    INPUT:
    - ``code`` -- string; MAGMA code to evaluate
    - ``strip`` -- boolean (default: True); whether to strip whitespace from the output
    - ``columns`` -- integer (default: 0); number of columns for output formatting

    OUTPUT:
    - A string containing the output from the MAGMA online calculator.

    LIMITATIONS: 
    - The code must evaluate in at most 120 seconds.
    - There is a limitation on the amount of RAM available.

    EXAMPLES::

        sage: magma_free_eval("Factorization(9290348092384)") # optional - internet
        [ <2, 5>, <290323377887, 1> ]
    """

    from urllib.parse import urlencode
    from http import client as httplib
    from xml.dom.minidom import parseString
    
    server = "magma.maths.usyd.edu.au"
    processPath = "/xml/calculator.xml"
    refererPath = "/calc/"
    refererUrl = f"http://{server}{refererPath}"
    
    code = f"SetColumns({columns});\n{code}"
    params = urlencode({'input': code})
    headers = {
        "Content-type": "application/x-www-form-urlencoded",
        "Accept": "text/html, application/xml, application/xhtml+xml",
        "Referer": refererUrl
    }
    
    conn = httplib.HTTPConnection(server)
    conn.request("POST", processPath, params, headers)
    response = conn.getresponse()
    results = response.read()
    conn.close()
    
    xmlDoc = parseString(results)
    res = []
    resultsNodeList = xmlDoc.getElementsByTagName('results')
    if len(resultsNodeList) > 0:
        resultsNode = resultsNodeList[0]
        lines = resultsNode.getElementsByTagName('line')
        for line in lines:
            for textNode in line.childNodes:
                res.append(textNode.data)
    
    res = "\n".join(res)
    if strip:
        res = res.strip()
    
    return MagmaExpr(res)

class MagmaFree:
    """
    Interface to the free online MAGMA calculator.

    Evaluate MAGMA code without requiring MAGMA to be installed on your computer
    by using the free online MAGMA calculator.

    LIMITATIONS: 
    - The code must evaluate in at most 120 seconds.
    - There is a limitation on the amount of RAM available.

    INPUT:
    - ``code`` -- string; MAGMA code to evaluate
    - ``strip`` -- boolean (default: True); whether to strip whitespace from the output
    - ``columns`` -- integer (default: 0); number of columns for output formatting

    OUTPUT:
    - A string containing the output from the MAGMA online calculator.

    EXAMPLES::

        sage: magma_free("Factorization(9290348092384)") # optional - internet
        [ <2, 5>, <290323377887, 1> ]
    """

    def __call__(self, code, strip=True, columns=0):
        """
        Evaluate MAGMA code using the free online calculator.

        See the class documentation for more details and examples.
        """
        return magma_free_eval(code, strip=strip, columns=columns)

magma_free = MagmaFree()
