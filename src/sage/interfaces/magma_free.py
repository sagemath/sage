"Interface to the free online MAGMA calculator"

# ****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************


class MagmaExpr(str):
    def __repr__(self) -> str:
        return str(self)


def magma_free_eval(code: str, strip=True, columns=0):
    """
    Use the free online MAGMA calculator to evaluate the given
    input code and return the answer as a string.

    .. WARNING::

        The code must evaluate in at most 120 seconds
        and there is a limitation on the amount of RAM.
        Otherwise, some :exc:`TimeoutError` will be raised.

    EXAMPLES::

        sage: magma_free("Factorization(9290348092384)")  # optional - internet
        [ <2, 5>, <290323377887, 1> ]
    """
    from urllib.parse import urlencode
    from http import client as httplib
    from xml.dom.minidom import parseString

    server = "magma.maths.usyd.edu.au"
    processPath = "/xml/calculator.xml"
    refererPath = "/calc/"
    refererUrl = "http://%s%s" % (server, refererPath)
    code = "SetColumns(%s);\n" % columns + code
    params = urlencode({'input': code})
    headers = {"Content-type": "application/x-www-form-urlencoded",
               "Accept": "Accept: text/html, application/xml, application/xhtml+xml", "Referer": refererUrl}
    conn = httplib.HTTPConnection(server)
    conn.request("POST", processPath, params, headers)
    response = conn.getresponse()
    results = response.read()
    conn.close()

    if b"Timeout" in results:
        raise TimeoutError('timeout from the server')

    xmlDoc = parseString(results)
    res: list[str] = []
    resultsNodeList = xmlDoc.getElementsByTagName('results')
    if resultsNodeList:
        resultsNode = resultsNodeList[0]
        lines = resultsNode.getElementsByTagName('line')
        for line in lines:
            res.extend(textNode.data for textNode in line.childNodes)
    return MagmaExpr("\n".join(res))


class MagmaFree:
    """
    Evaluate MAGMA code without requiring that MAGMA be installed
    on your computer by using the free online MAGMA calculator.

    EXAMPLES::

        sage: magma_free("Factorization(9290348092384)")  # optional - internet
        [ <2, 5>, <290323377887, 1> ]
    """
    def eval(self, x, **kwds):
        return magma_free_eval(x)

    def __call__(self, code, strip=True, columns=0):
        return magma_free_eval(code, strip=strip, columns=columns)


magma_free = MagmaFree()
