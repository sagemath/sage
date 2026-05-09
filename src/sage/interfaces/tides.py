r"""
Interface to TIDES

This module contains tools to write the .c files needed for TIDES [TIDES]_ .

Tides is an integration engine based on the Taylor method. It is implemented
as a c library. The user must translate its initial value problem (IVP) into a
pair of .c files that will then be compiled and linked against the TIDES
library. The resulting binary will produce the desired output. The tools in this
module can be used to automate the generation of these files from the symbolic
expression of the differential equation.

::

    # ************************************************************************
    #  Copyright (C) 2014 Miguel Marco <mmarco@unizar.es>, Marcos Rodriguez
    #   <marcos@uunizar.es>
    #
    #  Distributed under the terms of the GNU General Public License (GPL):
    #
    #                  https://www.gnu.org/licenses/
    # ************************************************************************

AUTHORS:

- Miguel Marco (06-2014) - Implementation of tides solver

- Marcos Rodriguez (06-2014) - Implementation of tides solver

- Alberto Abad (06-2014) - tides solver

- Roberto Barrio (06-2014) - tides solver

REFERENCES:

- [ABBR2012]_

- [TIDES]_
"""

from sage.rings.real_mpfr import RealField
from sage.misc.lazy_import import lazy_import
lazy_import("sage.calculus.all", "symbolic_expression")
from sage.misc.flatten import flatten
from sage.ext.fast_callable import fast_callable
from sage.rings.semirings.non_negative_integer_semiring import NN
from sage.functions.log import log, exp
from sage.functions.other import floor, ceil
from sage.misc.functional import sqrt


def subexpressions_list(f, pars=None):
    """
    Construct the lists with the intermediate steps on the evaluation of the
    function.

    INPUT:

    - ``f`` -- a symbolic function of several components

    - ``pars`` -- list of the parameters that appear in the function
      this should be the symbolic constants that appear in f but are not
      arguments

    OUTPUT:

    - a list of the intermediate subexpressions that appear in the evaluation
      of f.

    - a list with the operations used to construct each of the subexpressions.
      each element of this list is a tuple, formed by a string describing the
      operation made, and the operands.

    For the trigonometric functions, some extra expressions will be added.
    These extra expressions will be used later to compute their derivatives.

    EXAMPLES::

        sage: from sage.interfaces.tides import subexpressions_list
        sage: var('x,y')
        (x, y)
        sage: f(x,y) = [x^2+y, cos(x)/log(y)]
        sage: subexpressions_list(f)
        ([x^2, x^2 + y, sin(x), cos(x), log(y), cos(x)/log(y)],
        [('mul', x, x),
        ('add', y, x^2),
        ('sin', x),
        ('cos', x),
        ('log', y),
        ('div', log(y), cos(x))])

    ::

        sage: f(a)=[cos(a), arctan(a)]
        sage: from sage.interfaces.tides import subexpressions_list
        sage: subexpressions_list(f)
        ([sin(a), cos(a), a^2, a^2 + 1, arctan(a)],
        [('sin', a), ('cos', a), ('mul', a, a), ('add', 1, a^2), ('atan', a)])

    ::

        sage: from sage.interfaces.tides import subexpressions_list
        sage: var('s,b,r')
        (s, b, r)
        sage: f(t,x,y,z)= [s*(y-x),x*(r-z)-y,x*y-b*z]
        sage: subexpressions_list(f,[s,b,r])
        ([-y,
        x - y,
        s*(x - y),
        -s*(x - y),
        -z,
        r - z,
        (r - z)*x,
        -y,
        (r - z)*x - y,
        x*y,
        b*z,
        -b*z,
        x*y - b*z],
        [('mul', -1, y),
        ('add', -y, x),
        ('mul', x - y, s),
        ('mul', -1, s*(x - y)),
        ('mul', -1, z),
        ('add', -z, r),
        ('mul', x, r - z),
        ('mul', -1, y),
        ('add', -y, (r - z)*x),
        ('mul', y, x),
        ('mul', z, b),
        ('mul', -1, b*z),
        ('add', -b*z, x*y)])

    ::

        sage: var('x, y')
        (x, y)
        sage: f(x,y)=[exp(x^2+sin(y))]
        sage: from sage.interfaces.tides import *
        sage: subexpressions_list(f)
        ([x^2, sin(y), cos(y), x^2 + sin(y), e^(x^2 + sin(y))],
        [('mul', x, x),
        ('sin', y),
        ('cos', y),
        ('add', sin(y), x^2),
        ('exp', x^2 + sin(y))])
    """
    from sage.functions.trig import sin, cos, arcsin, arctan, arccos
    variables = f[0].arguments()
    if not pars:
        parameters = []
    else:
        parameters = pars
    varpar = list(parameters) + list(variables)
    F = symbolic_expression([i(*variables) for i in f]).function(*varpar)
    lis = flatten([fast_callable(i,vars=varpar).op_list() for i in F], max_level=1)
    stack = []
    const = []
    stackcomp = []
    detail = []
    for i in lis:
        if i[0] == 'load_arg':
            stack.append(varpar[i[1]])
        elif i[0] == 'ipow':
            if i[1] in NN:
                basis = stack[-1]
                for j in range(i[1]-1):
                    a = stack.pop(-1)
                    detail.append(('mul', a, basis))
                    stack.append(a*basis)
                    stackcomp.append(stack[-1])
            else:
                detail.append(('pow',stack[-1],i[1]))
                stack[-1] = stack[-1]**i[1]
                stackcomp.append(stack[-1])

        elif i[0] == 'load_const':
            const.append(i[1])
            stack.append(i[1])
        elif i == 'mul':
            a = stack.pop(-1)
            b = stack.pop(-1)
            detail.append(('mul', a, b))
            stack.append(a*b)
            stackcomp.append(stack[-1])

        elif i == 'div':
            a = stack.pop(-1)
            b = stack.pop(-1)
            detail.append(('div', a, b))
            stack.append(b/a)
            stackcomp.append(stack[-1])

        elif i == 'add':
            a = stack.pop(-1)
            b = stack.pop(-1)
            detail.append(('add',a,b))
            stack.append(a+b)
            stackcomp.append(stack[-1])

        elif i == 'pow':
            a = stack.pop(-1)
            b = stack.pop(-1)
            detail.append(('pow', b, a))
            stack.append(b**a)
            stackcomp.append(stack[-1])

        elif i[0] == 'py_call' and str(i[1]) == 'log':
            a = stack.pop(-1)
            detail.append(('log', a))
            stack.append(log(a))
            stackcomp.append(stack[-1])

        elif i[0] == 'py_call' and str(i[1]) == 'exp':
            a = stack.pop(-1)
            detail.append(('exp', a))
            stack.append(exp(a))
            stackcomp.append(stack[-1])

        elif i[0] == 'py_call' and str(i[1]) == 'sin':
            a = stack.pop(-1)
            detail.append(('sin', a))
            detail.append(('cos', a))
            stackcomp.append(sin(a))
            stackcomp.append(cos(a))
            stack.append(sin(a))

        elif i[0] == 'py_call' and str(i[1]) == 'cos':
            a = stack.pop(-1)
            detail.append(('sin', a))
            detail.append(('cos', a))
            stackcomp.append(sin(a))
            stackcomp.append(cos(a))
            stack.append(cos(a))

        elif i[0] == 'py_call' and str(i[1]) == 'tan':
            a = stack.pop(-1)
            b = sin(a)
            c = cos(a)
            detail.append(('sin', a))
            detail.append(('cos', a))
            detail.append(('div', b, c))
            stackcomp.append(b)
            stackcomp.append(c)
            stackcomp.append(b/c)
            stack.append(b/c)

        elif i[0] == 'py_call' and str(i[1]) == 'arctan':
            a = stack.pop(-1)
            detail.append(('mul', a, a))
            detail.append(('add', 1, a*a))
            detail.append(('atan', a))
            stackcomp.append(a*a)
            stackcomp.append(1+a*a)
            stackcomp.append(arctan(a))
            stack.append(arctan(a))

        elif i[0] == 'py_call' and str(i[1]) == 'arcsin':
            a = stack.pop(-1)
            detail.append(('mul', a, a))
            detail.append(('mul', -1, a*a))
            detail.append(('add', 1, -a*a))
            detail.append(('pow', 1 - a*a, 0.5))
            detail.append(('asin', a))
            stackcomp.append(a*a)
            stackcomp.append(-a*a)
            stackcomp.append(1-a*a)
            stackcomp.append(sqrt(1-a*a))
            stackcomp.append(arcsin(a))
            stack.append(arcsin(a))

        elif i[0] == 'py_call' and str(i[1]) == 'arccos':
            a = stack.pop(-1)
            detail.append(('mul', a, a))
            detail.append(('mul', -1, a*a))
            detail.append(('add', 1, -a*a))
            detail.append(('pow', 1 - a*a, 0.5))
            detail.append(('mul', -1, sqrt(1-a*a)))
            detail.append(('acos', a))
            stackcomp.append(a*a)
            stackcomp.append(-a*a)
            stackcomp.append(1-a*a)
            stackcomp.append(sqrt(1-a*a))
            stackcomp.append(-sqrt(1-a*a))
            stackcomp.append(arccos(a))
            stack.append(arccos(a))

        elif i[0] == 'py_call' and 'sqrt' in str(i[1]):
            a = stack.pop(-1)
            detail.append(('pow', a, 0.5))
            stackcomp.append(sqrt(a))
            stack.append(sqrt(a))

        elif i == 'neg':
            a = stack.pop(-1)
            detail.append(('mul', -1, a))
            stack.append(-a)
            stackcomp.append(-a)

    return stackcomp,detail


def remove_repeated(l1, l2):
    """
    Given two lists, remove the repeated elements in l1, and the elements
    in l2 that are on the same position.
    positions.

    EXAMPLES::

        sage: from sage.interfaces.tides import (subexpressions_list, remove_repeated)
        sage: f(a)=[1 + a^2, arcsin(a)]
        sage: l1, l2 = subexpressions_list(f)
        sage: l1, l2
        ([a^2, a^2 + 1, a^2, -a^2, -a^2 + 1, sqrt(-a^2 + 1), arcsin(a)],
        [('mul', a, a),
        ('add', 1, a^2),
        ('mul', a, a),
        ('mul', -1, a^2),
        ('add', 1, -a^2),
        ('pow', -a^2 + 1, 0.5),
        ('asin', a)])
        sage: remove_repeated(l1, l2)
        sage: l1, l2
        ([a^2, a^2 + 1, -a^2, -a^2 + 1, sqrt(-a^2 + 1), arcsin(a)],
        [('mul', a, a),
        ('add', 1, a^2),
        ('mul', -1, a^2),
        ('add', 1, -a^2),
        ('pow', -a^2 + 1, 0.5),
        ('asin', a)])
    """
    for i in range(len(l1)-1):
        j = i+1
        while j < len(l1):
            if str(l1[j]) == str(l1[i]):
                l1.pop(j)
                l2.pop(j)
            else:
                j += 1


def remove_constants(l1, l2):
    """
    Given two lists, remove the entries in the first that are real constants,
    and also the corresponding elements in the second one.

    EXAMPLES::

        sage: from sage.interfaces.tides import subexpressions_list, remove_constants
        sage: f(a)=[1+cos(7)*a]
        sage: l1, l2 = subexpressions_list(f)
        sage: l1, l2
        ([sin(7), cos(7), a*cos(7), a*cos(7) + 1],
        [('sin', 7), ('cos', 7), ('mul', cos(7), a), ('add', 1, a*cos(7))])
        sage: remove_constants(l1,l2)
        sage: l1, l2
        ([a*cos(7), a*cos(7) + 1], [('mul', cos(7), a), ('add', 1, a*cos(7))])
    """
    i = 0
    while i < len(l1):
        if l1[i] in RealField():
            l1.pop(i)
            l2.pop(i)
        else:
            i += 1


def genfiles_mintides(integrator, driver, f, ics, initial, final, delta,
                      tolrel=1e-16, tolabs=1e-16, output=''):
    r"""
    Generate the needed files for the min_tides library.

    INPUT:

    - ``integrator`` -- the name of the integrator file

    - ``driver`` -- the name of the driver file

    - ``f`` -- the function that determines the differential equation

    - ``ics`` -- list or tuple with the initial conditions

    - ``initial`` -- the initial time for the integration

    - ``final`` -- the final time for the integration

    - ``delta`` -- the step of the output

    - ``tolrel`` -- the relative tolerance

    - ``tolabs`` -- the absolute tolerance

    - ``output`` -- the name of the file that the compiled integrator will write to

    This function creates two files, integrator and driver, that can be used
    later with the min_tides library [TIDES]_.


    TESTS::

        sage: from sage.interfaces.tides import genfiles_mintides
        sage: import os
        sage: import shutil
        sage: from sage.misc.temporary_file import tmp_dir
        sage: tempdir = tmp_dir()
        sage: intfile = os.path.join(tempdir, 'integrator.c')
        sage: drfile = os.path.join(tempdir ,'driver.c')
        sage: var('t,x,y,X,Y')
        (t, x, y, X, Y)
        sage: f(t,x,y,X,Y)=[X, Y, -x/(x^2+y^2)^(3/2), -y/(x^2+y^2)^(3/2)]
        sage: genfiles_mintides(intfile, drfile, f, [1,0, 0, 0.2], 0, 10, 0.1, output = 'out')
        sage: fileint = open(intfile)
        sage: l = fileint.readlines()
        sage: fileint.close()
        sage: l[5]
        '    #include "minc_tides.h"\n'
        sage: l[15]
        '    double XX[TT+1][MO+1];\n'
        sage: l[25]
        '\n'
        sage: l[35]
        '\t\tXX[1][i+1] = XX[3][i] / (i+1.0);\n'
        sage: filedr = open(drfile)
        sage: l = filedr.readlines()
        sage: filedr.close()
        sage: l[6]
        '    #include "minc_tides.h"\n'
        sage: l[15]
        '    double tolrel, tolabs, tini, tend, dt;\n'
        sage: l[25]
        '\ttolrel = 9.9999999999999998e-17 ;\n'
        sage: shutil.rmtree(tempdir)

    Check that issue :issue:`17179` is fixed (handle expressions like `\\pi`)::

        sage: from sage.interfaces.tides import genfiles_mintides
        sage: import os
        sage: import shutil
        sage: from sage.misc.temporary_file import tmp_dir
        sage: tempdir = tmp_dir()
        sage: intfile = os.path.join(tempdir, 'integrator.c')
        sage: drfile = os.path.join(tempdir ,'driver.c')
        sage: var('t,x,y,X,Y')
        (t, x, y, X, Y)
        sage: f(t,x,y,X,Y)=[X, Y, -x/(x^2+y^2)^(3/2), -y/(x^2+y^2)^(3/2)]
        sage: genfiles_mintides(intfile, drfile, f, [pi, 0, 0, 0.2], 0, 10, 0.1, output = 'out')
        sage: fileint = open(intfile)
        sage: l = fileint.readlines()
        sage: fileint.close()
        sage: l[30]
        '\t\tXX[8][i] = pow_mc_c(XX[7],-1.5000000000000000,XX[8], i);\n'
        sage: filedr = open(drfile)
        sage: l = filedr.readlines()
        sage: filedr.close()
        sage: l[18]
        '    \tv[0] = 3.1415926535897931 ; \n'
        sage: shutil.rmtree(tempdir)
    """
    RR = RealField()

    l1, l2 = subexpressions_list(f)

    remove_repeated(l1, l2)
    remove_constants(l1, l2)
    l0 = [str(l) for l in l1]
    # generate the corresponding c lines

    l3 = []
    var = f[0].arguments()
    lv = [str(v) for v in var]
    for i in l2:
        oper = i[0]
        if oper in ["log", "exp", "sin", "cos"]:
            a = i[1]
            if a in var:
                l3.append((oper, 'XX[{}]'.format(lv.index(str(a)))))
            elif a in l1:
                l3.append((oper, 'XX[{}]'.format(l0.index(str(a))+len(var))))

        else:
            a = i[1]
            b = i[2]
            consta = False
            constb = False

            if str(a) in lv:
                aa = 'XX[{}]'.format(lv.index(str(a)))
            elif str(a) in l0:
                aa = 'XX[{}]'.format(l0.index(str(a))+len(var))
            else:
                consta = True
                aa = RR(a).str()
            if str(b) in lv:
                bb = 'XX[{}]'.format(lv.index(str(b)))
            elif str(b) in l0:
                bb = 'XX[{}]'.format(l0.index(str(b))+len(var))
            else:
                constb = True
                bb = RR(b).str()
            if consta:
                oper += '_c'
                if not oper == 'div':
                    bb, aa = aa, bb
            elif constb:
                oper += '_c'
            l3.append((oper, aa, bb))

    n = len(var)
    res = []
    for i in range(len(l3)):
        el = l3[i]
        string = "XX[{}][i] = ".format(i + n)
        if el[0] == 'add':
            string += el[1] + "[i] + " + el[2] + "[i];"
        elif el[0] == 'add_c':
            string += "(i==0)? {}+".format(el[2]) + el[1] + "[0] : " + el[1] + "[i];"
        elif el[0] == 'mul':
            string += "mul_mc("+el[1]+","+el[2]+",i);"
        elif el[0] == 'mul_c':
            string += el[2] + "*" + el[1] + "[i];"
        elif el[0] == 'pow_c':
            string += "pow_mc_c("+el[1]+","+el[2]+",XX[{}], i);".format(i+n)
        elif el[0] == 'div':
            string += "div_mc("+el[2]+","+el[1]+",XX[{}], i);".format(i+n)
        elif el[0] == 'div_c':
            string += "inv_mc("+el[2]+","+el[1]+",XX[{}], i);".format(i+n)
        elif el[0] == 'log':
            string += "log_mc(" + el[1] + ",XX[{}], i);".format(i+n)
        elif el[0] == 'exp':
            string += "exp_mc(" + el[1] + ",XX[{}], i);".format(i+n)
        elif el[0] == 'sin':
            string += "sin_mc(" + el[1] + ",XX[{}], i);".format(i+n+1)
        elif el[0] == 'cos':
            string += "cos_mc(" + el[1] + ",XX[{}], i);".format(i+n-1)

        res.append(string)

    l0 = lv + l0
    indices = [l0.index(str(i(*var))) + n for i in f]
    for i in range(1, n):
        res.append("XX[{}][i+1] = XX[{}][i] / (i+1.0);".format(i,indices[i-1]-n))

    code = res

    outfile = open(integrator, 'a')
    auxstring = """
    /****************************************************************************
    This file has been created by Sage for its use with TIDES
    *****************************************************************************/

    #include "minc_tides.h"

    void    mincseries(double t,double *v, double *p, double **XVAR,int ORDER, int MO)
    {
    int VAR,PAR,TT,i,j, inext;
    """
    outfile.write(auxstring)

    outfile.write("\tVAR = {};\n".format(n))
    outfile.write("\tPAR = {};\n".format(0))
    outfile.write("\tTT = {};\n".format(len(res)))

    auxstring = """

    double XX[TT+1][MO+1];

    for(j=0; j<=TT; j++)
        for(i=0; i<=ORDER; i++)
            XX[j][i] = 0.e0;
    XX[0][0] = t;
    XX[0][1] = 1.e0;
    for(i=1;i<=VAR;i++) {
        XX[i][0] = v[i-1];
    }

    for(i=0;i<ORDER;i++) {
    """
    outfile.write(auxstring)
    outfile.writelines(["\t\t" + i + "\n" for i in code])

    outfile.write('\t}\n')
    outfile.write('\n')
    outfile.write('\tfor(j=0; j<=VAR; j++)\n')
    outfile.write('\t\tfor(i=0; i<=ORDER; i++)\n')
    outfile.write('\t\t\tXVAR[i][j] = XX[j][i];\n')
    outfile.write('}\n')
    outfile.write('\n')

    outfile = open(driver, 'a')

    auxstring = """
    /****************************************************************************
        Driver file of the minc_tides program
        This file has been automatically created by Sage
    *****************************************************************************/

    #include "minc_tides.h"

    int main() {

        int  i, VARS, PARS;


    VARS = %s ;
    PARS = 1;
    double tolrel, tolabs, tini, tend, dt;
    double v[VARS], p[PARS];

    """ % (n-1)
    outfile.write(auxstring)
    for i in range(len(ics)):
        outfile.write('\tv[{}] = {} ; \n'.format(i, RR(ics[i]).str()))
    outfile.write('\ttini = {} ;\n'.format(RR(initial).str()))
    outfile.write('\ttend = {} ;\n'.format(RR(final).str()))
    outfile.write('\tdt   = {} ;\n'.format(RR(delta).str()))
    outfile.write('\ttolrel = {} ;\n'.format(RR(tolrel).str()))
    outfile.write('\ttolabs = {} ;\n'.format(RR(tolabs).str()))
    outfile.write('\textern char ofname[500];')
    outfile.write('\tstrcpy(ofname, "' + output + '");\n')
    outfile.write('\tminc_tides(v,VARS,p,PARS,tini,tend,dt,tolrel,tolabs);\n')
    outfile.write('\treturn 0; \n }')
    outfile.close()


def genfiles_mpfr(integrator, driver, f, ics, initial, final, delta,
                  parameters=None, parameter_values=None, dig=20, tolrel=1e-16,
                  tolabs=1e-16, output=''):
    r"""
        Generate the needed files for the mpfr module of the tides library.

    INPUT:

    - ``integrator`` -- the name of the integrator file

    - ``driver`` -- the name of the driver file

    - ``f`` -- the function that determines the differential equation

    - ``ics`` -- list or tuple with the initial conditions

    - ``initial`` -- the initial time for the integration

    - ``final`` -- the final time for the integration

    - ``delta`` -- the step of the output

    - ``parameters`` -- the variables inside the function that should be treated
      as parameters

    - ``parameter_values`` -- the values of the parameters for the particular
      initial value problem

    - ``dig`` -- the number of digits of precision that will be used in the integration

    - ``tolrel`` -- the relative tolerance

    - ``tolabs`` -- the absolute tolerance

    - ``output`` -- the name of the file that the compiled integrator will write to

    This function creates two files, integrator and driver, that can be used
    later with the tides library ([TIDES]_).


    TESTS::

        sage: from tempfile import mkdtemp
        sage: from sage.interfaces.tides import genfiles_mpfr
        sage: import os
        sage: import shutil
        sage: from sage.misc.temporary_file import tmp_dir
        sage: tempdir = tmp_dir()
        sage: intfile = os.path.join(tempdir, 'integrator.c')
        sage: drfile = os.path.join(tempdir ,'driver.c')
        sage: var('t,x,y,X,Y')
        (t, x, y, X, Y)
        sage: f(t,x,y,X,Y)=[X, Y, -x/(x^2+y^2)^(3/2), -y/(x^2+y^2)^(3/2)]
        sage: genfiles_mpfr(intfile, drfile, f, [1,0, 0, 0.2], 0, 10, 0.1, output = 'out', dig = 50)
        sage: fileint = open(intfile)
        sage: l = fileint.readlines()
        sage: fileint.close()
        sage: l[5]
        '    #include "mp_tides.h"\n'
        sage: l[15]
        '\tstatic int PARAMETERS = 0;\n'
        sage: l[25]
        '\t\tmpfrts_var_t(itd, link[5], var[3], i);\n'
        sage: l[30]
        '\t\tmpfrts_pow_t_c(itd, link[2], "-1.500000000000000000000000000000000000000000000000000", link[3], i);\n'
        sage: l[35]
        '\n'
        sage: l[36]
        '    }\n'
        sage: l[37]
        '    write_mp_solution();\n'
        sage: filedr = open(drfile)
        sage: l = filedr.readlines()
        sage: filedr.close()
        sage: l[6]
        '    #include "mpfr.h"\n'
        sage: l[16]
        '    int nfun = 0;\n'
        sage: l[26]
        '\tmpfr_set_str(v[2], "0.000000000000000000000000000000000000000000000000000", 10, TIDES_RND);\n'
        sage: l[30]
        '\tmpfr_init2(tolabs, TIDES_PREC); \n'
        sage: l[34]
        '\tmpfr_init2(tini, TIDES_PREC); \n'
        sage: l[40]
        '\tmp_tides_delta(function_iteration, NULL, nvar, npar, nfun, v, p, tini, dt, nipt, tolrel, tolabs, NULL, fd);\n'
        sage: shutil.rmtree(tempdir)

    Check that issue :issue:`17179` is fixed (handle expressions like `\\pi`)::

        sage: from sage.interfaces.tides import genfiles_mpfr
        sage: import os
        sage: import shutil
        sage: from sage.misc.temporary_file import tmp_dir
        sage: tempdir = tmp_dir()
        sage: intfile = os.path.join(tempdir, 'integrator.c')
        sage: drfile = os.path.join(tempdir ,'driver.c')
        sage: var('t,x,y,X,Y')
        (t, x, y, X, Y)
        sage: f(t,x,y,X,Y)=[X, Y, -x/(x^2+y^2)^(3/2), -y/(x^2+y^2)^(3/2)]
        sage: genfiles_mpfr(intfile, drfile, f, [pi, 0, 0, 0.2], 0, 10, 0.1, output = 'out', dig = 50)
        sage: fileint = open(intfile)
        sage: l = fileint.readlines()
        sage: fileint.close()
        sage: l[30]
        '\t\tmpfrts_pow_t_c(itd, link[2], "-1.500000000000000000000000000000000000000000000000000", link[3], i);\n'
        sage: filedr = open(drfile)
        sage: l = filedr.readlines()
        sage: filedr.close()
        sage: l[24]
        '\tmpfr_set_str(v[0], "3.141592653589793238462643383279502884197169399375101", 10, TIDES_RND);\n'
        sage: shutil.rmtree(tempdir)
    """
    if parameters is None:
        parameters = []
    if parameter_values is None:
        parameter_values = []
    RR = RealField(ceil(dig * 3.3219))
    l1, l2 = subexpressions_list(f, parameters)
    remove_repeated(l1, l2)
    remove_constants(l1, l2)
    l3 = []
    var = f[0].arguments()
    l0 = [str(l) for l in l1]
    lv = [str(v) for v in var]
    lp = [str(p) for p in parameters]
    for i in l2:
        oper = i[0]
        if oper in ["log", "exp", "sin", "cos", "atan", "asin", "acos"]:
            a = i[1]
            if str(a) in lv:
                l3.append((oper, 'var[{}]'.format(lv.index(str(a)))))
            elif str(a) in lp:
                l3.append((oper, 'par[{}]'.format(lp.index(str(a)))))
            else:
                l3.append((oper, 'link[{}]'.format(l0.index(str(a)))))

        else:
            a = i[1]
            b = i[2]
            sa = str(a)
            sb = str(b)
            consta = False
            constb = False

            if sa in lv:
                aa = 'var[{}]'.format(lv.index(sa))
            elif sa in l0:
                aa = 'link[{}]'.format(l0.index(sa))
            elif sa in lp:
                aa = 'par[{}]'.format(lp.index(sa))
            else:
                consta = True
                aa = RR(a).str()
            if sb in lv:
                bb = 'var[{}]'.format(lv.index(sb))
            elif sb in l0:
                bb = 'link[{}]'.format(l0.index(sb))
            elif sb in lp:
                bb = 'par[{}]'.format(lp.index(sb))
            else:
                constb = True
                bb = RR(b).str()
            if consta:
                oper += '_c'
                if not oper == 'div':
                    bb, aa = aa,bb
            elif constb:
                oper += '_c'
            l3.append((oper, aa, bb))

    n = len(var)
    code = []

    l0 = lv + l0
    indices = [l0.index(str(i(*var))) + n for i in f]
    for i in range(1, n):
        aux = indices[i - 1] - n
        if aux < n:
            code.append(f'mpfrts_var_t(itd, var[{aux}], var[{i}], i);')
        else:
            code.append(f'mpfrts_var_t(itd, link[{aux-n}], var[{i}], i);')

    for i in range(len(l3)):
        el = l3[i]
        string = "mpfrts_"
        if el[0] == 'add':
            string += 'add_t(itd, ' + el[1] + ', ' + el[2] + f', link[{i}], i);'
        elif el[0] == 'add_c':
            string += 'add_t_c(itd, "' + el[2] + '", ' + el[1] + f', link[{i}], i);'
        elif el[0] == 'mul':
            string += 'mul_t(itd, ' + el[1] + ', ' + el[2] + f', link[{i}], i);'
        elif el[0] == 'mul_c':
            string += 'mul_t_c(itd, "' + el[2] + '", ' + el[1] + f', link[{i}], i);'
        elif el[0] == 'pow_c':
            string += 'pow_t_c(itd, ' + el[1] + ', "' + el[2] + f'", link[{i}], i);'
        elif el[0] == 'div':
            string += 'div_t(itd, ' + el[2] + ', ' + el[1] + f', link[{i}], i);'
        elif el[0] == 'div_c':
            string += 'div_t_cv(itd, "' + el[2] + '", ' + el[1] + f', link[{i}], i);'
        elif el[0] == 'log':
            string += 'log_t(itd, ' + el[1] + f', link[{i}], i);'
        elif el[0] == 'exp':
            string += 'exp_t(itd, ' + el[1] + f', link[{i}], i);'
        elif el[0] == 'sin':
            string += 'sin_t(itd, ' + el[1] + f', link[{i+1}], link[{i}], i);'
        elif el[0] == 'cos':
            string += 'cos_t(itd, ' + el[1] + f', link[{i-1}], link[{i}], i);'
        elif el[0] == 'atan':
            indarg = l0.index(str(1+l2[i][1]**2))-n
            string += 'atan_t(itd, ' + el[1] + f', link[{indarg}], link[{i}], i);'
        elif el[0] == 'asin':
            indarg = l0.index(str(sqrt(1-l2[i][1]**2)))-n
            string += 'asin_t(itd, ' + el[1] + f', link[{indarg}], link[{i}], i);'
        elif el[0] == 'acos':
            indarg = l0.index(str(-sqrt(1-l2[i][1]**2)))-n
            string += 'acos_t(itd, ' + el[1] + f', link[{indarg}], link[{i}], i);'
        code.append(string)

    VAR = n - 1
    PAR = len(parameters)
    TT = len(code) + 1 - VAR

    auxstring1 = """
    /****************************************************************************
    This file has been created by Sage for its use with TIDES
    *****************************************************************************/

    #include "mp_tides.h"

    long  function_iteration(iteration_data *itd, mpfr_t t, mpfr_t v[], mpfr_t p[], int ORDER, mpfr_t *cvfd)
    {

    int i;
    int NCONST = 0;
    mpfr_t ct[0];
    """

    auxstring2 = """
    }
    write_mp_solution();
    clear_vpl();
    clear_cts();
    return NUM_COLUMNS;
}
    """

    with open(integrator, 'a') as outfile:
        outfile.write(auxstring1)

        outfile.write(f"\n\tstatic int VARIABLES = {VAR};\n")
        outfile.write(f"\tstatic int PARAMETERS = {PAR};\n")
        outfile.write(f"\tstatic int LINKS = {TT};\n")
        outfile.write('\tstatic int   FUNCTIONS        = 0;\n')
        outfile.write('\tstatic int   POS_FUNCTIONS[1] = {0};\n')
        outfile.write('\n\tinitialize_mp_case();\n')
        outfile.write('\n\tfor(i=0;  i<=ORDER; i++) {\n')
        for i in code:
            outfile.write('\t\t' + i + '\n')

        outfile.write(auxstring2)

    npar = len(parameter_values)
    outfile = open(driver, 'a')

    auxstring3 = """
    /****************************************************************************
    Driver file of the mp_tides program
    This file has been created automatically by Sage
    *****************************************************************************/

    #include "mpfr.h"
    #include "mp_tides.h"
    long  function_iteration(iteration_data *itd, mpfr_t t, mpfr_t v[], mpfr_t p[], int ORDER, mpfr_t *cvfd);

    int main() {

        int i;



    int nfun = 0;
    """
    outfile.write(auxstring3)
    outfile.write('\tset_precision_digits({});'.format(dig))
    outfile.write('\n\tint npar = {};\n'.format(npar))
    outfile.write('\tmpfr_t p[npar];\n')
    outfile.write('\tfor(i=0; i<npar; i++) mpfr_init2(p[i], TIDES_PREC);\n')

    for i in range(npar):
        outfile.write('\tmpfr_set_str(p[{}], "{}", 10, TIDES_RND);\n'.format(i,RR(parameter_values[i]).str()))
    outfile.write('\tint nvar = {};\n\tmpfr_t v[nvar];\n'.format(VAR))
    outfile.write('\tfor(i=0; i<nvar; i++) mpfr_init2(v[i], TIDES_PREC);\n')
    for i in range(len(ics)):
        outfile.write('\tmpfr_set_str(v[{}], "{}", 10, TIDES_RND);\n'.format(i,RR(ics[i]).str()))
    outfile.write('\tmpfr_t tolrel, tolabs;\n')
    outfile.write('\tmpfr_init2(tolrel, TIDES_PREC); \n')
    outfile.write('\tmpfr_init2(tolabs, TIDES_PREC); \n')
    outfile.write('\tmpfr_set_str(tolrel, "{}", 10, TIDES_RND);\n'.format(RR(tolrel).str()))
    outfile.write('\tmpfr_set_str(tolabs, "{}", 10, TIDES_RND);\n'.format(RR(tolabs).str()))

    outfile.write('\tmpfr_t tini, dt; \n')
    outfile.write('\tmpfr_init2(tini, TIDES_PREC); \n')
    outfile.write('\tmpfr_init2(dt, TIDES_PREC); \n')

    outfile.write('\tmpfr_set_str(tini, "{}", 10, TIDES_RND);;\n'.format(RR(initial).str()))
    outfile.write('\tmpfr_set_str(dt, "{}", 10, TIDES_RND);\n'.format(RR(delta).str()))
    outfile.write('\tint nipt = {};\n'.format(floor((final-initial)/delta)))
    outfile.write('\tFILE* fd = fopen("' + output + '", "w");\n')
    outfile.write('\tmp_tides_delta(function_iteration, NULL, nvar, npar, nfun, v, p, tini, dt, nipt, tolrel, tolabs, NULL, fd);\n')
    outfile.write('\tfclose(fd);\n\treturn 0;\n}')
    outfile.close()
