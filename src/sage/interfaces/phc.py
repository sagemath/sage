r"""
Interface to PHC.

PHC computes numerical information about systems of polynomials over
the complex numbers.

PHC implements polynomial homotopy methods to exploit structure in
order to better approximate all isolated solutions. The package also
includes extra tools to handle positive dimensional solution
components.

AUTHORS:

- PHC was written by J. Verschelde, R. Cools, and many others (?)

- William Stein and Kelly ?? -- first version of interface to PHC

- Marshall Hampton -- second version of interface to PHC

- Marshall Hampton and Alex Jokela -- third version, path tracking
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2008 Marshall Hampton <hamptonio@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import os
import re
import pexpect
from random import random

from sage.misc.temporary_file import tmp_filename
from sage.rings.real_mpfr import RR
from sage.rings.cc import CC
from sage.rings.integer import Integer
from sage.misc.lazy_import import lazy_import
lazy_import("sage.plot.line", "line")
lazy_import("sage.plot.point", "point")


def get_solution_dicts(output_file_contents, input_ring, get_failures=True):
    """
    Return a list of dictionaries of variable:value (key:value)
    pairs.  Only used internally; see the solution_dict function in
    the PHC_Object class definition for details.

    INPUT:

    - ``output_file_contents`` -- phc solution output as a string
    - ``input_ring`` -- a PolynomialRing that variable names can be coerced into

    OUTPUT: list of dictionaries of solutions

    EXAMPLES::

        sage: from sage.interfaces.phc import *
        sage: R2.<x1,x2> = PolynomialRing(QQ,2)
        sage: test_sys = [(x1-1)^5-x2, (x2-1)^5-1]
        sage: sol = phc.blackbox(test_sys, R2)             # optional -- phc
        sage: test = get_solution_dicts(sol.output_file_contents,R2)  # optional -- phc
        sage: str(sum([q[x1].real() for q in test]))[0:4]  # optional -- phc
        '25.0'
    """
    output_list = output_file_contents.splitlines()
    solution_dicts = []
    for solution_line in range(len(output_list) - 1, -1, -1):
        if output_list[solution_line].find('THE SOLUTIONS') == 0:
            break
    try:
        var_number = int(output_list[solution_line + 2].split(' ')[1])
        # sol_number = int(output_list[solution_line+2].split(' ')[0])
    except IndexError:
        var_number = int(output_list[solution_line + 1].split(' ')[1])
        # sol_number = int(output_list[solution_line+1].split(' ')[0])
    for i in range(solution_line + 1, len(output_list)):
        if output_list[i].count('the solution for t') == 1:
            if output_list[i - 3].count('success') > 0 or get_failures:
                temp_dict = {}
                for j in range(1, var_number + 1):
                    rawsplit = output_list[i + j].split(': ')[1].split(' ')
                    for extras in range(rawsplit.count('')):
                        rawsplit.remove('')
                    temp_var = output_list[i + j].split(': ')[0].replace(' ', '')
                    temp_dict[input_ring(temp_var)] = CC(rawsplit[0],
                                                         rawsplit[1])
                solution_dicts.append(temp_dict)
    return solution_dicts


def get_classified_solution_dicts(output_file_contents, input_ring, get_failures=True):
    """
    Return a dictionary of lists of dictionaries of variable:value (key:value)
    pairs.  Only used internally; see the classified_solution_dict function in
    the PHC_Object class definition for details.

    INPUT:

    - ``output_file_contents`` -- phc solution output as a string
    - ``input_ring`` -- a PolynomialRing that variable names can be coerced into

    OUTPUT: a dictionary of lists if dictionaries of solutions, classifies by type

    EXAMPLES::

        sage: from sage.interfaces.phc import *
        sage: R2.<x1,x2> = PolynomialRing(QQ,2)
        sage: test_sys = [(x1-2)^5-x2, (x2-1)^5-1]
        sage: sol = phc.blackbox(test_sys, R2)          # optional -- phc
        sage: sol_classes = get_classified_solution_dicts(sol.output_file_contents,R2)  # optional -- phc
        sage: len(sol_classes['real'])            # optional -- phc
        1
    """
    output_list = output_file_contents.splitlines()
    solution_dicts = {}
    solution_types = ['complex', 'real', 'failure']
    for sol_type in solution_types:
        solution_dicts[sol_type] = []
    for solution_line in range(len(output_list) - 1, -1, -1):
        if output_list[solution_line].find('THE SOLUTIONS') == 0:
            break
    var_number = int(output_list[solution_line + 2].split(' ')[1])
    # sol_number = int(output_list[solution_line+2].split(' ')[0])
    for i in range(solution_line + 1, len(output_list)):
        if output_list[i].count('the solution for t') == 1:
            phc_type = output_list[i + var_number + 1].split(' = ')[-1]
            if phc_type.find('complex') != -1:
                phc_type = 'complex'
            elif phc_type.find('real') != -1:
                phc_type = 'real'
            else:
                phc_type = 'failure'
            temp_dict = {}
            for j in range(1, var_number + 1):
                rawsplit = output_list[i + j].split(': ')[1].split(' ')
                for extras in range(rawsplit.count('')):
                    rawsplit.remove('')
                temp_var = output_list[i + j].split(': ')[0].replace(' ', '')
                if phc_type == 'real':
                    temp_dict[input_ring(temp_var)] = RR(rawsplit[0])
                else:
                    temp_dict[input_ring(temp_var)] = CC(rawsplit[0],
                                                         rawsplit[1])
            solution_dicts[phc_type].append(temp_dict)
    return solution_dicts


def get_variable_list(output_file_contents):
    """
    Return the variables, as strings, in the order in which PHCpack has processed them.

    EXAMPLES::

        sage: from sage.interfaces.phc import *
        sage: R2.<x1,x2> = PolynomialRing(QQ,2)
        sage: test_sys = [(x1-2)^5-x2, (x2-1)^5-1]
        sage: sol = phc.blackbox(test_sys, R2)             # optional -- phc
        sage: get_variable_list(sol.output_file_contents)  # optional -- phc
        ['x1', 'x2']
    """
    output_list = output_file_contents.splitlines()
    for solution_line in range(len(output_list) - 1, -1, -1):
        if output_list[solution_line].find('THE SOLUTIONS') == 0:
            break
    var_number = int(output_list[solution_line + 2].split(' ')[1])
    varlist = []
    for var_ind in range(var_number):
        var = output_list[solution_line + 8 + var_ind].split(' ')[1]
        varlist.append(var)
    return varlist


class PHC_Object:

    def __init__(self, output_file_contents, input_ring):
        """
        A container for data from the PHCpack program - lists of float
        solutions, etc.  Currently the file contents are kept as a string;
        for really large outputs this would be bad.

        INPUT:

        - ``output_file_contents`` -- the string output of PHCpack
        - ``input_ring`` -- for coercion of the variables into the desired ring

        EXAMPLES::

            sage: from sage.interfaces.phc import phc
            sage: R2.<x,y> = PolynomialRing(QQ,2)
            sage: start_sys = [(x-1)^2+(y-1)-1, x^2+y^2-1]
            sage: sol = phc.blackbox(start_sys, R2)  # optional -- phc
            sage: str(sum([x[0] for x in sol.solutions()]).real())[0:3]  # optional -- phc
            '2.0'
        """
        self.output_file_contents = output_file_contents
        self.input_ring = input_ring

    def save_as_start(self, start_filename=None, sol_filter=''):
        """
        Save a solution as a phcpack start file.  The usual output is
        just as a string, but it can be saved to a file as well.  Even
        if saved to a file, it still returns the output string.

        EXAMPLES::

            sage: from sage.interfaces.phc import phc
            sage: R2.<x,y> = PolynomialRing(QQ,2)
            sage: start_sys = [x^3-y^2,y^5-1]
            sage: sol = phc.blackbox(start_sys, R2)  # optional -- phc
            sage: start_save = sol.save_as_start()   # optional -- phc
            sage: end_sys = [x^7-2,y^5-x^2]          # optional -- phc
            sage: sol = phc.start_from(start_save, end_sys, R2)  # optional -- phc
            sage: len(sol.solutions())               # optional -- phc
            15
        """
        start_data = ''
        output_list = self.output_file_contents.splitlines()
        for a_line in output_list:
            if a_line.find('ROOT COUNTS') != -1 or a_line.find('START SOLUTIONS') != -1:
                break
            else:
                start_data += a_line + '\n'
        for index in range(len(output_list) - 1, 0, -1):
            a_line = output_list[index]
            if a_line.find('THE SOLUTIONS') != -1:
                found_solutions = index
                break
        start_data += output_list[found_solutions] + '\n\n'
        try:
            var_number = int(output_list[found_solutions + 1].split(' ')[1])
        except Exception:
            # bad error handling
            var_number = int(output_list[found_solutions + 2].split(' ')[1])
        sol_count = 0
        sol_data = ''
        for i in range(found_solutions + 2, len(output_list)):
            if output_list[i].count('the solution for t') == 1 and output_list[i + 1 + var_number].find(sol_filter) != -1:
                phc_type = output_list[i + var_number + 1].split(' = ')[-1]
                if phc_type.find('no solution') == -1:
                    sol_count += 1
                    for ind2 in range(i - 3, i + var_number + 2):
                        sol_data += output_list[ind2] + '\n'
        jan_bar = '=' * 75 + '\n'
        sol_data += jan_bar
        start_data += str(sol_count) + ' ' + str(var_number) + '\n'
        start_data += jan_bar + sol_data
        if start_filename is not None:
            with open(start_filename, 'w') as start_file:
                start_file.write(start_data)
        return start_data

    def classified_solution_dicts(self):
        """
        Return a dictionary of lists of dictionaries of solutions.

        Its not as crazy as it sounds; the keys are the types of solutions as
        classified by phcpack: regular vs. singular, complex vs. real

        INPUT:

        - None

        OUTPUT: a dictionary of lists of dictionaries of solutions

        EXAMPLES::

            sage: from sage.interfaces.phc import phc
            sage: R.<x,y> = PolynomialRing(CC,2)
            sage: p_sys = [x^10-y,y^2-1]
            sage: sol = phc.blackbox(p_sys,R)         # optional -- phc
            sage: classifieds = sol.classified_solution_dicts()          # optional -- phc
            sage: str(sum([q[y] for q in classifieds['real']]))[0:3]     # optional -- phc
            '2.0'
        """
        try:
            return self.__classified_sols
        except AttributeError:
            pass
        classified_sols = get_classified_solution_dicts(self.output_file_contents, self.input_ring)
        self.__classified_sols = classified_sols
        return classified_sols

    def solution_dicts(self, get_failures=False):
        """
        Return a list of solutions in dictionary form: variable:value.

        INPUT:

        - ``self`` -- for access to self_out_file_contents, the string
          of raw PHCpack output
        - ``get_failures`` -- boolean (default: ``False``); the default
          is to not process failed homotopies.  These either lie on
          positive-dimensional components or at infinity.

        OUTPUT:

        - solution_dicts: a list of dictionaries.  Each dictionary
          element is of the form variable:value, where the variable
          is an element of the input_ring, and the value is in
          ComplexField.

        EXAMPLES::

            sage: from sage.interfaces.phc import *
            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: fs = [x^2-1,y^2-x,z^2-y]
            sage: sol = phc.blackbox(fs,R)          # optional -- phc
            sage: s_list = sol.solution_dicts()     # optional -- phc
            sage: s_list.sort()                     # optional -- phc
            sage: s_list[0]                         # optional -- phc
            {y: 1.00000000000000, z: -1.00000000000000, x: 1.00000000000000}
        """
        try:
            return self.__solution_dicts
        except AttributeError:
            pass
        solution_dicts = get_solution_dicts(self.output_file_contents, self.input_ring, get_failures=get_failures)
        self.__solution_dicts = solution_dicts
        return solution_dicts

    def solutions(self, get_failures=False):
        """
        Return a list of solutions in the ComplexField.

        Use the variable_list function to get the order of variables used by
        PHCpack, which is usually different than the term order of the
        input_ring.

        INPUT:

        - ``self`` -- for access to self_out_file_contents, the string
          of raw PHCpack output
        - ``get_failures`` -- boolean (default: ``False``); the default
          is to not process failed homotopies.  These either lie on
          positive-dimensional components or at infinity.

        OUTPUT: solutions: a list of lists of ComplexField-valued solutions

        EXAMPLES::

            sage: from sage.interfaces.phc import *
            sage: R2.<x1,x2> = PolynomialRing(QQ,2)
            sage: test_sys = [x1^5-x1*x2^2-1, x2^5-x1*x2-1]
            sage: sol = phc.blackbox(test_sys, R2)          # optional -- phc
            sage: len(sol.solutions())                      # optional -- phc
            25
        """
        try:
            return self.__solutions
        except AttributeError:
            pass
        solution_dicts = get_solution_dicts(self.output_file_contents, self.input_ring, get_failures=get_failures)
        self.__solution_dicts = solution_dicts
        solutions = [sol_dict.values() for sol_dict in solution_dicts]
        self.__solutions = solutions
        return solutions

    def variable_list(self):
        """
        Return the variables, as strings, in the order in which
        PHCpack has processed them.

        EXAMPLES::

            sage: from sage.interfaces.phc import *
            sage: R2.<x1,x2> = PolynomialRing(QQ,2)
            sage: test_sys = [x1^5-x1*x2^2-1, x2^5-x1*x2-1]
            sage: sol = phc.blackbox(test_sys, R2)          # optional -- phc
            sage: sol.variable_list()                       # optional -- phc
            ['x1', 'x2']
        """
        try:
            return self.__var_list
        except AttributeError:
            pass
        var_list = get_variable_list(self.output_file_contents)
        self.__var_list = var_list
        return var_list


class PHC:
    """
    A class to interface with PHCpack, for computing numerical
    homotopies and root counts.

    EXAMPLES::

        sage: from sage.interfaces.phc import phc
        sage: R.<x,y> = PolynomialRing(CDF,2)
        sage: testsys = [x^2 + 1, x*y - 1]
        sage: phc.mixed_volume(testsys)        # optional -- phc
        2
        sage: v = phc.blackbox(testsys, R)     # optional -- phc
        sage: sols = v.solutions()             # optional -- phc
        sage: sols.sort()                      # optional -- phc
        sage: sols                             # optional -- phc
        [[-1.00000000000000*I, 1.00000000000000*I], [1.00000000000000*I, -1.00000000000000*I]]
        sage: sol_dict = v.solution_dicts()    # optional -- phc
        sage: x_sols_from_dict = [d[x] for d in sol_dict]    # optional -- phc
        sage: x_sols_from_dict.sort(); x_sols_from_dict      # optional -- phc
        [-1.00000000000000*I, 1.00000000000000*I]
        sage: residuals = [[test_equation.change_ring(CDF).subs(sol) for test_equation in testsys] for sol in v.solution_dicts()]  # optional -- phc
        sage: residuals                        # optional -- phc
        [[0, 0], [0, 0]]
    """

    def _output_from_command_list(self, command_list, polys, verbose=False):
        """
        A pexpect interface to phcpack, given a command list for
        interactive dialogs. The input file is supplied from the
        polynomial list, output file is also supplied.  This is
        only used as a building block for the interface.

        INPUT:

        - ``command_list`` -- list of commands to phc
        - ``polys`` -- a polynomial system as a list of polynomials

        OUTPUT: an output string from phc

        EXAMPLES::

            sage: from sage.interfaces.phc import *
            sage: R2.<x,y> = PolynomialRing(QQ,2)
            sage: start_sys = [(x-1)^2+(y-1)-1, x^2+y^2-1]  # optional -- phc
            sage: a = phc._output_from_command_list(['phc -m','4','n','n','n'], start_sys)  # optional -- phc
        """
        # Get temporary file names (these will be in SAGE_HOME/.sage/tmp/pid)
        input_filename = tmp_filename()
        output_filename = tmp_filename()

        # Get the input polynomial text
        input = self._input_file(polys)
        if verbose:
            print("Writing the input file to %s" % input_filename)
        with open(input_filename, 'w') as file:
            file.write(input)

        if verbose:
            print("The following file will be the input polynomial file to phc.")
            print(input)

        # Create a phc process
        child_phc = pexpect.spawn(command_list[0])
        # feed it the commands
        child_phc.sendline('y')
        child_phc.sendline(input_filename)
        child_phc.sendline(output_filename)
        for command_string in command_list[1:]:
            if verbose:
                print(command_string)
            child_phc.sendline(command_string)
        child_phc.expect('results')
        read_stuff = child_phc.read()
        if verbose:
            print(read_stuff)
        child_phc.close()
        if not os.path.exists(output_filename):
            raise RuntimeError("The output file does not exist; something went wrong running phc.")

        # Delete the input file
        os.unlink(input_filename)

        # Return the output filename
        return output_filename

    def _input_file(self, polys):
        r"""
        This is used internally to implement the PHC interface.

        INPUT:

        - ``polys`` -- list of polynomials in a Sage polynomial ring
          over a field that embeds into the complex numbers

        OUTPUT: a PHC input file (as a text string) that describes these
        polynomials

        EXAMPLES::

            sage: from sage.interfaces.phc import *
            sage: R2.<x,y> = PolynomialRing(QQ,2)
            sage: start_sys = [(x-1)^2+(y-1)-1, x^2+y^2-1]
            sage: phc._input_file(start_sys)         # optional -- phc
            '2\nx^2 - 2*x + y - 1;\nx^2 + y^2 - 1;\n'
        """
        if not isinstance(polys, (list, tuple)):
            raise TypeError('polys must be a list or tuple')
        s = '%s\n' % len(polys)
        for f in polys:
            s += f._repr_() + ';\n'     # note the semicolon *terminators*
        return s

    def _parse_path_file(self, input_filename, verbose=False):
        """
        Take a phpack output file containing path tracking information
        and parses it into a list of lists of dictionaries - i.e. a
        list of solutions paths, where each solution path is a list of
        dictionaries of variable and homotopy parameter values.

        INPUT:

        - ``input_filename`` -- file must have path-tracking information

        OUTPUT: list of lists of dictionaries, described above

        EXAMPLES::

            sage: from sage.interfaces.phc import *
            sage: R2.<x,y> = PolynomialRing(QQ,2)
            sage: start_sys = [x^5-y^2,y^5-1]
            sage: sol = phc.blackbox(start_sys, R2)    # optional -- phc
            sage: start_save = sol.save_as_start()     # optional -- phc
            sage: end_sys = [x^5-2,y^5-x^2]            # optional -- phc
            sage: path_track_filename = phc._path_track_file(start_save, end_sys, R2, c_skew = .001)  # optional -- phc
            sage: sol_paths = phc._parse_path_file(path_track_filename)   # optional -- phc
            sage: len(sol_paths)  # optional -- phc
            25
        """
        if not os.path.exists(input_filename):
            raise RuntimeError("The file containing output from phc (" + input_filename + ") cannot be found")

        fh = open(input_filename)
        line_idx = 0
        begin = 0
        count = 0
        solutions_dicts = []
        steps_dicts = []

        # regular expressions for matching certain output types
        var_cnt_regex = re.compile('^ +([0-9]+)')
        output_regex = re.compile('^OUTPUT INFORMATION DURING')
        t_regex = re.compile(r'(^t +: +(-{0,1}[0-9]+\.[0-9]+E[-+][0-9]+) +(-{0,1}[0-9]+\.[0-9]+E[-+][0-9]+)$)', re.IGNORECASE)
        sols_regex = re.compile(r'(^ *(([a-z]|[0-9])+) +: +(-?[0-9]+\.[0-9]+E[-+][0-9]+) +(-?[0-9]+\.[0-9]+E[-+][0-9]+)$)', re.IGNORECASE)
        complete_regex = re.compile('^TIMING INFORMATION')

        breakfast = False
        a_line = fh.readline()
        end_test = ''
        while a_line:
            # processing....
            a_line = a_line.replace("\n", '')
            if line_idx == 0:
                m = var_cnt_regex.match(a_line)
                if m:
                    count = Integer(m.group(1))
            if count > 0:
                m = output_regex.match(a_line)
                if m:
                    begin = 1
                if begin:
                    m = t_regex.match(a_line)
                    if m:
                        # put the t-values into a dict
                        # m.group(2) contains the real val
                        # m.group(3) contains the imaginary val
                        # fh_w.write( "T=> G1(" + m.group(2) + '),G2(' + m.group(3) + ")\n")
                        # read off two lines - this should be 'm' and 'the solution for t :'
                        a_line = fh.readline()
                        end_test = a_line  # store this to check for end of solution
                        a_line = fh.readline()
                        t_val = CC(m.group(2), m.group(3))
                        temp_dict = {}
                        temp_dict["t"] = t_val
                        for i in range(count):
                            a_line = fh.readline()
                            m = sols_regex.match(a_line)
                            if m:
                                # m.group(2) contains our var name
                                # m.group(4) contains our real val
                                # m.group(5) contains our imaginary val
                                temp_dict[m.group(2)] = CC(m.group(4),
                                                           m.group(5))
                        steps_dicts.append(temp_dict)
                    # check if its the end of a solution
                    if end_test.find('Length of path') != -1:
                        if verbose:
                            print("recording sol")
                        if steps_dicts != []:
                            solutions_dicts.append(steps_dicts)
                        steps_dicts = []
                    m = complete_regex.match(a_line)
                    if m:
                        breakfast = True
            if breakfast:
                break
            line_idx += 1
            a_line = fh.readline()
        fh.close()
        return solutions_dicts

    def _path_track_file(self, start_filename_or_string, polys, input_ring, c_skew=0.001, verbose=False):
        """
        Return the filename which contains path tracking output.

        EXAMPLES::

            sage: from sage.interfaces.phc import *
            sage: R2.<x,y> = PolynomialRing(QQ,2)
            sage: start_sys = [x^6-y^2,y^5-1]
            sage: sol = phc.blackbox(start_sys, R2)        # optional -- phc
            sage: start_save = sol.save_as_start()         # optional -- phc
            sage: end_sys = [x^7-2,y^5-x^2]                # optional -- phc
            sage: path_track_filename = phc._path_track_file(start_save, end_sys, R2, c_skew = .001)  # optional -- phc
            sage: sol_paths = phc._parse_path_file(path_track_filename)        # optional -- phc
            sage: len(sol_paths)        # optional -- phc
            30
        """
        # Probably unnecessarily redundant from the start_from function
        if start_filename_or_string.find('THE SOLUTIONS') != -1:
            start_filename = tmp_filename()
            with open(start_filename, 'w') as start_file:
                start_file.write(start_filename_or_string)
        elif os.path.exists(start_filename_or_string):
            start_filename = start_filename_or_string
        else:
            raise RuntimeError("There is something wrong with your start string or filename")

        return self._output_from_command_list(['phc', '0', '0', 'A', start_filename, 'y', '1', '0', 'n', 'k', '2', 'a', '1', str(c_skew), '0', '0', '2'], polys, verbose=verbose)

    def path_track(self, start_sys, end_sys, input_ring, c_skew=.001, saved_start=None):
        """
        This function computes homotopy paths between the solutions of
        ``start_sys`` and ``end_sys``.

        INPUT:

        - ``start_sys`` -- a square polynomial system, given as a list of
          polynomials
        - ``end_sys`` -- same type as ``start_sys``
        - ``input_ring`` -- for coercion of the variables into the desired ring
        - ``c_skew`` -- (optional) the imaginary part of homotopy multiplier;
          nonzero values are often necessary to avoid intermediate path
          collisions
        - ``saved_start`` -- (optional) a phc output file; if not given, start
          system solutions are computed via the ``phc.blackbox`` function

        OUTPUT: list of paths as dictionaries, with the keys variables and
        t-values on the path

        EXAMPLES::

            sage: from sage.interfaces.phc import *
            sage: R2.<x,y> = PolynomialRing(QQ,2)
            sage: start_sys = [x^6-y^2,y^5-1]
            sage: sol = phc.blackbox(start_sys, R2)        # optional -- phc
            sage: start_save = sol.save_as_start()         # optional -- phc
            sage: end_sys = [x^7-2,y^5-x^2]                # optional -- phc
            sage: sol_paths = phc.path_track(start_sys, end_sys, R2, saved_start = start_save)  # optional -- phc
            sage: len(sol_paths)        # optional -- phc
            30
        """
        if not saved_start:
            sol = phc.blackbox(start_sys, input_ring)
            saved_start = sol.save_as_start()
        path_track_filename = phc._path_track_file(saved_start, end_sys, input_ring=input_ring, c_skew=c_skew)
        sol_paths = phc._parse_path_file(path_track_filename)
        os.unlink(path_track_filename)
        return sol_paths

    def plot_paths_2d(self, start_sys, end_sys, input_ring, c_skew=.001, endpoints=True, saved_start=None, rand_colors=False):
        """
        Return a graphics object of solution paths in the complex plane.

        INPUT:

        - ``start_sys`` -- a square polynomial system, given as a list of
          polynomials
        - ``end_sys`` -- same type as start_sys
        - ``input_ring`` -- for coercion of the variables into the desired ring
        - ``c_skew`` -- (optional) the imaginary part of homotopy multiplier;
          nonzero values are often necessary to avoid intermediate path
          collisions
        - ``endpoints`` -- (optional) whether to draw in the ends of paths as
          points
        - ``saved_start`` -- (optional) a phc output file; if not given, start
          system solutions are computed via the ``phc.blackbox`` function

        OUTPUT: lines and points of solution paths

        EXAMPLES::

            sage: from sage.interfaces.phc import *
            sage: from sage.structure.sage_object import SageObject
            sage: R2.<x,y> = PolynomialRing(QQ,2)
            sage: start_sys = [x^5-y^2,y^5-1]
            sage: sol = phc.blackbox(start_sys, R2)    # optional -- phc
            sage: start_save = sol.save_as_start()     # optional -- phc
            sage: end_sys = [x^5-25,y^5-x^2]           # optional -- phc
            sage: testing = phc.plot_paths_2d(start_sys, end_sys, R2)  # optional -- phc
            sage: type(testing)                        # optional -- phc (normally use plot here)
            <class 'sage.plot.graphics.Graphics'>
        """
        paths = phc.path_track(start_sys, end_sys, input_ring, c_skew=c_skew, saved_start=saved_start)
        path_lines = []
        sol_pts = []
        if rand_colors:
            r_color = {}
            for a_var in input_ring.gens():
                var_name = str(a_var)
                r_color[var_name] = (random(), random(), random())
        for a_sol in paths:
            for a_var in input_ring.gens():
                var_name = str(a_var)
                temp_line = []
                for data in a_sol:
                    temp_line.append([data[var_name].real(), data[var_name].imag()])
                if rand_colors:
                    path_lines.append(line(temp_line, rgbcolor=r_color[var_name]))
                else:
                    path_lines.append(line(temp_line))
        if endpoints:
            sol_pts = []
            for a_sol in paths:
                for a_var in input_ring.gens():
                    var_name = str(a_var)
                    sol_pts.append(point([a_sol[0][var_name].real(), a_sol[0][var_name].imag()]))
                    sol_pts.append(point([a_sol[-1][var_name].real(), a_sol[-1][var_name].imag()]))
            return sum(sol_pts) + sum(path_lines)
        else:
            return sum(path_lines)

    def mixed_volume(self, polys, verbose=False):
        """
        Compute the mixed volume of the polynomial system given by the input polys.

        INPUT:

        - ``polys`` -- list of multivariate polynomials (elements of a multivariate
          polynomial ring).
        - ``verbose`` -- print lots of verbose information about what this function does

        OUTPUT: the mixed volume

        EXAMPLES::

            sage: from sage.interfaces.phc import *
            sage: R2.<x,y,z> = PolynomialRing(QQ,3)
            sage: test_sys = [(x+y+z)^2-1,x^2-x,y^2-1]
            sage: phc.mixed_volume(test_sys)                # optional -- phc
            4
        """
        output_filename = self._output_from_command_list(['phc -m', '4', 'n', 'n', 'n'], polys, verbose=verbose)

        with open(output_filename) as out:
            out.read()
        # All done
        out_lines = out.split('\n')
        for a_line in out_lines:
            # the two conditions below are necessary because of changes in output format
            if a_line.find('The mixed volume equals :') == 0 or a_line.find('common mixed volume :') == 0:
                if verbose:
                    print('found line: ' + a_line)
                mixed_vol = Integer(a_line.split(':')[1])
                break

        try:
            return mixed_vol
        except NameError:
            raise RuntimeError("Mixed volume not found in output; something went wrong running phc.")

    def start_from(self, start_filename_or_string, polys, input_ring, path_track_file=None, verbose=False):
        """
        This computes solutions starting from a phcpack solution file.

        INPUT:

        - ``start_filename_or_string`` -- the filename for a phcpack start system,
          or the contents of such a file as a string.  Variable names must match
          the inputring variables.  The value of the homotopy variable t should
          be 1, not 0.
        - ``polys`` -- list of multivariate polynomials (elements of a multivariate
          polynomial ring).
        - input_ring: for coercion of the variables into the desired ring.
        - path_track_file: whether to save path-tracking information
        - ``verbose`` -- print lots of verbose information about what this function does

        OUTPUT: a solution in the form of a PHCObject

        EXAMPLES::

            sage: from sage.interfaces.phc import *
            sage: R2.<x,y> = PolynomialRing(QQ,2)
            sage: start_sys = [x^6-y^2,y^5-1]
            sage: sol = phc.blackbox(start_sys, R2)        # optional -- phc
            sage: start_save = sol.save_as_start()         # optional -- phc
            sage: end_sys = [x^7-2,y^5-x^2]                # optional -- phc
            sage: sol = phc.start_from(start_save, end_sys, R2)  # optional -- phc
            sage: len(sol.solutions())                     # optional -- phc
            30
        """
        input_filename = tmp_filename()
        output_filename = tmp_filename()

        if start_filename_or_string.find('THE SOLUTIONS') != -1:
            start_filename = tmp_filename()
            with open(start_filename, 'w') as start_file:
                start_file.write(start_filename_or_string)
        elif os.path.exists(start_filename_or_string):
            start_filename = start_filename_or_string
        else:
            raise RuntimeError("There is something wrong with your start string or filename")

        # Get the input polynomial text
        input = self._input_file(polys)
        if verbose:
            print("Writing the input file to %s" % input_filename)
        with open(input_filename, 'w') as f:
            f.write(input)

        if verbose:
            print("The following file will be the input polynomial file to phc.")
            print(input)

        # Create a phc process
        child_phc = pexpect.spawn('phc')
        child_phc.sendline('y')
        child_phc.sendline(input_filename)
        child_phc.sendline(output_filename)
        child_phc.sendline('0')
        child_phc.sendline('0')
        child_phc.expect('Nonlinear Reduction')
        child_phc.sendline('A')
        child_phc.sendline(start_filename)
        child_phc.sendline('y')
        child_phc.sendline('1')
        child_phc.sendline('0')
        if verbose:
            phc_dialog = child_phc.read(size=40)
            print(phc_dialog)
        child_phc.sendline('n')
        child_phc.sendline('0')
        if verbose:
            child_phc.expect('CURRENT CONTINUATION')
            phc_dialog = child_phc.read(size=40)
            print(phc_dialog)
        child_phc.sendline('0')
        if path_track_file is None:
            child_phc.sendline('0')
        else:
            child_phc.sendline('2')
        child_phc.expect('results')
        dots = child_phc.read()
        if verbose:
            print("should be . : " + dots)

        # close down the process:
        child_phc.close()
        if not os.path.exists(output_filename):
            raise RuntimeError("The output file does not exist; something went wrong running phc.")

        # Read the output produced by PHC
        with open(output_filename) as f:
            out = f.read()

        # Delete the temporary files
        os.unlink(output_filename)
        os.unlink(input_filename)

        # All done
        return PHC_Object(out, input_ring)

    def blackbox(self, polys, input_ring, verbose=False):
        """
        Return as a string the result of running PHC with the given polynomials
        under blackbox mode (the '-b' option).

        INPUT:

        - ``polys`` -- list of multivariate polynomials (elements of a multivariate
          polynomial ring).
        - ``input_ring`` -- for coercion of the variables into the desired ring
        - ``verbose`` -- print lots of verbose information about what this function does

        OUTPUT: a PHC_Object object containing the phcpack output string

        EXAMPLES::

            sage: from sage.interfaces.phc import *
            sage: R2.<x,y> = PolynomialRing(QQ,2)
            sage: start_sys = [x^6-y^2,y^5-1]
            sage: sol = phc.blackbox(start_sys, R2)        # optional -- phc
            sage: len(sol.solutions())                     # optional -- phc
            30
        """

        # Get three temporary file names (these will be in SAGE_HOME/.sage/tmp/pid)
        input_filename = tmp_filename()
        output_filename = input_filename + ".phc"
        log_filename = tmp_filename()

        # Get the input polynomial text
        input = self._input_file(polys)
        if verbose:
            print("Writing the input file to %s" % input_filename)
        with open(input_filename, 'w') as f:
            f.write(input)

        if verbose:
            print("The following file will be the input polynomial file to phc.")
            print(input)

        # Create the phc command line>
        cmd = 'phc -b %s %s' % (input_filename, output_filename)

        if verbose:
            print("The phc command line is:")
            print(cmd)

        # Do it -- make the system call.
        e = os.system(cmd)

        # Was there an error?
        if e:
            from sage.misc.sage_ostools import have_program
            if not have_program('phc'):
                print(str(os.system('which phc')) + '  PHC needs to be installed and in your path')
                raise RuntimeError
            # todo -- why? etc.
            with open(log_filename) as f:
                msg = f.read()
            raise RuntimeError(msg + "\nError running phc.")

        if not os.path.exists(output_filename):
            raise RuntimeError("The output file does not exist; something went wrong running phc.")

        # Read the output produced by PHC
        with open(output_filename) as f:
            out = f.read()

        # All done
        return PHC_Object(out, input_ring)


################################

# The unique phc interface instance.
phc = PHC()
