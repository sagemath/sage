r"""
Interface to MATLAB

According to their website, MATLAB is "a high-level language and
interactive environment that enables you to perform computationally
intensive tasks faster than with traditional programming languages
such as C, C++, and Fortran."

The commands in this section only work if you have the "matlab"
interpreter installed and available in your PATH. It's not
necessary to install any special Sage packages.

EXAMPLES::

    sage: matlab.eval('2+2')                 # optional - matlab
    '\nans =\n\n     4\n'

::

    sage: a = matlab(10)                     # optional - matlab
    sage: a**10                              # optional - matlab
       1.0000e+10

AUTHORS:

- William Stein (2006-10-11)

Tutorial
--------

EXAMPLES::

    sage: # optional - matlab
    sage: matlab('4+10')
    14
    sage: matlab('date')
    18-Oct-2006
    sage: matlab('5*10 + 6')
    56
    sage: matlab('(6+6)/3')
    4
    sage: matlab('9')^2
    81
    sage: a = matlab(10); b = matlab(20); c = matlab(30)
    sage: avg = (a+b+c)/3 ; avg
    20
    sage: parent(avg)
    Matlab

::

    sage: # optional - matlab
    sage: my_scalar = matlab('3.1415')
    sage: my_scalar
    3.1415
    sage: my_vector1 = matlab('[1,5,7]')
    sage: my_vector1
    1     5     7
    sage: my_vector2 = matlab('[1;5;7]')
    sage: my_vector2
    1
    5
    7
    sage: my_vector1 * my_vector2
    75

::

    sage: # optional - matlab
    sage: row_vector1 = matlab('[1 2 3]')
    sage: row_vector2 = matlab('[3 2 1]')
    sage: matrix_from_row_vec = matlab('[%s; %s]'%(row_vector1.name(), row_vector2.name()))
    sage: matrix_from_row_vec
    1     2     3
    3     2     1

::

    sage: # optional - matlab
    sage: column_vector1 = matlab('[1;3]')
    sage: column_vector2 = matlab('[2;8]')
    sage: matrix_from_col_vec = matlab('[%s %s]'%(column_vector1.name(), column_vector2.name()))
    sage: matrix_from_col_vec
    1     2
    3     8

::

    sage: my_matrix = matlab('[8, 12, 19; 7, 3, 2; 12, 4, 23; 8, 1, 1]')    # optional - matlab
    sage: my_matrix                                      # optional - matlab
         8    12    19
         7     3     2
        12     4    23
         8     1     1

::

    sage: combined_matrix = matlab('[%s, %s]'%(my_matrix.name(), my_matrix.name()))                                        # optional - matlab
    sage: combined_matrix                               # optional - matlab
     8    12    19     8    12    19
     7     3     2     7     3     2
    12     4    23    12     4    23
     8     1     1     8     1     1

::

    sage: tm = matlab('0.5:2:10')                       # optional - matlab
    sage: tm                                            # optional - matlab
    0.5000    2.5000    4.5000    6.5000    8.5000

::

    sage: # optional - matlab
    sage: my_vector1 = matlab('[1,5,7]')
    sage: my_vector1(1)
    1
    sage: my_vector1(2)
    5
    sage: my_vector1(3)
    7

Matrix indexing works as follows::

    sage: my_matrix = matlab('[8, 12, 19; 7, 3, 2; 12, 4, 23; 8, 1, 1]')     # optional - matlab
    sage: my_matrix(3,2)                                # optional - matlab
    4

Setting using parenthesis cannot work (because of how the Python
language works). Use square brackets or the set function::

    sage: my_matrix = matlab('[8, 12, 19; 7, 3, 2; 12, 4, 23; 8, 1, 1]')    # optional - matlab
    sage: my_matrix.set(2,3, 1999)                          # optional - matlab
    sage: my_matrix                                         # optional - matlab
               8          12          19
               7           3        1999
              12           4          23
               8           1           1
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL).
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

import os

from .expect import Expect, ExpectElement
from sage.misc.instancedoc import instancedoc


class Matlab(Expect):
    """
    Interface to the Matlab interpreter.

    EXAMPLES::

        sage: # optional - matlab
        sage: a = matlab('[ 1, 1, 2; 3, 5, 8; 13, 21, 33 ]')
        sage: b = matlab('[ 1; 3; 13]')
        sage: c = a * b
        sage: print(c)
            30
           122
           505
    """
    def __init__(self, maxread=None, script_subdirectory=None,
                 logfile=None, server=None, server_tmpdir=None):
        Expect.__init__(self,
                        name='matlab',
                        prompt='>> ',
                        command="matlab -nodisplay",
                        server=server,
                        server_tmpdir=server_tmpdir,
                        script_subdirectory=script_subdirectory,
                        restart_on_ctrlc=False,
                        verbose_start=False,
                        logfile=logfile,
                        eval_using_file_cutoff=100)

    def __reduce__(self):
        return reduce_load_Matlab, tuple([])

    def _read_in_file_command(self, filename):
        """
        Return the command used to read in and execute a file in Matlab.

        EXAMPLES::

            sage: matlab._read_in_file_command('/tmp/matlab_file')
            "eval(fileread('/tmp/matlab_file'));"

        Here is an indirect doctest to check that it does indeed
        work::

            sage: m = identity_matrix(ZZ, 10)
            sage: sm = matlab.sage2matlab_matrix_string(m)
            sage: m = matlab(sm)  # optional - matlab
        """
        return "eval(fileread('{0}'));".format(filename)

    def _quit_string(self):
        return 'quit;'

    def _install_hints(self):
        return """
        You must obtain the program MATLAB in order to use MATLAB
        from Sage.   You can read all about MATLAB at
                  http://www.mathworks.com/

        You might have to buy MATLAB or get away with setting up a remote connection to a server running Maple. Type
   print(matlab._install_hints_ssh())
for hints on how to do that).
        """

    def _start(self):
        Expect._start(self)

    def whos(self):
        return self.eval('whos')

#    pdehaye/20070819: This is no obsolete, see Expect._get_tmpfile_from_server and Expect._send_tmpfile_to_server

#    def get_via_file(self, var_name):
#        t = self._temp_file(var_name)
#        self.eval('save -text "%s" %s'%(t,var_name))
#        r = open(t).read()
#        os.unlink(t)
#        return r.strip('\n')

#    def set_via_file(self, var_name, x):
#        t = self._temp_file(var_name)
#        open(t,'w').write(x)
#        print('load "%s" %s'%(t, var_name))
#        self.eval('load "%s" %s'%(t, var_name))
#        #os.unlink(t)

    def set(self, var, value):
        """
        Set the variable var to the given value.
        """
        cmd = '{0}={1};'.format(var, value)
        out = self.eval(cmd)
        if out.find("error") != -1:
            raise TypeError("Error executing code in Matlab\nCODE:\n\t{0}\nMatlab ERROR:\n\t{1}".format(cmd, out))

    def get(self, var):
        """
        Get the value of the variable var.

        EXAMPLES::

            sage: s = matlab.eval('a = 2') # optional - matlab
            sage: matlab.get('a')               # optional - matlab
            '     2'
        """
        s = self.eval('{0}'.format(var))
        return self.strip_answer(s)

    def strip_answer(self, s):
        r"""
        Return the string s with Matlab's answer prompt removed.

        EXAMPLES::

            sage: s = '\nans =\n\n     2\n'
            sage: matlab.strip_answer(s)
            '     2'
        """
        i = s.find('=')
        return s[i+1:].strip('\n')

    def console(self):
        matlab_console()

    def version(self):
        return matlab_version()[1:]

    def chdir(self, directory):
        """
        Change MATLAB's current working directory.

        EXAMPLES::

            sage: matlab.chdir('/')          # optional - matlab
            sage: matlab.pwd()               # optional - matlab
            /
        """
        self.eval("cd('{0}')".format(directory))

    def sage2matlab_matrix_string(self, A):
        """
        Return a matlab matrix from a Sage matrix.

        INPUT:

        - ``A`` -- Sage matrix with entries in the rationals or reals

        OUTPUT: string that evaluates to a Matlab matrix

        EXAMPLES::

            sage: M33 = MatrixSpace(QQ,3,3)
            sage: A = M33([1,2,3,4,5,6,7,8,0])
            sage: matlab.sage2matlab_matrix_string(A)   # optional - matlab
            '[1, 2, 3; 4, 5, 6; 7, 8, 0]'

        AUTHOR:

        - David Joyner and William Stein
        """
        return str(A.rows()).replace('), (', '; ').replace('(', '').replace(')', '')

    def _object_class(self):
        return MatlabElement


@instancedoc
class MatlabElement(ExpectElement):
    def __getitem__(self, n):
        raise RuntimeError("Use parenthesis for MATLAB matrices instead.")

    def _matrix_(self, R):
        r"""
        Return Sage matrix from this matlab element.

        EXAMPLES::

            sage: # optional - matlab
            sage: A = matlab('[1,2;3,4]')
            sage: matrix(ZZ, A)
            [1 2]
            [3 4]
            sage: A = matlab('[1,2;3,4.5]')
            sage: matrix(RR, A)
            [1.00000000000000 2.00000000000000]
            [3.00000000000000 4.50000000000000]

            sage: a = matlab('eye(50)')         # optional - matlab
            sage: matrix(RR, a)                 # optional - matlab
            50 x 50 dense matrix over Real Field with 53 bits of precision
        """
        from sage.matrix.constructor import matrix
        matlab = self.parent()
        entries = matlab.strip_answer(matlab.eval("mat2str({0})".format(self.name())))
        entries = entries.strip()[1:-1].replace(';', ' ')
        entries = [R(_) for _ in entries.split(' ')]
        nrows, ncols = map(int, str(self.size()).strip().split())
        m = matrix(R, nrows, ncols, entries)
        return m

    def set(self, i, j, x):
        P = self._check_valid()
        z = P(x)
        P.eval('{0}({1},{2}) = {3}'.format(self.name(), i, j, z.name()))


# An instance
matlab = Matlab()


def reduce_load_Matlab():
    return matlab


def matlab_console():
    """
    This requires that the optional matlab program be installed and in
    your PATH, but no optional Sage packages need be installed.

    EXAMPLES::

        sage: matlab_console()                # optional - matlab; not tested
                                       < M A T L A B >
                           Copyright 1984-2006 The MathWorks, Inc.
        ...
        >> 2+3

    ans =

    5

    quit

    Typing quit exits the matlab console and returns you to Sage.
    matlab, like Sage, remembers its history from one session to
    another.
    """
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%matlab magics instead.')
    os.system('matlab -nodisplay')


def matlab_version():
    """
    Return the version of Matlab installed.

    EXAMPLES::

        sage: matlab_version()    # random; optional - matlab
        '7.2.0.283 (R2006a)'
    """
    return str(matlab('version')).strip()
