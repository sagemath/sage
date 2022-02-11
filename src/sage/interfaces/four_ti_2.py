r"""
Interface to 4ti2

http://www.4ti2.de

You must have the 4ti2 Sage package installed on your computer
for this interface to work.

Use ``sage -i 4ti2`` to install the package.

AUTHORS:

- Mike Hansen (2009): Initial version.

- Bjarke Hammersholt Roune (2009-06-26): Added Groebner, made code
  usable as part of the Sage library and added documentation and some
  doctests.

- Marshall Hampton (2011): Minor fixes to documentation.
"""

#*****************************************************************************
#       Copyright (C) 2009 Mike Hansen <mhansen@gmail.com>
#       Copyright (C) 2009 Bjarke Hammersholt Roune <www.broune.com>
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer_ring import ZZ
from sage.features.four_ti_2 import FourTi2Executable

import os


class FourTi2(object):
    r"""
    This object defines an interface to the program 4ti2. Each command
    4ti2 has is exposed as one method.
    """
    def __init__(self, directory=None):
        r"""
        Initialize this object.

        INPUT:

        - ``directory`` -- 4ti2 only deals with files, and this is the
          directory that Sage will write input files to and run 4ti2
          in. Use an appropriate temporary directory if the value is
          ``None``.

        EXAMPLES::

            sage: from sage.interfaces.four_ti_2 import FourTi2
            sage: f = FourTi2("/tmp/")
            sage: f.directory()
            '/tmp/'
        """
        self._directory = directory

    ##################
    # Input / output #
    ##################

    def directory(self):
        r"""
        Return the directory where the input files for 4ti2 are
        written by Sage and where 4ti2 is run.

        EXAMPLES::

            sage: from sage.interfaces.four_ti_2 import FourTi2
            sage: f = FourTi2("/tmp/")
            sage: f.directory()
            '/tmp/'
        """
        from sage.misc.temporary_file import tmp_dir
        if self._directory is None:
            # we have to put this here rather than in the __init__
            # method since apparently importing sage.misc.misc does not
            # work until Sage is done starting up.
            self._directory = tmp_dir()
        return self._directory

    def temp_project(self):
        r"""
        Return an input project file name that has not been used yet.

        EXAMPLES::

            sage: from sage.interfaces.four_ti_2 import four_ti_2
            sage: four_ti_2.temp_project()
            'project_...'
        """
        n = 0
        while True:
            project = "project_%s" % n
            touch_file = os.path.join(self.directory(),project) + '.touch'
            if not os.path.exists(touch_file):
                break
            n += 1
        f = open(touch_file, 'w')
        f.write(' ')
        f.close()
        return project

    def write_matrix(self, mat, filename):
        r"""
        Write the matrix ``mat`` to the file ``filename`` in 4ti2 format.

        INPUT:

        - ``mat`` -- A matrix of integers or something that can be
          converted to that.
        - ``filename`` -- A file name not including a path.

        EXAMPLES::

            sage: from sage.interfaces.four_ti_2 import four_ti_2
            sage: four_ti_2.write_matrix([[1,2],[3,4]], "test_file")
        """
        from sage.matrix.constructor import matrix
        from sage.structure.element import is_Matrix
        if not is_Matrix(mat):
            mat = matrix(ZZ, mat)
        if mat.base_ring() != ZZ:
            mat = mat.change_ring(ZZ)

        self.write_array(mat, mat.nrows(), mat.ncols(), filename)

    def write_single_row(self, row, filename):
        r"""
        Write the list ``row`` to the file ``filename`` in 4ti2 format
        as a matrix with one row.

        INPUT:

        - ``row`` -- A list of integers.
        - ``filename`` -- A file name not including a path.

        EXAMPLES::

            sage: from sage.interfaces.four_ti_2 import four_ti_2
            sage: four_ti_2.write_single_row([1,2,3,4], "test_file")
        """
        self.write_array([row], 1, len(row), filename)

    def write_array(self, array, nrows, ncols, filename):
        r"""
        Write the matrix ``array`` of integers (can be represented as
        a list of lists) to the file ``filename`` in directory
        ``directory()`` in 4ti2 format. The matrix must have ``nrows``
        rows and ``ncols`` columns.

        INPUT:

        - ``array`` -- A matrix of integers. Can be represented as a list
          of lists.
        - ``nrows`` -- The number of rows in ``array``.
        - ``ncols`` -- The number of columns in ``array``.
        - ``file`` -- A file name not including a path.

        EXAMPLES::

            sage: from sage.interfaces.four_ti_2 import four_ti_2
            sage: four_ti_2.write_array([[1,2,3],[3,4,5]], 2, 3, "test_file")
        """
        f = open(os.path.join(self.directory(), filename), 'w')
        f.write("%s %s\n"%(nrows, ncols))
        for row in array:
            f.write(" ".join(map(str, row)))
            f.write("\n")
        f.close()

    def read_matrix(self, filename):
        r"""
        Read a matrix in 4ti2 format from the file ``filename`` in
        directory ``directory()``.

        INPUT:

        - ``filename`` -- The name of the file to read from.

        OUTPUT:

        The data from the file as a matrix over `\ZZ`.

        EXAMPLES::

            sage: from sage.interfaces.four_ti_2 import four_ti_2
            sage: four_ti_2.write_matrix([[1,2,3],[3,4,6]], "test_file")
            sage: four_ti_2.read_matrix("test_file")
            [1 2 3]
            [3 4 6]
        """
        from sage.matrix.constructor import matrix
        try:
            f = open(os.path.join(self.directory(), filename))
            lines = f.readlines()
            f.close()
        except IOError:
            return matrix(ZZ, 0, 0)

        nrows, ncols = map(ZZ, lines.pop(0).strip().split())
        return matrix(ZZ, nrows, ncols,
                      [[ZZ(_) for _ in line.strip().split()] for line in lines
                       if line.strip() != ""])

    def _process_input(self, kwds):
        r"""
        kwds is a dict, and the values are written to files with
        extension given by the keys, except for the keys ``self``
        and ``project``.

        This interesting method is intended to be called as the first
        thing going on in a method implementing some action of 4ti2,
        where the value of ``locals()`` is passed as the dict, thus
        achieving to write out many project files to the right places
        just by giving the parameters of the method names that are the
        extension of the corresponding files.

        Nothing is written if the value is None. Otherwise the value
        is written as a matrix to the file given by the value of the
        key ``'project'`` with extension given by the key.

        INPUT:

        - kwds -- A dict controlling what data is written to what files.

        OUTPUT:

        The value of the key ``project``.

        EXAMPLES::

            sage: from sage.interfaces.four_ti_2 import four_ti_2
            sage: pr = four_ti_2._process_input(
            ....:     {'project': "test_file",
            ....:      'self': None,
            ....:      'tst': [[1,2,3],[3,4,5]]})
            sage: four_ti_2.read_matrix("test_file.tst")
            [1 2 3]
            [3 4 5]
        """
        # Get the project
        project = kwds.get('project', None)
        if project is None:
            project = self.temp_project()

        for ext, value in kwds.items():
            if value is None:
                continue
            if ext == "project" or ext == "self":
                continue

            if (isinstance(value, list) and
                not (value and isinstance(value[0], list))):
                self.write_single_row(value, project + "." + ext)
            else:
                self.write_matrix(value, project + "." + ext)

        return project

    ############
    # Commands #
    ############

    def call(self, command, project, verbose=True, *, options=()):
        r"""
        Run the 4ti2 program ``command`` on the project named
        ``project`` in the directory ``directory()``.

        INPUT:

        - command -- The 4ti2 program to run.
        - project -- The file name of the project to run on.
        - verbose -- Display the output of 4ti2 if ``True``.
        - options -- A list of strings to pass to the program.

        EXAMPLES::

            sage: from sage.interfaces.four_ti_2 import four_ti_2
            sage: four_ti_2.write_matrix([[6,10,15]], "test_file")
            sage: four_ti_2.call("groebner", "test_file", False) # optional - 4ti2
            sage: four_ti_2.read_matrix("test_file.gro") # optional - 4ti2
            [-5  0  2]
            [-5  3  0]
        """
        import subprocess
        feature = FourTi2Executable(command)
        feature.require()
        executable = feature.executable
        options = " ".join(options)
        cmd = f'{executable} {options} {project}'
        if verbose is False:
            cmd += " > /dev/null 2> /dev/null"
        subprocess.call(cmd, shell=True, cwd=self.directory())

    def zsolve(self, mat=None, rel=None, rhs=None, sign=None, lat=None, project=None):
        r"""
        Run the 4ti2 program ``zsolve`` on the parameters. See
        ``http://www.4ti2.de/`` for details.

        EXAMPLES::

            sage: from sage.interfaces.four_ti_2 import four_ti_2
            sage: A = [[1,1,1],[1,2,3]]
            sage: rel = ['<', '<']
            sage: rhs = [2, 3]
            sage: sign = [1,0,1]
            sage: four_ti_2.zsolve(A, rel, rhs, sign) # optional - 4ti2
            [
                     [ 1 -1  0]
                     [ 0 -1  0]
            [0 0 1]  [ 0 -3  2]
            [1 1 0]  [ 1 -2  1]
            [0 1 0], [ 0 -2  1], []
            ]
            sage: four_ti_2.zsolve(lat=[[1,2,3],[1,1,1]]) # optional - 4ti2
            [
                         [1 2 3]
            [0 0 0], [], [1 1 1]
            ]

        """
        project = self._process_input(locals())
        self.call('zsolve', project, options=['-q'])
        return [self.read_matrix(project+'.'+ext) for ext in
                ['zinhom', 'zhom', 'zfree']]

    def qsolve(self, mat=None, rel=None, sign=None, project=None):
        r"""
        Run the 4ti2 program ``qsolve`` on the parameters. See
        ``http://www.4ti2.de/`` for details.

        EXAMPLES::

            sage: from sage.interfaces.four_ti_2 import four_ti_2
            sage: A = [[1,1,1],[1,2,3]]
            sage: four_ti_2.qsolve(A) # optional - 4ti2
            [[], [ 1 -2  1]]
        """
        project = self._process_input(locals())
        self.call('qsolve', project, options=['-q', '-parbitrary'])
        return [self.read_matrix(project+'.'+ext) for ext in
                ['qhom', 'qfree']]

    def rays(self, mat=None, project=None):
        r"""
        Run the 4ti2 program ``rays`` on the parameters. See
        ``http://www.4ti2.de/`` for details.

        EXAMPLES::

            sage: from sage.interfaces.four_ti_2 import four_ti_2
            sage: four_ti_2.rays(four_ti_2._magic3x3()) # optional - 4ti2
            [0 2 1 2 1 0 1 0 2]
            [1 0 2 2 1 0 0 2 1]
            [1 2 0 0 1 2 2 0 1]
            [2 0 1 0 1 2 1 2 0]
        """
        project = self._process_input(locals())
        self.call('rays', project, options=['-q', '-parbitrary'])
        return self.read_matrix(project+'.ray')

    def hilbert(self, mat=None, lat=None, project=None):
        r"""
        Run the 4ti2 program ``hilbert`` on the parameters. See
        ``http://www.4ti2.de/`` for details.

        EXAMPLES::

            sage: from sage.interfaces.four_ti_2 import four_ti_2
            sage: four_ti_2.hilbert(four_ti_2._magic3x3()) # optional - 4ti2
            [2 0 1 0 1 2 1 2 0]
            [1 0 2 2 1 0 0 2 1]
            [0 2 1 2 1 0 1 0 2]
            [1 2 0 0 1 2 2 0 1]
            [1 1 1 1 1 1 1 1 1]
            sage: four_ti_2.hilbert(lat=[[1,2,3],[1,1,1]])   # optional - 4ti2
            [2 1 0]
            [0 1 2]
            [1 1 1]
        """
        project = self._process_input(locals())
        self.call('hilbert', project, options=['-q'])
        return self.read_matrix(project+'.hil')

    def graver(self, mat=None, lat=None, project=None):
        r"""
        Run the 4ti2 program ``graver`` on the parameters. See
        ``http://www.4ti2.de/`` for details.

        EXAMPLES::

            sage: from sage.interfaces.four_ti_2 import four_ti_2
            sage: four_ti_2.graver([1,2,3]) # optional - 4ti2
            [ 2 -1  0]
            [ 3  0 -1]
            [ 1  1 -1]
            [ 1 -2  1]
            [ 0  3 -2]
            sage: four_ti_2.graver(lat=[[1,2,3],[1,1,1]])  # optional - 4ti2
            [ 1  0 -1]
            [ 0  1  2]
            [ 1  1  1]
            [ 2  1  0]
        """
        project = self._process_input(locals())
        self.call('graver', project, options=['-q'])
        return self.read_matrix(project+'.gra')

    def ppi(self, n):
        r"""
        Run the 4ti2 program ``ppi`` on the parameters. See
        ``http://www.4ti2.de/`` for details.

        EXAMPLES::

            sage: from sage.interfaces.four_ti_2 import four_ti_2
            sage: four_ti_2.ppi(3) # optional - 4ti2
            [-2  1  0]
            [ 0 -3  2]
            [-1 -1  1]
            [-3  0  1]
            [ 1 -2  1]

        """
        self.call('ppi', f'{n} 2> /dev/null')
        return self.read_matrix('ppi%s.gra'%n)

    def circuits(self, mat=None, project=None):
        r"""
        Run the 4ti2 program ``circuits`` on the parameters. See
        ``http://www.4ti2.de/`` for details.

        EXAMPLES::

            sage: from sage.interfaces.four_ti_2 import four_ti_2
            sage: four_ti_2.circuits([1,2,3]) # optional - 4ti2
            [ 0  3 -2]
            [ 2 -1  0]
            [ 3  0 -1]
        """
        project = self._process_input(locals())
        self.call('circuits', project, options=['-q', '-parbitrary'])
        return self.read_matrix(project+'.cir')

    def minimize(self, mat=None, lat=None):
        r"""
        Run the 4ti2 program ``minimize`` on the parameters. See
        ``http://www.4ti2.de/`` for details.

        EXAMPLES::

            sage: from sage.interfaces.four_ti_2 import four_ti_2
            sage: four_ti_2.minimize() # optional - 4ti2
            Traceback (most recent call last):
            ...
            NotImplementedError: 4ti2 command 'minimize' not implemented in Sage.
        """
        raise NotImplementedError("4ti2 command 'minimize' not implemented "
                                   "in Sage.")

    def groebner(self, mat=None, lat=None, project=None):
        r"""
        Run the 4ti2 program ``groebner`` on the parameters. This
        computes a Toric Groebner basis of a matrix. See
        ``http://www.4ti2.de/`` for details.

        EXAMPLES::

            sage: from sage.interfaces.four_ti_2 import four_ti_2
            sage: A = [6,10,15]
            sage: four_ti_2.groebner(A) # optional - 4ti2
            [-5  0  2]
            [-5  3  0]
            sage: four_ti_2.groebner(lat=[[1,2,3],[1,1,1]]) # optional - 4ti2
            [-1  0  1]
            [ 2  1  0]
        """
        project = self._process_input(locals())
        self.call('groebner', project, options=['-q', '-parbitrary'])
        return self.read_matrix(project+'.gro')

    def _magic3x3(self):
        r"""
        Return a matrix used for testing this class.

        EXAMPLES::

            sage: from sage.interfaces.four_ti_2 import four_ti_2
            sage: four_ti_2._magic3x3() # optional - 4ti2
            [ 1  1  1 -1 -1 -1  0  0  0]
            [ 1  1  1  0  0  0 -1 -1 -1]
            [ 0  1  1 -1  0  0 -1  0  0]
            [ 1  0  1  0 -1  0  0 -1  0]
            [ 1  1  0  0  0 -1  0  0 -1]
            [ 0  1  1  0 -1  0  0  0 -1]
            [ 1  1  0  0 -1  0 -1  0  0]

        """
        from sage.matrix.constructor import matrix
        return matrix \
            (ZZ, 7, 9,
             [[1, 1, 1, -1, -1, -1,  0,  0,  0],
              [1, 1, 1,  0,  0,  0, -1, -1, -1],
              [0, 1, 1, -1,  0,  0, -1,  0,  0],
              [1, 0, 1,  0, -1,  0,  0, -1,  0],
              [1, 1, 0,  0,  0, -1,  0,  0, -1],
              [0, 1, 1,  0, -1,  0,  0,  0, -1],
              [1, 1, 0,  0, -1,  0, -1,  0,  0]])

# The instance that should be used outside this file.
four_ti_2 = FourTi2()
