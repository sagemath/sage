r"""
Interface to several Rubik's cube solvers.

The first is by Michael Reid, and tries to find an optimal solution given
the cube's state, and may take a long time.
See http://www.math.ucf.edu/~reid/Rubik/optimal_solver.html

The second is by Eric Dietz, and uses a standard (?) algorithm to
solve the cube one level at a time. It is extremely fast, but often
returns a far from optimal solution.
See https://web.archive.org/web/20121212175710/http://www.wrongway.org/?rubiksource

The third is by Dik Winter and implements Kociemba's algorithm which
finds reasonable solutions relatively quickly, and if it is kept running
will eventually find the optimal solution.

AUTHOR:

   -- Optimal was written by Michael Reid <reid@math.ucf.edu> (2004)
   -- Cubex was written by Eric Dietz <root@wrongway.org> (2003)
   -- Kociemba was written by Dik T. Winter <dik.winter@cwi.nl> (1993)
   -- Initial interface by Robert Bradshaw (2007-08)
"""

########################################################################
#   Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
########################################################################

import pexpect
import time
from . import cleaner

from sage.cpython.string import bytes_to_str
from sage.groups.perm_gps.cubegroup import index2singmaster



# Can't seem to find consistency in letter ordering
# between us and them... These are copied from the source.
optimal_solver_tokens = ["UF", "UR", "UB", "UL", \
                        "DF", "DR", "DB", "DL", \
                        "FR", "FL", "BR", "BL", \
                        "FU", "RU", "BU", "LU", \
                        "FD", "RD", "BD", "LD", \
                        "RF", "LF", "RB", "LB", \
                        "UFR", "URB", "UBL", "ULF", \
                        "DRF", "DFL", "DLB", "DBR", \
                        "FRU", "RBU", "BLU", "LFU", \
                        "RFD", "FLD", "LBD", "BRD", \
                        "RUF", "BUR", "LUB", "FUL", \
                        "FDR", "LDF", "BDL", "RDB"]

# The input format.
optimal_solver_format = "UF UR UB UL DF DR DB DL FR FL BR BL UFR URB UBL ULF DRF DFL DLB DBR"


class SingNot:
    """
    This class is to resolve difference between various Singmaster notation.

    Case is ignored, and the second and third letters may be swapped.

    EXAMPLES::

        sage: from sage.interfaces.rubik import SingNot
        sage: SingNot("acb") == SingNot("ACB")
        True
        sage: SingNot("acb") == SingNot("bca")
        False
    """
    def __init__(self, s):
        self.rep = s
        self.canonical = (s[0] + "".join(sorted(s[1:]))).lower()

    def __eq__(self, other):
        return isinstance(other, SingNot) and other.canonical == self.canonical

    def __repr__(self):
        return self.rep

    def __hash__(self):
        return hash(self.canonical)


# This is our list
singmaster_list = [''] + [SingNot(index2singmaster(i + 1)) for i in range(48)]


class OptimalSolver:
    """
    Interface to Michael Reid's optimal Rubik's Cube solver.
    """
    __cmd = "optimal"

    def __init__(self, verbose=False, wait=True):
        self.verbose = verbose
        self.start()
        if wait:
            print("Initializing tables...")
            self.ready()
            print("Done.")

    def start(self):
        child = pexpect.spawn(self.__cmd)
        cleaner.cleaner(child.pid, self.__cmd)
        child.timeout = None
        self.child = child
        self._ready = False

    def stop(self):
        if self.child:
            self.child.sendline(chr(3))  # send ctrl-c
            self.child.sendline(chr(4))  # send ctrl-d
            self.child.close(True)
            self.child = None

    def ready(self):
        if not self._ready:
            self.child.expect('enter cube')
            self._ready = True

    def __call__(self, facets):
        return self.solve(facets)

    def solve(self, facets):
        """
        The initial startup and precomputation are substantial...

        .. TODO:: Let it keep searching once it found a solution?

        EXAMPLES::

            sage: from sage.interfaces.rubik import *    # optional - rubiks
            sage: solver = DikSolver()                   # optional - rubiks
            sage: solver = OptimalSolver()  # optional - rubiks # long time (28s on sage.math, 2012)
            Initializing tables...
            Done.
            sage: C = RubiksCube("R U")                  # optional - rubiks
            sage: solver.solve(C.facets())               # optional - rubiks
            'R  U'
            sage: C = RubiksCube("R U F L B D")          # optional - rubiks
            sage: solver.solve(C.facets())               # optional - rubiks
            'R  U  F  L  B  D'
            sage: C = RubiksCube("R2 D2")                # optional - rubiks
            sage: solver.solve(C.facets())               # optional - rubiks
            'R2 D2'
        """
        self.ready()
        self.child.sendline(self.format_cube(facets))
        self.child.expect(r"([LRUDBF'2 ]+)\s+\((\d+)q\*?, (\d+)f\*?\)")
        self.child.sendline(chr(3))  # send ctrl-c
        return bytes_to_str(self.child.match.groups()[0]).strip()

    def format_cube(self, facets):
        L = []
        optimal_solver_list = [SingNot(x) for x in optimal_solver_tokens]
        for f in optimal_solver_format.split(" "):
            ix = facets[singmaster_list.index(SingNot(f)) - 1]
            facet = singmaster_list[ix]
            L.append(optimal_solver_list[optimal_solver_list.index(facet)])
        return " ".join(str(f) for f in L)


move_map = {
    "LD":"L'",
    "LU":"L",
    "RD":"R",
    "RU":"R'",
    "FA":"F",
    "FC":"F'",
    "BA":"B'",
    "BC":"B",
    "UR":"U",
    "UL":"U'",
    "DR":"D'",
    "DL":"D"
}


class CubexSolver:

    __cmd = "cubex"

    def __call__(self, facets):
        return self.solve(facets)

    def solve(self, facets):
        """
        EXAMPLES::

            sage: from sage.interfaces.rubik import *      # optional - rubiks
            sage: C = RubiksCube("R U")                    # optional - rubiks
            sage: CubexSolver().solve(C.facets())          # optional - rubiks
            'R U'
            sage: C = RubiksCube("R U F L B D")            # optional - rubiks
            sage: sol = CubexSolver().solve(C.facets()); sol  # optional - rubiks
            "U' L' L' U L U' L U D L L D' L' D L' D' L D L' U' L D' L' U L' B' U' L' U B L D L D' U' L' U L B L B' L' U L U' L' F' L' F L' F L F' L' D' L' D D L D' B L B' L B' L B F' L F F B' L F' B D' D' L D B' B' L' D' B U' U' L' B' D' F' F' L D F'"
            sage: RubiksCube(sol) == C                     # optional - rubiks
            True
            sage: C = RubiksCube("R2 F'")                  # optional - rubiks
            sage: CubexSolver().solve(C.facets())          # optional - rubiks
            "R' R' F'"
            sage: C = RubiksCube().scramble()              # optional - rubiks
            sage: sol = CubexSolver().solve(C.facets())    # optional - rubiks
            sage: C == RubiksCube(sol)                     # optional - rubiks
            True
        """
        s = self.format_cube(facets)
        child = pexpect.spawn(self.__cmd+" "+s)
        ix = child.expect(['210.*?:', r'^5\d+(.*)'])
        if ix == 0:
            child.expect(['211', pexpect.EOF])
            moves = bytes_to_str(child.before).strip().replace(',', '').split(' ')
            return " ".join(move_map[m] for m in reversed(moves))
        else:
            s = child.after
            while child.expect([r'^5\d+', pexpect.EOF]) == 0:
                s += child.after
            raise ValueError(bytes_to_str(s))

    def format_cube(self, facets):
        colors = sum([[i]*8 for i in range(1, 7)], [])
        facet_colors = [0] * 54
        for i in range(48):
            f = facets[i]-1
            f += (f+4) // 8 # to compensate for the centers
            facet_colors[f] = colors[i]
        for i in range(6):
            facet_colors[i*9+4] = i+1
        return "".join(str(c) for c in facet_colors)


class DikSolver:

    __cmd = "dikcube"

    def __call__(self, facets):
        return self.solve(facets)

    def solve(self, facets, timeout=10, extra_time=2):
        """
        EXAMPLES::

            sage: from sage.interfaces.rubik import *   # optional - rubiks
            sage: C = RubiksCube().move("R U")          # optional - rubiks
            sage: DikSolver().solve(C.facets())         # optional - rubiks
            'R U'
            sage: C = RubiksCube().move("R U F L B D")  # optional - rubiks
            sage: DikSolver().solve(C.facets())         # optional - rubiks
            'R U F L B D'
            sage: C = RubiksCube().move("R2 F'")        # optional - rubiks
            sage: DikSolver().solve(C.facets())         # optional - rubiks
            "R2 F'"
        """
        cube_str = self.format_cube(facets)
        child = pexpect.spawn(self.__cmd + " -p")
        child.expect('Initialization done!')
        child.sendline(cube_str)

        # We use send(chr(4)) instead of sendeof in this case, since
        # child.sendoef() when run in the background with the Dik solver
        # sends a SIGTTOU which suspends the process -- this is very bad.
        # This is only a temporary workaround, and does not fix the problem
        # on OS X.  The Dik C program itself will need to be fixed.
        # See trac #1683. (TODO)      -- willem jp, wstein, mabshoff
        child.send(chr(4))
        #child.sendeof()

        ix = child.expect(['Solution[^\n]*:', pexpect.EOF, pexpect.TIMEOUT], timeout=timeout)
        if ix == 0:
            child.expect(['[^\n]+'])
            sol = child.after.strip()
            start_time = time.time()
            while extra_time > time.time() - start_time:
                ix = child.expect(['Solution[^\n]*:', pexpect.EOF, pexpect.TIMEOUT], timeout=extra_time - int(time.time() - start_time))
                if ix == 0:
                    child.expect(['[^\n]+'])
                    sol = child.after.strip()
                else:
                    extra_time = 0
            # format the string into our notation
            child.close(True)
            sol = bytes_to_str(sol)
            return ' '.join(self.rot_map[m[0]] + str(4 - int(m[1]))
                            for m in reversed(sol.split(' '))).replace('1', '').replace('3', "'")
        elif ix == 1:
            # invalid format
            child.close(True)
            raise ValueError(bytes_to_str(child.before))
        else:
            child.close(True)
            raise RuntimeError("timeout")

    def format_cube(self, facets):
        colors = sum([[i]*8 for i in range(1,7)], [])
        facet_colors = [0] * 54
        for i in range(48):
            f = self.facet_map.index(facets[i])
            facet_colors[f] = colors[i]
        # now do the centers
        facet_colors[4] = 1
        facet_colors[49] = 6
        for i in range(2,6):
            facet_colors[16+i*3] = i
        return "".join(str(c) for c in facet_colors)

    facet_map = [      1,  2,  3,                                \
                       4,  0,  5,                                \
                       6,  7,  8,                                \
           9, 10, 11, 17, 18, 19, 25, 26, 27, 33, 34, 35,        \
          12,  0, 13, 20,  0, 21, 28,  0, 29, 36,  0, 37,        \
          14, 15, 16, 22, 23, 24, 30, 31, 32, 38, 39, 40,        \
                      41, 42, 43,                                \
                      44,  0, 45,                                \
                      46, 47, 48,                                \
            ]


    # to compensate for different face naming
    rot_map = dict(zip("BLURDF", "ULFRBD"))


#    facet_map = [
#                      1,  2,  3,
#                      4,      6,
#                      7,  8,  9,
#          10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
#          23,     25, 26,     28,     29, 30, 31,     33,
#          34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45,
#                      46, 47, 48,
#                      49,     51,
#                      52, 53, 54,
#    ]
