r"""
This module contains Cython code for a backtracking algorithm to solve Sudoku puzzles.

Once Cython implements closures and the ``yield`` keyword is possible, this can be moved into the ``sage.games.sudoku`` module, as part of the ``Sudoku.backtrack`` method, and this module can be banned.
"""


def backtrack_all(n, puzzle):
    r"""
    A routine to compute all the solutions to a Sudoku puzzle.

    INPUT:

        - ``n`` -- the size of the puzzle, where the array is an `n^2\times n^2` grid

        - ``puzzle`` -- list of the entries of the puzzle (1-based), in row-major order

    OUTPUT:

        A list of solutions, where each solution is a (1-based) list similar to ``puzzle``.

    TESTS:

    This is just a cursory test here, since eventually this code will move.
    See the `backtrack` method of the `Sudoku` class in the
    `sage.games.sudoku` module for more enlightening examples. ::

        sage: from sage.games.sudoku_backtrack import backtrack_all
        sage: c = [0, 0, 0, 0, 1, 0, 9, 0, 0, 8, 0, 0, 4, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 3, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 2, 0, 4, 0, 0, 0, 0, 0, 0, 0, 5, 8, 0, 6, 0, 0, 0, 0, 1, 3, 0, 7, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0]
        sage: print(backtrack_all(3, c))
        [[6, 5, 4, 3, 1, 2, 9, 8, 7, 8, 3, 1, 4, 7, 9, 5, 2, 6, 2, 9, 7, 6, 8, 5, 4, 1, 3, 4, 7, 2, 5, 3, 8, 6, 9, 1, 3, 8, 5, 1, 9, 6, 2, 7, 4, 9, 1, 6, 7, 2, 4, 3, 5, 8, 5, 6, 8, 9, 4, 7, 1, 3, 2, 7, 4, 3, 2, 5, 1, 8, 6, 9, 1, 2, 9, 8, 6, 3, 7, 4, 5]]

    ALGORITHM:

    We traverse a search tree depth-first.  Each level of the tree corresponds to a location in the grid, listed in row-major order.  At each location we maintain a list of the symbols which may be used in that location as follows.

    A location has "peers", which are the locations in the same row, column or box (sub-grid).  As symbols are chosen (or fixed initially) at a location, they become ineligible for use at a peer.  We track this in the ``available`` array where at each location each symbol has a count of how many times it has been made ineligible.  As this counter transitions between 0 and 1, the number of eligible symbols at a location is tracked in the ``card`` array.  When the number of eligible symbols at any location becomes 1, then we know that *must* be the symbol employed in that location.  This then allows us to further update the eligible symbols at the peers of that location.  When the number of the eligible symbols at any location becomes 0, then we know that we can prune the search tree.

    So at each new level of the search tree, we propagate as many fixed symbols as we can, placing them into a two-ended queue (``fixed`` and ``fixed_symbol``) that we process until it is empty or we need to prune.  All this recording of ineligible symbols and numbers of eligible symbols has to be unwound as we backup the tree, though.

    The notion of propagating singleton cells forward comes from an essay by Peter Norvig [sudoku:norvig]_.
    """
    cdef:
        # Arrays sizes are n^4, and n^2, with 3n^2-2n-1 for second slot of peers, n = 4
        int i, j, count, level, apeer
        int nsquare, npeers, nboxes
        int grid_row, grid_col, grid_corner
        int peers[256][39]
        int box[256]
        int available[256][16]
        int card[256]
        int hint, symbol, abox
        int feasible
        int nfixed[256]
        int fixed[256][256]
        int fixed_symbol[256][256]

        int process, asymbol, alevel, process_level, process_symbol

    # sanity check on size (types)
    # n is "base" dimension
    # nsquare is size of a grid
    # nboxes is total number of entries
    nsquare = n*n
    nboxes = nsquare * nsquare
    npeers = 3*n*n-2*n-1  # 2(n^2-1)+n^2-2n+1

    # "Peers" of a box are boxes in the same column, row or grid
    # Like the conflict graph when expressed as a graph coloring problem
    for level in range(nboxes):
        # location as row and column in square
        # grids are numbered similarly, in row-major order
        row = level // nsquare
        col = level %  nsquare
        grid_corner = (row - (row % n))*nsquare + (col - (col % n))
        grid_row = row // n
        grid_col = col // n
        count = -1
        # Peers' levels in same grid, but not the box itself
        for i in range(n):
            for j in range(n):
                grid_level = grid_corner + i*nsquare + j
                if grid_level != level:
                    count += 1
                    peers[level][count] = grid_level
        # Peers' levels in the same row, but not in grid
        for i in range(nsquare):
            if (i // 3 != grid_col):
                count += 1
                peers[level][count] = row*nsquare + i
        # Peers' levels in same column, but not in grid
        for i in range(nsquare):
            if (i // 3 != grid_row):
                count += 1
                peers[level][count] = col + i*nsquare

    # Initialize data structures
    # Make every symbol available initially for a box
    # And make set cardinality the size of symbol set
    for level in range(nboxes):
        box[level] = -1
        card[level] = nsquare
        for j in range(nsquare):
            available[level][j] = 0

    # For nonzero entries of input puzzle
    # (1) Convert to zero-based indexing
    # (2) Make a set of size 1 available initially
    for level in range(nboxes):
        # location as row and column in square
        # grids are numbered similarly, in row-major order
        hint = puzzle[level] - 1
        if hint != -1:
            # Limit symbol set at the hint's location to a singleton
            for j in range(nsquare):
                available[level][j] = 1
            available[level][hint] = 0
            card[level] = 1
            #  Remove hint from all peers' available symbols
            #  Track cardinality as sets adjust
            for i in range(npeers):
                apeer = peers[level][i]
                available[apeer][hint] += 1
                if available[apeer][hint] == 1:
                    card[apeer] -= 1

    # Start backtracking
    solutions = []
    level = 0
    box[level] = -1
    while (level > -1):
        symbol = box[level]
        if (symbol != -1):
            # restore symbols to peers
            for i in range(nfixed[level]):
                alevel = fixed[level][i]
                asymbol = fixed_symbol[level][i]
                for j in range(npeers):
                    abox = peers[alevel][j]
                    available[abox][asymbol] -= 1
                    if available[abox][asymbol] == 0:
                        card[abox] += 1
        # move sideways in search tree to next available symbol
        symbol +=  1
        while (symbol < nsquare) and (available[level][symbol] != 0):
            symbol += 1
        if symbol == nsquare:
            # fell off the end sideways, backup
            level = level - 1
        else:
            box[level] = symbol
            # Remove elements of sets, adjust cardinalities
            # Descend in search tree if no empty sets created
            # Can't break early at an empty set
            #   or we will confuse restore that happens immediately
            feasible = True
            fixed[level][0] = level
            fixed_symbol[level][0] = symbol
            count = 0
            process = -1
            while (process < count) and feasible:
                process += 1
                process_level = fixed[level][process]
                process_symbol = fixed_symbol[level][process]
                for i in range(npeers):
                    abox = peers[process_level][i]
                    available[abox][process_symbol] += 1
                    if available[abox][process_symbol] == 1:
                        card[abox] -= 1
                        if card[abox] == 0:
                            feasible = False
                        if card[abox] == 1:
                            count += 1
                            fixed[level][count] = abox
                            asymbol = 0
                            while (available[abox][asymbol] != 0):
                                asymbol += 1
                            fixed_symbol[level][count] = asymbol
            nfixed[level] = process+1
            if feasible:
                if level == nboxes - 1:
                    # Have a solution to save, stay at this bottom-most level
                    # Once Cython implements closures, a yield can go here
                    solutions.append([box[i]+1 for i in range(nboxes)])
                else:
                    level = level + 1
                    box[level] = -1
    return solutions
