# sage.doctest: needs sage.modules
r"""
Knot mosaics

A knot mosaic is a square array of elementary tiles whose endpoints
match along shared edges.  Mosaics were introduced by Lomonaco and
Kauffman [LK2008]_ as a combinatorial model for quantum knot systems.
They have since been used to study mosaic number [LHLO2014]_,
tabulation [LLPP2018]_, virtual knots [GH2020]_, and surface-link
variants [CK2024]_.  Mosaics give a combinatorial way to encode knots
and links, and can be converted to the planar diagram codes used by
:class:`sage.knots.link.Link`.

EXAMPLES::

    sage: from sage.knots.mosaic import Mosaic
    sage: H = Mosaic([[0, 2, 1, 0],
    ....:             [2, 9, 10, 1],
    ....:             [3, 10, 10, 4],
    ....:             [0, 3, 4, 0]])
    sage: H
    Mosaic of dimension 4
    sage: H.is_suitably_connected()
    True
    sage: H.find_crossings()
    [(1, 1), (1, 2), (2, 1), (2, 2)]
    sage: H.number_of_crossings()
    4
    sage: H.number_of_components()
    2
    sage: H.pd_code()
    [[3, 5, 4, 8], [5, 1, 6, 4], [7, 3, 8, 2], [6, 1, 7, 2]]

The matrix representation can be recovered as a Sage matrix::

    sage: H.matrix()
    [ 0  2  1  0]
    [ 2  9 10  1]
    [ 3 10 10  4]
    [ 0  3  4  0]

REFERENCES:

.. [LK2008] Samuel J. Lomonaco and Louis H. Kauffman,
   *Quantum knots and mosaics*, Quantum Information Processing 7
   (2008), 85-115. :doi:`10.1007/s11128-008-0076-7`

.. [LHLO2014] Hwa Jeong Lee, Kyungpyo Hong, Ho Lee, and Seungsang Oh,
   *Mosaic number of knots*, Journal of Knot Theory and its
   Ramifications 23 (2014), no. 13, 1450069.
   :doi:`10.1142/S0218216514500692`

.. [LLPP2018] Hwa Jeong Lee, Lewis D. Ludwig, Joseph Paat, and
   Amanda Peiffer, *Knot mosaic tabulation*, Involve 11 (2018),
   no. 1, 13-26. :doi:`10.2140/involve.2018.11.13`

.. [GH2020] Sandy Ganzell and Allison Henrich, *Virtual mosaic knot
   theory*, Journal of Knot Theory and its Ramifications 29 (2020),
   no. 14, 2050091. :doi:`10.1142/S0218216520500911`

.. [CK2024] Seonmi Choi and Jieon Kim, *Mosaics for immersed
   surface-links*, Topology and its Applications 353 (2024), 108961.
   :doi:`10.1016/j.topol.2024.108961`

AUTHORS:

- Andrew Tawfeek (2026-05-09): initial version
"""

# ****************************************************************************
#       Copyright (C) 2026 Andrew Tawfeek <atawfeek@uw.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from math import inf
from random import choice, randrange

oo = inf


TILE_CONNECTIONS = {
    0: (),
    1: ("left", "down"),
    2: ("right", "down"),
    3: ("up", "right"),
    4: ("left", "up"),
    5: ("left", "right"),
    6: ("up", "down"),
    7: (("down", "left"), ("up", "right")),
    8: (("down", "right"), ("left", "up")),
    9: (("down", "up"), ("left", "right")),
    10: (("left", "right"), ("down", "up")),
}

FOUR_POINT_TILES = {7, 8, 9, 10}
CROSSING_TILES = {9, 10}

TILE_ZOOM_MAPS = {
    0: ((0, 0, 0), (0, 0, 0), (0, 0, 0)),
    1: ((0, 0, 0), (5, 1, 0), (0, 6, 0)),
    2: ((0, 0, 0), (0, 2, 5), (0, 6, 0)),
    3: ((0, 6, 0), (0, 3, 5), (0, 0, 0)),
    4: ((0, 6, 0), (5, 4, 0), (0, 0, 0)),
    5: ((0, 0, 0), (5, 5, 5), (0, 0, 0)),
    6: ((0, 6, 0), (0, 6, 0), (0, 6, 0)),
    7: ((0, 3, 1), (1, 0, 3), (3, 1, 0)),
    8: ((2, 4, 0), (4, 0, 2), (0, 2, 4)),
    9: ((0, 6, 0), (5, 9, 5), (0, 6, 0)),
    10: ((0, 6, 0), (5, 10, 5), (0, 6, 0)),
}

TILE_9_UP_DOWN_ZOOM = ((2, 8, 1), (7, 10, 7), (3, 8, 4))

TILES_GOING_UP = {3, 4, 6, 7, 8, 9, 10}
TILES_GOING_DOWN = {1, 2, 6, 7, 8, 9, 10}
TILES_GOING_LEFT = {1, 4, 5, 7, 8, 9, 10}
TILES_GOING_RIGHT = {2, 3, 5, 7, 8, 9, 10}

OPPOSITE_DIRECTIONS = {
    "up": "down",
    "down": "up",
    "left": "right",
    "right": "left",
}


def _flatten(lst):
    """
    Flatten nested tuples of directions.
    """
    result = []
    for item in lst:
        if isinstance(item, (list, tuple)):
            result.extend(_flatten(item))
        else:
            result.append(item)
    return result


def _as_rows(mosaic_matrix):
    """
    Return ``mosaic_matrix`` as a tuple of tuples of integers.
    """
    if isinstance(mosaic_matrix, Mosaic):
        return mosaic_matrix._matrix
    if hasattr(mosaic_matrix, "rows"):
        rows = mosaic_matrix.rows()
    else:
        rows = mosaic_matrix

    matrix_rows = tuple(tuple(int(entry) for entry in row) for row in rows)
    if not matrix_rows:
        raise ValueError("a mosaic must have positive dimension")

    row_lengths = {len(row) for row in matrix_rows}
    if len(row_lengths) != 1:
        raise ValueError("all rows of a mosaic must have the same length")
    if len(matrix_rows) != row_lengths.pop():
        raise ValueError("a mosaic must be square")

    for row in matrix_rows:
        for entry in row:
            if entry not in TILE_CONNECTIONS:
                raise ValueError("mosaic tiles must be integers between 0 and 10")
    return matrix_rows


def opposite(direction):
    """
    Return the opposite direction.

    EXAMPLES::

        sage: from sage.knots.mosaic import opposite
        sage: opposite('left')
        'right'
    """
    try:
        return OPPOSITE_DIRECTIONS[direction]
    except KeyError:
        raise ValueError("direction must be 'up', 'down', 'left', or 'right'")


class MosaicTile:
    r"""
    A tile in a knot mosaic.

    INPUT:

    - ``tile`` -- integer between 0 and 10

    EXAMPLES::

        sage: from sage.knots.mosaic import MosaicTile
        sage: T = MosaicTile(9)
        sage: T.is_crossing()
        True
        sage: T.exit_path('left')
        'right'
    """
    def __init__(self, tile):
        tile = int(tile)
        if tile not in TILE_CONNECTIONS:
            raise ValueError("mosaic tiles must be integers between 0 and 10")

        self._tile = tile
        self._connections = TILE_CONNECTIONS[tile]
        self._orientation = []

    def __repr__(self):
        """
        Return a string representation of this tile.
        """
        return "Mosaic tile {}".format(self._tile)

    def number(self):
        """
        Return the tile number.

        EXAMPLES::

            sage: from sage.knots.mosaic import MosaicTile
            sage: MosaicTile(4).number()
            4
        """
        return self._tile

    def number_of_connection_points(self):
        """
        Return the number of boundary connection points of this tile.
        """
        if self._tile == 0:
            return 0
        if self._tile in range(1, 7):
            return 2
        return 4

    def number_of_strands(self):
        """
        Return the number of strands in this tile.
        """
        return self.number_of_connection_points() // 2

    def is_crossing(self):
        """
        Return whether this tile is a crossing.
        """
        return self._tile in CROSSING_TILES

    def connection_directions(self):
        """
        Return the connection directions of this tile.

        EXAMPLES::

            sage: from sage.knots.mosaic import MosaicTile
            sage: MosaicTile(7).connection_directions()
            (('down', 'left'), ('up', 'right'))
        """
        return self._connections

    def directions(self):
        """
        Return the flattened connection directions of this tile.
        """
        return tuple(_flatten(self._connections))

    def is_going(self, direction):
        """
        Return whether this tile has a connection in ``direction``.
        """
        return direction in self.directions()

    def exit_path(self, direction):
        """
        Return the exit direction for a strand entering from ``direction``.
        """
        if direction not in self.directions():
            raise ValueError("direction is not a connection direction")

        if self.number_of_strands() == 1:
            if direction == self._connections[0]:
                return self._connections[1]
            return self._connections[0]

        for strand in self._connections:
            if direction in strand:
                if strand[0] == direction:
                    return strand[1]
                return strand[0]

    def zoom(self, only_up_down=False):
        """
        Return the `3 \times 3` zoom replacement for this tile.
        """
        if self._tile == 9 and only_up_down:
            return TILE_9_UP_DOWN_ZOOM
        return TILE_ZOOM_MAPS[self._tile]

    def orient(self, direction):
        """
        Record an orientation through this tile.
        """
        if direction not in self.directions():
            raise ValueError("direction is not a connection direction")
        self._orientation.append(direction)

    def show(self, ax=None, resolution=5, color="blue"):
        """
        Draw this tile using matplotlib.

        INPUT:

        - ``ax`` -- a matplotlib axes object (default: ``None``)
        - ``resolution`` -- positive number (default: `5`)
        - ``color`` -- string (default: ``'blue'``)
        """
        import matplotlib.pyplot as plt
        from matplotlib import patches

        standalone = ax is None
        if standalone:
            _, ax = plt.subplots(figsize=(resolution, resolution))

        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(True)

        line_width = resolution
        tile = self._tile

        def arc(cx, cy, radius, theta1, theta2):
            ax.add_patch(patches.Arc((cx, cy), 2 * radius, 2 * radius,
                                     theta1=theta1, theta2=theta2,
                                     linewidth=line_width, fill=False,
                                     color=color))

        def line(x1, y1, x2, y2):
            ax.plot([x1, x2], [y1, y2], color=color,
                    linewidth=line_width, solid_capstyle="butt")

        if tile == 1:
            arc(0, 0, 0.5, 0, 90)
        elif tile == 2:
            arc(1, 0, 0.5, 90, 180)
        elif tile == 3:
            arc(1, 1, 0.5, 180, 270)
        elif tile == 4:
            arc(0, 1, 0.5, 270, 360)
        elif tile == 5:
            line(0, 0.5, 1, 0.5)
        elif tile == 6:
            line(0.5, 0, 0.5, 1)
        elif tile == 7:
            arc(0, 0, 0.5, 0, 90)
            arc(1, 1, 0.5, 180, 270)
        elif tile == 8:
            arc(1, 0, 0.5, 90, 180)
            arc(0, 1, 0.5, 270, 360)
        elif tile == 9:
            line(0, 0.5, 1, 0.5)
            line(0.5, 0, 0.5, 0.35)
            line(0.5, 0.65, 0.5, 1)
        elif tile == 10:
            line(0.5, 0, 0.5, 1)
            line(0, 0.5, 0.35, 0.5)
            line(0.65, 0.5, 1, 0.5)

        return ax


class Mosaic:
    r"""
    A knot mosaic.

    INPUT:

    - ``mosaic_matrix`` -- a square matrix, or a square list of lists,
      whose entries are tile numbers from 0 to 10
    """
    def __init__(self, mosaic_matrix):
        self._matrix = _as_rows(mosaic_matrix)
        self._size = len(self._matrix)

    def __repr__(self):
        """
        Return a string representation of this mosaic.
        """
        return "Mosaic of dimension {}".format(self._size)

    def __eq__(self, other):
        """
        Compare two mosaics.
        """
        if not isinstance(other, Mosaic):
            try:
                other = Mosaic(other)
            except (TypeError, ValueError):
                return False
        return self._matrix == other._matrix

    def __hash__(self):
        """
        Return a hash of this mosaic.
        """
        return hash(self._matrix)

    def __iter__(self):
        """
        Iterate over rows of this mosaic.
        """
        return iter(self._matrix)

    def __getitem__(self, key):
        """
        Return an entry or row of this mosaic.
        """
        if isinstance(key, tuple):
            i, j = key
            return self._matrix[i][j]
        return self._matrix[key]

    def size(self):
        """
        Return the dimension of this mosaic.
        """
        return self._size

    def rows(self):
        """
        Return the rows of this mosaic as tuples.
        """
        return self._matrix

    def matrix(self):
        """
        Return the matrix representation of this mosaic.
        """
        from sage.matrix.constructor import matrix
        from sage.rings.integer_ring import ZZ

        return matrix(ZZ, self._matrix, sparse=False)

    def show(self, resolution=5, color="blue"):
        """
        Draw this mosaic using matplotlib.
        """
        import matplotlib.pyplot as plt

        _, axes = plt.subplots(self._size, self._size,
                               figsize=(resolution, resolution),
                               gridspec_kw={"wspace": 0.15,
                                            "hspace": 0.15})
        if self._size == 1:
            axes = [[axes]]

        for i in range(self._size):
            for j in range(self._size):
                ax = axes[i][j] if self._size > 1 else axes[0][0]
                MosaicTile(self._matrix[i][j]).show(ax=ax, color=color)

    def directions(self, i, j):
        """
        Return the connection directions of the tile in position ``(i, j)``.
        """
        return MosaicTile(self._matrix[i][j]).directions()

    def is_suitably_connected(self):
        """
        Return whether all adjacent tile edges match.
        """
        for i in range(self._size):
            for j in range(self._size):
                tile = MosaicTile(self._matrix[i][j])

                if tile.is_going("up"):
                    if i == 0 or not MosaicTile(self._matrix[i - 1][j]).is_going("down"):
                        return False
                if tile.is_going("down"):
                    if i == self._size - 1 or not MosaicTile(self._matrix[i + 1][j]).is_going("up"):
                        return False
                if tile.is_going("left"):
                    if j == 0 or not MosaicTile(self._matrix[i][j - 1]).is_going("right"):
                        return False
                if tile.is_going("right"):
                    if (j == self._size - 1 or
                            not MosaicTile(self._matrix[i][j + 1]).is_going("left")):
                        return False
        return True

    def zoom(self, only_up_down=False):
        """
        Return the `3`-fold zoom of this mosaic.

        If ``only_up_down`` is ``True``, every 9-tile is replaced by the
        zoom pattern whose central crossing is a 10-tile.
        """
        rows = []
        for row in self._matrix:
            block_rows = [[], [], []]
            for tile in row:
                block = MosaicTile(tile).zoom(only_up_down=only_up_down)
                for i in range(3):
                    block_rows[i].extend(block[i])
            rows.extend(block_rows)
        return Mosaic(rows)

    def find_crossings(self):
        """
        Return the coordinates of the crossing tiles.
        """
        crossings = []
        for i, row in enumerate(self._matrix):
            for j, tile in enumerate(row):
                if tile in CROSSING_TILES:
                    crossings.append((i, j))
        return crossings

    def number_of_crossings(self):
        """
        Return the number of crossings in this mosaic.
        """
        return len(self.find_crossings())

    def exit_path(self, i, j, direction):
        """
        Return the next tile and exit direction after entering a tile.
        """
        exit_dir = MosaicTile(self._matrix[i][j]).exit_path(direction)
        next_positions = {
            "up": ((i - 1, j), "up"),
            "down": ((i + 1, j), "down"),
            "left": ((i, j - 1), "left"),
            "right": ((i, j + 1), "right"),
        }
        return next_positions[exit_dir]

    def shift(self, i, j, dictionary=False):
        """
        Return adjacent connected tile coordinates.

        If ``dictionary`` is ``True``, return a dictionary from directions
        to coordinates.
        """
        if not self.is_suitably_connected():
            raise ValueError("the mosaic must be suitably connected")

        directions = self.directions(i, j)
        directions_dict = {}
        if "up" in directions:
            directions_dict["up"] = (i - 1, j)
        if "down" in directions:
            directions_dict["down"] = (i + 1, j)
        if "left" in directions:
            directions_dict["left"] = (i, j - 1)
        if "right" in directions:
            directions_dict["right"] = (i, j + 1)

        if dictionary:
            return directions_dict
        return list(directions_dict.values())

    def walk(self, crossing, direction, path_list=False, tangent=False):
        """
        Walk from a crossing in the given direction to the next crossing.
        """
        all_crossings = self.find_crossings()
        if crossing not in all_crossings:
            raise ValueError("the starting tile must be a crossing")

        pos_x, pos_y = crossing
        direction_deltas = {
            "up": (-1, 0),
            "down": (1, 0),
            "left": (0, -1),
            "right": (0, 1),
        }

        current_direction = direction
        dx, dy = direction_deltas[current_direction]
        pos_x += dx
        pos_y += dy

        path = [crossing, (pos_x, pos_y)]

        while self._matrix[pos_x][pos_y] not in CROSSING_TILES:
            entrance = opposite(current_direction)
            current_direction = MosaicTile(self._matrix[pos_x][pos_y]).exit_path(entrance)
            dx, dy = direction_deltas[current_direction]
            pos_x += dx
            pos_y += dy
            path.append((pos_x, pos_y))

        incidence = opposite(current_direction)

        if path_list:
            return path
        if tangent:
            return (pos_x, pos_y), opposite(incidence)
        return (pos_x, pos_y), incidence

    def strand_of(self, tile, direction=None, direction_tracking=False,
                  verbose=False):
        """
        Trace the complete strand through ``tile``.
        """
        tile_type = self._matrix[tile[0]][tile[1]]
        if tile_type == 0:
            return []

        if direction is None:
            directions = self.directions(tile[0], tile[1])
            direction = opposite(choice(directions))

        start_tile = tile
        start_direction = direction
        path = []

        tile, direction = self.exit_path(start_tile[0], start_tile[1],
                                         opposite(start_direction))
        path.append((tile, direction))

        while not (tile == start_tile and direction == start_direction):
            tile, direction = self.exit_path(tile[0], tile[1],
                                             opposite(direction))
            path.append((tile, direction))

        if verbose:
            direction_tracking = True
            for step in path:
                print("Went {} into tile {}.".format(step[1], step[0]))

        if direction_tracking:
            return path
        return [tile for tile, direction in path]

    def strand_matrix(self):
        """
        Return the matrix of the number of strands in each tile.
        """
        rows = []
        for row in self._matrix:
            rows.append([MosaicTile(tile).number_of_strands() for tile in row])
        return Mosaic(rows).matrix()

    def _strand_count_rows(self):
        """
        Return strand counts as mutable rows.
        """
        return [[MosaicTile(tile).number_of_strands() for tile in row]
                for row in self._matrix]

    def strand_orientation_at(self, tile, previous_tile):
        """
        Return the induced orientation through a tile.
        """
        if previous_tile[0] < tile[0]:
            return "down"
        if previous_tile[0] > tile[0]:
            return "up"
        if previous_tile[1] < tile[1]:
            return "right"
        return "left"

    def strands(self):
        """
        Return the strands of this mosaic.
        """
        strand_list = []
        nonvisited = self._strand_count_rows()
        nonempty_tiles = []
        for i, row in enumerate(self._matrix):
            for j, tile in enumerate(row):
                if tile != 0:
                    nonempty_tiles.append((i, j))

        for tile in nonempty_tiles:
            if nonvisited[tile[0]][tile[1]] > 0:
                strand = self.strand_of(tile)
                strand_list.append(strand)
                for strand_tile in strand:
                    nonvisited[strand_tile[0]][strand_tile[1]] -= 1

        unique = []
        seen = set()
        for strand in strand_list:
            key = tuple(sorted(strand))
            if key not in seen:
                seen.add(key)
                unique.append(strand)
        return unique

    def number_of_components(self):
        """
        Return the number of connected components.
        """
        if not self.is_suitably_connected():
            raise ValueError("the mosaic must be suitably connected")
        return len(self.strands())

    def is_unknot(self):
        """
        Return whether this one-component mosaic is detected as the unknot.
        """
        if not self.is_suitably_connected():
            raise ValueError("the mosaic must be suitably connected")

        number_of_components = self.number_of_components()
        if number_of_components != 1:
            raise ValueError("is_unknot only works for knots")

        if self.number_of_crossings() == 0:
            return True

        try:
            import spherogram
        except ImportError as exc:
            raise ImportError("is_unknot requires spherogram") from exc

        knot = spherogram.Link(self.pd_code())
        knot.simplify("global")
        return int(knot.knot_floer_homology()["total_rank"]) == 1

    def local_frames(self):
        """
        Return the upward and rightward adjacent tiles at each crossing.
        """
        frames = []
        for crossing in self.find_crossings():
            shift_dict = self.shift(crossing[0], crossing[1], True)
            frames.append((shift_dict["up"], shift_dict["right"]))
        return frames

    def flip(self):
        """
        Return this mosaic flipped upside-down.
        """
        flip_map = {1: 4, 4: 1, 2: 3, 3: 2, 7: 8, 8: 7}
        flipped = []
        for row in reversed(self._matrix):
            flipped.append([flip_map.get(tile, tile) for tile in row])
        return Mosaic(flipped)

    def potential_tiles(self, i, j):
        """
        Return the tile numbers compatible with neighboring connections.
        """
        necessary_connections = []
        top_boundary = False
        bottom_boundary = False
        left_boundary = False
        right_boundary = False

        if i == 0:
            top_boundary = True
        elif self.directions(i - 1, j) == ():
            top_boundary = True
        elif "down" in self.directions(i - 1, j):
            necessary_connections.append("up")
        else:
            top_boundary = True

        if i == self._size - 1:
            bottom_boundary = True
        elif "up" in self.directions(i + 1, j):
            necessary_connections.append("down")

        if j == 0:
            left_boundary = True
        elif self.directions(i, j - 1) == ():
            left_boundary = True
        elif "right" in self.directions(i, j - 1):
            necessary_connections.append("left")
        else:
            left_boundary = True

        if j == self._size - 1:
            right_boundary = True
        elif "left" in self.directions(i, j + 1):
            necessary_connections.append("right")

        tile_set = []
        for tile in range(11):
            tile_directions = set(MosaicTile(tile).directions())
            if set(necessary_connections).issubset(tile_directions):
                tile_set.append(tile)

        if top_boundary:
            tile_set = [tile for tile in tile_set
                        if tile not in TILES_GOING_UP]
        if bottom_boundary:
            tile_set = [tile for tile in tile_set
                        if tile not in TILES_GOING_DOWN]
        if left_boundary:
            tile_set = [tile for tile in tile_set
                        if tile not in TILES_GOING_LEFT]
        if right_boundary:
            tile_set = [tile for tile in tile_set
                        if tile not in TILES_GOING_RIGHT]

        boundary_tile = (top_boundary or bottom_boundary or
                         left_boundary or right_boundary)
        if necessary_connections == [] and not boundary_tile:
            tile_set = [0]

        return tile_set

    def combine_components(self, tile=None, _depth=0):
        """
        Return a mosaic obtained by locally combining components.
        """
        if _depth > 5000:
            raise ValueError("could not satisfy constraints after 5000 attempts")

        matrix_rows = [list(row) for row in self._matrix]
        strands = self.strands()

        if tile is None:
            longest_strand = max(strands, key=len)
            tile = choice(longest_strand)

        strand = self.strand_of(tile)
        strand_counts = self._strand_count_rows()

        for strand_tile in strand:
            strand_counts[strand_tile[0]][strand_tile[1]] -= 1

        for strand_tile in strand:
            if strand_counts[strand_tile[0]][strand_tile[1]] == 1:
                i, j = strand_tile
                tile_type = matrix_rows[i][j]
                if tile_type == 7:
                    matrix_rows[i][j] = 8
                elif tile_type == 8:
                    matrix_rows[i][j] = 7
                elif tile_type in CROSSING_TILES:
                    matrix_rows[i][j] = choice([7, 8])

                return Mosaic(matrix_rows).combine_components(strand_tile,
                                                              _depth + 1)

        strand_set = set(strand)
        for i in range(self._size):
            for j in range(self._size):
                if (i, j) not in strand_set:
                    matrix_rows[i][j] = 0

        return Mosaic(matrix_rows)

    def oriented_gauss_code(self):
        """
        Return an oriented Gauss code compatible with :class:`Link`.
        """
        def pick_starting_tile():
            strand_matrix = self._strand_count_rows()
            for i in range(self._size):
                for j in range(self._size):
                    if strand_matrix[i][j] == 1:
                        return (i, j)
            raise ValueError("could not find a one-strand starting tile")

        def crossing_handedness(tile_type, orientation_pair):
            sorted_pair = sorted(orientation_pair)

            if tile_type == 9:
                if sorted_pair in [["right", "up"], ["down", "left"]]:
                    return 1
                if sorted_pair in [["down", "right"], ["left", "up"]]:
                    return -1
            else:
                if sorted_pair in [["left", "up"], ["down", "right"]]:
                    return 1
                if sorted_pair in [["down", "left"], ["right", "up"]]:
                    return -1

        def over_under(tile_type, orientation, numeric=False):
            if tile_type == 9:
                positioning = "under" if orientation in ["up", "down"] else "over"
            else:
                positioning = "over" if orientation in ["up", "down"] else "under"

            if numeric:
                return 1 if positioning == "over" else -1
            return positioning

        path = self.strand_of(pick_starting_tile())
        path = list(enumerate(path))

        crossings = self.find_crossings()
        appearances = []

        for crossing in crossings:
            for index, tile in path:
                if tile == crossing:
                    appearances.append((self._matrix[crossing[0]][crossing[1]],
                                        crossing, index, path[index - 1][1]))

        appearances.sort()

        orientations = []
        for appearance in appearances:
            tile, coord, index, previous_coord = appearance
            entrance = self.strand_orientation_at(coord, previous_coord)
            orientations.append((index, tile, coord, entrance,
                                 over_under(tile, entrance, numeric=True)))

        crossings = list(dict.fromkeys([crossing for index, tile, crossing,
                                        orientation, positioning
                                        in orientations]))

        crossing_orientations = []
        for crossing in crossings:
            tile_type = self._matrix[crossing[0]][crossing[1]]
            orientation_pair = [entrance for index, tile, crossing0, entrance,
                                positioning in orientations
                                if crossing0 == crossing]
            crossing_orientations.append(crossing_handedness(tile_type,
                                                             orientation_pair))

        orientations.sort()

        code_filter = [
            (crossings.index(crossing) + 1) * positioning
            for index, tile, crossing, entrance, positioning in orientations
        ]

        return [[code_filter], crossing_orientations]

    def pd_code(self):
        """
        Return the planar diagram code of this mosaic.
        """
        if self.number_of_crossings() == 0:
            raise ValueError("a crossing-free mosaic has no PD code")

        crossings = sorted(self.find_crossings())
        all_visits = {crossing: [] for crossing in crossings}
        arc_counter = 0
        ccw_order = ["right", "up", "left", "down"]

        remaining = {}
        for i in range(self._size):
            for j in range(self._size):
                tile_val = self._matrix[i][j]
                if tile_val == 0:
                    continue
                tile = MosaicTile(tile_val)
                if tile.number_of_strands() == 1:
                    remaining[(i, j)] = [frozenset(tile.connection_directions())]
                else:
                    remaining[(i, j)] = [frozenset(s)
                                         for s in tile.connection_directions()]

        while remaining:
            start = min(remaining.keys())
            strand_candidates = sorted(tuple(sorted(s))
                                       for s in remaining[start])
            chosen = strand_candidates[0]
            start_direction = opposite(chosen[0])

            walk = self.strand_of(start, direction=start_direction,
                                  direction_tracking=True)

            crossing_visits = []
            walk_length = len(walk)
            for k in range(walk_length):
                tile_t = tuple(walk[k][0])
                entry_dir = walk[k][1]
                exit_dir = walk[(k + 1) % walk_length][1]
                strand_used = frozenset({opposite(entry_dir), exit_dir})

                if tile_t in remaining:
                    try:
                        remaining[tile_t].remove(strand_used)
                    except ValueError:
                        pass
                    if not remaining[tile_t]:
                        del remaining[tile_t]

                if tile_t in all_visits:
                    crossing_visits.append((tile_t, entry_dir))

            number_of_visits = len(crossing_visits)
            if number_of_visits == 0:
                continue

            strand_arcs = [arc_counter + i + 1
                           for i in range(number_of_visits)]
            arc_counter += number_of_visits

            for idx, (tile, direc) in enumerate(crossing_visits):
                arc_in = strand_arcs[(idx - 1) % number_of_visits]
                arc_out = strand_arcs[idx]
                entry_port = opposite(direc)

                tile_type = self._matrix[tile[0]][tile[1]]
                if tile_type == 9:
                    is_over = direc in ("left", "right")
                else:
                    is_over = direc in ("up", "down")

                all_visits[tile].append({
                    "entry_port": entry_port,
                    "is_over": is_over,
                    "arc_in": arc_in,
                    "arc_out": arc_out,
                })

        pd_code = []
        for crossing in crossings:
            visits = all_visits[crossing]
            if len(visits) != 2:
                raise ValueError("crossing {} has {} visits"
                                 .format(crossing, len(visits)))

            port_to_arc = {}
            for visit in visits:
                port_to_arc[visit["entry_port"]] = visit["arc_in"]
                port_to_arc[opposite(visit["entry_port"])] = visit["arc_out"]
            if len(port_to_arc) != 4:
                raise ValueError("crossing {} does not have four arcs"
                                 .format(crossing))

            under = next(visit for visit in visits if not visit["is_over"])
            a_port = under["entry_port"]
            a_idx = ccw_order.index(a_port)

            pd_tuple = [port_to_arc[ccw_order[(a_idx + i) % 4]]
                        for i in range(4)]
            pd_code.append(pd_tuple)

        return pd_code

    def link(self):
        """
        Return the link represented by this mosaic.
        """
        from sage.knots.link import Link

        return Link(self.pd_code())


def random_mosaic(dimension, suitably_connected=True, number_of_crossings=None,
                  number_of_components=None, unknot=None):
    """
    Return a random mosaic satisfying the requested constraints.

    INPUT:

    - ``dimension`` -- positive integer
    - ``suitably_connected`` -- boolean (default: ``True``)
    - ``number_of_crossings`` -- integer or ``None`` (default: ``None``)
    - ``number_of_components`` -- integer or ``None`` (default: ``None``)
    - ``unknot`` -- boolean or ``None`` (default: ``None``)
    """
    dimension = int(dimension)
    if dimension <= 0:
        raise ValueError("dimension must be positive")
    if unknot is not None and not suitably_connected:
        raise ValueError("unknot requires suitably_connected=True")

    for _ in range(5001):
        if suitably_connected:
            template = [[0 for j in range(dimension)]
                        for i in range(dimension)]
            for i in range(dimension):
                for j in range(dimension):
                    template[i][j] = choice(Mosaic(template).potential_tiles(i, j))
            mosaic = Mosaic(template)
        else:
            template = [[randrange(11) for j in range(dimension)]
                        for i in range(dimension)]
            mosaic = Mosaic(template)

        if (number_of_crossings is not None and
                mosaic.number_of_crossings() != number_of_crossings):
            continue

        if number_of_components is not None:
            if (not mosaic.is_suitably_connected() or
                    mosaic.number_of_components() != number_of_components):
                continue

        if unknot is not None:
            try:
                if mosaic.is_unknot() != unknot:
                    continue
            except ValueError:
                continue

        return mosaic

    raise ValueError("could not satisfy constraints after 5000 attempts")


def _is_infinity(value):
    """
    Return whether ``value`` represents positive infinity.
    """
    if value == oo:
        return True
    return str(value) in ("+Infinity", "Infinity", "+infinity", "inf")


def rational_tangle(value, flip=False):
    """
    Return a rational tangle mosaic for ``value``.

    INPUT:

    - ``value`` -- infinity, zero, or an integer
    - ``flip`` -- boolean (default: ``False``)
    """
    if _is_infinity(value):
        return Mosaic([[7]])
    if value == 0:
        return Mosaic([[8]])

    value = int(value)
    size = abs(value)
    rows = [[0 for j in range(size)] for i in range(size)]
    for i in range(size):
        rows[i][i] = 10 if value > 0 else 9

    if flip:
        for i in range(size - 1):
            rows[i][i + 1] = 1
        for i in range(1, size):
            rows[i][i - 1] = 3
    else:
        for i in range(size - 1):
            rows[i][i + 1] = 4
        for i in range(1, size):
            rows[i][i - 1] = 2
        rows.reverse()

    return Mosaic(rows)


def _tangle_connector(n, m, direction):
    """
    Return a connector block for tangle joins.
    """
    if direction == "bottom-right":
        return [[6] + [0 for j in range(m - 1)] for i in range(n - 1)] + [
            [4] + [0 for j in range(m - 1)]
        ]
    if direction == "top-left":
        return [[0 for j in range(m)] for i in range(n - 1)] + [
            [2] + [5 for j in range(m - 1)]
        ]
    raise ValueError("unknown connector direction")


def tangle_join(tangle_list):
    """
    Join two rational tangles.
    """
    if len(tangle_list) != 2:
        raise ValueError("tangle_join currently supports two tangles")

    tangle0 = rational_tangle(tangle_list[0])
    tangle1 = rational_tangle(tangle_list[1])
    tangle0_flipped = rational_tangle(tangle_list[0], flip=True)

    top_left = _tangle_connector(tangle1.size(), tangle0.size(), "top-left")
    bottom_right = _tangle_connector(tangle0.size(), tangle1.size(),
                                     "bottom-right")

    top = [top_left[i] + list(tangle1.rows()[i])
           for i in range(tangle1.size())]
    bottom = [list(tangle0_flipped.rows()[i]) + bottom_right[i]
              for i in range(tangle0.size())]

    return Mosaic(top + bottom)
