"""Define the DynamicPlanarMapShow class, which allows pretty visualization of planar maps."""

import math

# we use python's float type and math functions to avoid using the more accurate but much slower sage types
# eg. we want pi - 1 to be stored as its floating-point value instead of this formal expression
# this is why all numbers in this file are written with "r"s: type(1.0r+math.pi) = float, but type(1.0 + pi) = sage.symbolic.expression.Expression
# EDIT: no longer needed once the file is .py
import time
from typing import Any

import matplotlib.artist
import matplotlib.gridspec
import matplotlib.backend_bases
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot

from sage.all import Graph, Permutation
from sage.graphs.maps.labelled_map import LabelledMap, nx, plt


class Vector2D:
    """Simple internal Vector2D class."""

    def __init__(self, x=0.0, y=0.0):
        r"""
        Initialize a 2D Vector with the given coordinates.

        INPUT:
        - ``x`` -- float (default: ``0.0``)
        - ``y`` -- float (default: ``0.0``)

        EXAMPLES::

            sage: from sage.graphs.maps.dynamic_planar_map_show import Vector2D
            sage: v = Vector2D(1, -3.5)
            sage: print(v)
            (1.0; -3.5)
        """
        self.x = float(x)
        self.y = float(y)

    def normSq(self) -> float:
        r"""
        Return the squared norm of self.

        EXAMPLES::

            sage: from sage.graphs.maps.dynamic_planar_map_show import Vector2D
            sage: Vector2D(-3, 4).normSq()
            25.0
        """
        return self.x * self.x + self.y * self.y

    def norm(self) -> float:
        r"""
        Return the norm of self.

        EXAMPLES::

            sage: from sage.graphs.maps.dynamic_planar_map_show import Vector2D
            sage: Vector2D(-3, 4).norm()
            5.0
        """
        return math.sqrt(self.normSq())

    def normalized(self) -> "Vector2D":
        r"""
        Return a new Vector2D corresponding to the normalized vector of self (i.e. same direction, unit length).

        EXAMPLES::

            sage: from sage.graphs.maps.dynamic_planar_map_show import Vector2D
            sage: Vector2D(-3, 4).normalized()
            Vector2D(-0.6; 0.8)
        """
        if self.norm() == 0.0:
            # hopefully, this should cancel out if trying to use a normalized
            # null vector
            return Vector2D(1.0, 0.0)

        return self / self.norm()

    def force_float(self) -> None:
        r"""
        Convert the attributes x and y to float.

        EXAMPLES::

            sage: from sage.graphs.maps.dynamic_planar_map_show import Vector2D
            sage: v = Vector2D()
            sage: v.x += pi + sin(1)
            sage: print(v)
            (pi + sin(1); 0.0)
            sage: v.force_float()
            sage: print(v)
            (3.98306363839769; 0.0)

        NOTE:
            This function helps avoid costly computation when directly updating x and y, by using fast float operations instead of computing symbolic sage expressions.
        """
        self.x = float(self.x)
        self.y = float(self.y)

    def rotate90(self) -> "Vector2D":
        r"""
        Return a Vector2D corresponding to a rotation of 90 degrees (pi/2 rad) counterclockwise of self.

        EXAMPLES::

            sage: from sage.graphs.maps.dynamic_planar_map_show import Vector2D
            sage: Vector2D(1, 0).rotate90()
            Vector2D(-0.0; 1.0)
            sage: Vector2D(-3, 2).rotate90()
            Vector2D(-2.0; -3.0)
        """
        return Vector2D(-self.y, self.x)

    def dot(self, other: "Vector2D") -> float:
        r"""
        Return the dot product of self and other.

        INPUT:
        - ``other`` -- Vector2D

        EXAMPLES::

            sage: from sage.graphs.maps.dynamic_planar_map_show import Vector2D
            sage: Vector2D(3, 1).dot(Vector2D(-2, 2))
            -4.0
            sage: Vector2D(5, 0).dot(Vector2D(0, 3))
            0.0
        """
        return self.x * other.x + self.y * other.y

    def angle_towards(self, other: "Vector2D") -> float:
        r"""
        Return the angle from self to other (counterclockwise) in the range [0, 2pi[.

        INPUT:
        - ``other`` -- Vector2D

        EXAMPLES::

            sage: from sage.graphs.maps.dynamic_planar_map_show import Vector2D
            sage: Vector2D(1, 0).angle_towards(Vector2D(0, 1))
            1.5707963267948966
            sage: Vector2D(2, 2).angle_towards(Vector2D(-3, -3))
            3.141592653589793
        """
        a = math.atan2(self.x * other.y - self.y * other.x,
                       self.x * other.x + self.y * other.y)
        if a < 0.0:
            return a + 2.0 * math.pi
        return a

    def __add__(self, other: "Vector2D") -> "Vector2D":
        """Return self + other."""
        return Vector2D(self.x + other.x, self.y + other.y)

    def __sub__(self, other: "Vector2D") -> "Vector2D":
        """Return self - other."""
        return Vector2D(self.x - other.x, self.y - other.y)

    def __neg__(self) -> "Vector2D":
        """Return -self."""
        return Vector2D(-self.x, -self.y)

    def __mul__(self, sc: float) -> "Vector2D":
        """Return self * sc."""
        return Vector2D(self.x * sc, self.y * sc)

    def __rmul__(self, sc: float) -> "Vector2D":
        """Return self * sc."""
        return Vector2D(self.x * sc, self.y * sc)

    def __truediv__(self, sc: float) -> "Vector2D":
        """Return self / sc."""
        return Vector2D(self.x / sc, self.y / sc)

    def __str__(self) -> str:
        """Return a string representation of self."""
        return "(" + str(self.x) + "; " + str(self.y) + ")"

    def __repr__(self) -> str:
        """Return a string representation of self."""
        return "Vector2D" + str(self)

    def __eq__(self, other: Any) -> bool:
        """Return self == other."""
        if isinstance(other, Vector2D):
            return self.x == other.x and self.y == other.y
        return False

    def __lt__(self, other: "Vector2D") -> bool:        # lexicographical order
        """Return self <= other (in lexicographical order)."""
        return self.x < other.x or (self.x == other.x and self.y < other.y)

    def __hash__(self) -> int:
        """Return hash(self)."""
        return hash((self.x, self.y))


def check_segments_intersecting(p1: Vector2D, q1: Vector2D, p2: Vector2D, q2: Vector2D) -> bool:
    r"""
    Return True if the segment between the points [p1, q1] intersects the segment between points [p2, q2], False otherwise.

    INPUT:
    - ``p1`` -- Vector2D;
    - ``q1`` -- Vector2D;
    - ``p2`` -- Vector2D;
    - ``q2`` -- Vector2D;

    EXAMPLES::
        sage: from sage.graphs.maps.dynamic_planar_map_show import Vector2D, check_segments_intersecting
        sage: check_segments_intersecting(Vector2D(0, 0), Vector2D(5, 5), Vector2D(1, 3), Vector2D(3, -2))
        True
        sage: check_segments_intersecting(Vector2D(0, 0), Vector2D(5, 5), Vector2D(1, 3), Vector2D(3, 10))
        False
        sage: check_segments_intersecting(Vector2D(0, 0), Vector2D(5, 5), Vector2D(2,3), Vector2D(0,0))
        True
    """

    def onSegment(p: Vector2D, q: Vector2D, r: Vector2D) -> bool:
        "Check if q lies on [p, r] assuming p,q,r are colinear"
        return q.x <= max(p.x, r.x) and q.x >= min(
            p.x, r.x) and q.y <= max(p.y, r.y) and q.y >= min(p.y, r.y)

    def orient(p: Vector2D, q: Vector2D, r: Vector2D) -> int:
        "Returns the orientation of (p, q, r): 0 if colinear, 1 if clockwise, 2 if counterclockwise"

        val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y)

        if val == 0:
            return 0
        elif val > 0:
            return 1
        return 2

    o1, o2, o3, o4 = orient(p1, q1, p2), orient(
        p1, q1, q2), orient(p2, q2, p1), orient(p2, q2, q1)

    return (
        (o1 != o2 and o3 != o4) or
        (o1 == 0 and onSegment(p1, p2, q1)) or
        (o2 == 0 and onSegment(p1, q2, q1)) or
        (o3 == 0 and onSegment(p2, p1, q2)) or
        (o4 == 0 and onSegment(p2, q1, q2)))


def check_polygon_intersecting(segments: list[tuple[Vector2D, Vector2D]]) -> bool:
    r"""
    Return True if the polygon defined by the given segments (list of pairs Vector2D) intersects itself (except in its vertices).

    INPUT:
    - ``segments`` -- list[(Vector2D, Vector2D)];

    EXAMPLES::

        sage: from sage.graphs.maps.dynamic_planar_map_show import Vector2D, check_polygon_intersecting
        sage: p1, p2, p3, p4 = Vector2D(0,0), Vector2D(5,0), Vector2D(5,5), Vector2D(0,5)
        sage: check_polygon_intersecting([(p1,p2),(p2,p3),(p3,p1)])
        False
        sage: check_polygon_intersecting([(p1,p2),(p2,p3),(p3,p4),(p4,p1)])
        False
        sage: check_polygon_intersecting([(p1,p2),(p2,p4),(p4,p3),(p3,p1)])
        True


    NOTE:
        Complexity is O(n^2), where n is the number of segments; this could be improved to O(n log n) in a future implementation.
        Note that this function is only supposed to be used by DynamicPlanarMapShow.
    """

    # TODO: do this on O(n log n) instead of O(N^2)!

    for (p1, q1) in segments:
        for (p2, q2) in segments:
            # this case is already handled by the embedding check
            if p1 in (p2, q2) or q1 in (p2, q2):
                continue

            if check_segments_intersecting(p1, q1, p2, q2):
                return True

    return False


class DynamicPlanarMapShow:
    """Allows pretty visualization of planar maps."""

    # force expressions are heavily inspired by planarmap.js
    # (https://github.com/tgbudd/planarmap.js)

    # controls the vertex-vertex repulsion force between red nodes (real nodes)
    repulsionVVRRCoef = 3.0
    # controls the vertex-vertex repulsion force between a red node and a
    # white node (artificial node to split multiedges & loops)
    repulsionVVRWCoef = 2.0
    # controls the vertex-vertex repulsion force between white nodes
    repulsionVVWWCoef = 2.0

    # controls the torsion force between the consecutive edges around a vertex
    torsionCoef = 3.0
    torsionPower = 3.0                  # exponent of the angle in the torsion force
    # angles are divided by this quantity in the torsion force
    torsionScale = math.pi / 6

    springCoef = 3.0  # controls the strength of the spring force for each edge
    springLength = 1.0                 # controls the default length of an edge
    # if true, spring force is proportional to log(r / springLength); if
    # false, proportional to r - springLength
    useLogSpring = True

    # controls the strength of the repulsive force between nodes & edges
    repulsionVertexEdgeCoef = 3.0
    repulsionVertexEdgePower = 2.0      # exponent of the node-edge distance

    # controls the strength of the force pulling all nodes towards (0, 0)
    weakSpringCoef = 0.01

    # maximum distance a node is allowed to travel during a single iteration
    maxDispl = springLength * 2.0

    # first model
    fast_delta_t = 1.0                  # delta_t of the first iterations
    slow_delta_t = 0.05                 # delta_t of the last iterations
    decrease_power = 1.5                # delta_t = slow_delta_t * (frame / fullForceFrames) ** power in the second phase

    # second model
    default_delta_t = 1.0              # delta_t of the first iteration
    # delta_t will equal to default_delta_t * exp(frame / (exp_decrease *
    # n_edges))
    exp_decrease = 0.2

    use_exponential_model = False

    # consider that the optimum is found when all nodes move less than this
    # value in a single tick (when rescaled to [0,1])
    convergence_limit = 0.0005

    max_iter_tick = 3                  # max number of iterations during a single tick

    # print information about time use of different functions
    time_profile = True

    useVVForce = True
    useTorsionForce = True
    useSpringForce = True
    useWeakSpringForce = True
    useEEForce = True

    # for each valid frame, current_penalty decrease by one
    current_penalty = 0
    # for each wrong frame, current_penalty increase by this value
    wrong_move_penalty = 1

    validFrames = 0
    # start increasing current_delta_t -> delta_t after this number of valid
    # frames
    validFramesThreshold = 20

    def __init__(self, map: LabelledMap, break_down_num: int = 5):
        """
        Initializes DynamicPlanarMapShow object.

        INPUT:
            - ``map`` -- LabelledMap: the map to show
            - ``break_down_num`` -- int: each cycle and multiedge will be split into this number of smaller edges to allow smoother visualization

        EXAMPLES::

            sage: lm = LabelledMap(Permutation([(1,3), (2,), (4,)]), Permutation([(1,2), (3,4)]))
            sage: dyn_show = DynamicPlanarMapShow(lm)
        """
        # initialize the map

        self.map = map
        self.break_down_num = break_down_num

        vertices = map.sigma.to_cycles()

        real_n_vertices = len(vertices)  # Remove multiedges and loops

        self.real_n_vertices = real_n_vertices

        alpha = map.alpha
        sigma = map.sigma
        m = map.m

        def minmax(i: int, j: int) -> tuple[int, int]:
            "Returns (min(i,j), max(i,j)). Used to ensure edges always go from lowest to highest vertex id."""
            return min(i, j), max(i, j)

        # Map half-edge i to its corresponding vertex
        corres = [0] * (2 * m + 1)
        for i in range(1, len(vertices) + 1):
            for k in vertices[i - 1]:
                corres[k] = i

        # Three dictionaries to store half-edge IDs
        self.edge_labels_head = {}  # (i, j): half-edge from i to j
        self.edge_labels_tail = {}  # (i, j): half-edge from j to i
        self.edge_labels_middle = {}  # Used for loops & multiedges

        def rem(i: int) -> None:
            "Remove every occurrence of the value i in edge_labels_head and edge_labels_tail."
            for d in (self.edge_labels_head, self.edge_labels_tail):
                for (key, val) in list(d.items()):
                    if val == i:
                        del d[key]

        def break_down(i: int, break_down_num: int) -> None:
            """Add a new vertex v, and break down the edge whose half-edges are i & alpha(i) into ``break_down_num``
                edges (i, 2*m+1), (2*m+2, 2*m+3), .., (2*m+2*(break_down_num-1), alpha(i))."""
            nonlocal alpha, sigma, m

            # print ("breaking down", i)

            # self.edge_labels_middle[(corres[i] - 1, len(vertices))] = i# alpha(i)
            # self.edge_labels_middle[(
            #     corres[alpha(i)] - 1, len(vertices) + break_down_num - 2)] = alpha(i)#i

            self.edge_labels_middle[(corres[i] - 1, len(vertices) + break_down_num - 2)] = i  # alpha(i)
            self.edge_labels_middle[(
                corres[alpha(i)] - 1, len(vertices))] = alpha(i)  # i

            rem(i)
            rem(alpha(i))

            alpha_cycles = [(alpha(i), 2 * m + 1, i, 2 * m + 2 * (break_down_num - 1))] + \
                [(2 * k, 2 * k + 1)
                 for k in range(m + 1, m + break_down_num - 1)]
            # for some unknown reason, the typechecker assumes that Permutation needs two arguments
            alpha *= Permutation(alpha_cycles)  # type: ignore

            sigma_cycles = [(2 * k - 1, 2 * k)
                            for k in range(m + 1, m + break_down_num)]
            sigma *= Permutation(sigma_cycles)  # type: ignore

            for k in range(break_down_num - 1):
                corres.append(len(vertices) + 1)
                corres.append(len(vertices) + 1)

                vertices.append((2 * m + 1, 2 * m + 2))
                m += 1

        def break_loop(i: int):
            "Breaks the loop starting from i."
            break_down(i, max(self.break_down_num, 3))

        # For each loop a-a, add a new vertex v and replace the edge a-a
        # with two edges a-v, v-a.
        for i in range(1, 2 * m + 1):
            if corres[i] == corres[alpha(i)]:
                break_loop(i)

        # Handle each vertex and break down edges if needed.
        for v in range(1, len(vertices) + 1):
            seen_vertices = {}
            for i in vertices[v - 1]:
                if corres[alpha(i)] in seen_vertices:
                    if seen_vertices[corres[alpha(i)]] != -1:
                        duplicate_he = seen_vertices[corres[alpha(i)]]
                        seen_vertices[corres[alpha(i)]] = -1
                        break_down(duplicate_he, self.break_down_num)
                    break_down(i, self.break_down_num)
                else:
                    seen_vertices[corres[alpha(i)]] = i
                    if corres[i] <= real_n_vertices and corres[alpha(
                            i)] <= real_n_vertices:
                        if corres[i] < corres[alpha(
                                i)] and i not in self.edge_labels_middle.values():
                            self.edge_labels_head[minmax(
                                corres[i] - 1, corres[alpha(i)] - 1)] = i
                        elif corres[i] > corres[alpha(i)] and i not in self.edge_labels_middle.values():
                            self.edge_labels_tail[minmax(
                                corres[i] - 1, corres[alpha(i)] - 1)] = i

        # Build the graph embedding
        embedding = {
            i: list(corres[alpha(k)] for k in vertices[i - 1])[::-1]
            for i in range(1, len(vertices) + 1)
        }

        edges_num1 = [       # vertices must be numbered from 1 to nVertices included in networkx
            (corres[i], corres[alpha(i)]) for i in range(1, 2 * m + 1) if i < alpha(i)
        ]

        self.alpha, self.sigma = alpha, sigma

        self.nVertices = len(vertices)
        self.nEdges = m

        # initialize positions to a valid (but ugly) layout

        # but internally, it's more convenient to use zero-indexing
        self.edges = [(i - 1, j - 1) for (i, j) in edges_num1]
        self.faces = [[corres[j] - 1 for j in face]
                      for face in alpha.right_action_product(sigma).to_cycles()]
        self.embedding = [[j - 1 for j in embedding[i]]
                          for i in range(1, self.nVertices + 1)]

        g = Graph(edges_num1, loops=False, multiedges=False)
        g.set_embedding(embedding)

        self.G = nx.DiGraph()
        self.G.add_edges_from(minmax(i - 1, j - 1) for (i, j, _) in g.edges())

        layout = g.layout_planar() if map.genus() == 0 else g.layout()

        self.pos = [Vector2D(*layout[i + 1]) for i in range(self.nVertices)]
        # self.speed = [Vector2D(0, 0) for i in range(self.nVertices)]

        self.max_extent = 1

        self.anim_running = True
        self.done = False

        self.is_planar = map.genus() == 0

    def start(self, show_halfedges: bool | str = "auto", plt_show=True, frame_by_frame=False) -> None:
        """
        Dynamically show the map in a new matplotlib figure.

        INPUT:
            - ``show_halfedges`` -- bool or "auto"; whether to show halfedges (if "auto", it will be set to True if nEdges <= 10)
            - ``plt_show`` -- bool; whether to call plt.show(). if False, user needs to call plt.ion() before start() to show the map in a non-blocking way.
            - ``frame_by_frame`` -- bool; if False, automatically advance until convergence is found (see the Note).

        EXAMPLES::

            sage: lm = LabelledMap(Permutation([(1,3), (2,), (4,)]), Permutation([(1,2), (3,4)]))
            sage: dyn_show = DynamicPlanarMapShow(lm)
            sage: dyn_show.start()
            ...

        NOTE:
            Press Enter to pause or resume the computation, Space to compute a single frame if the animation is currently paused, and Q/Esc to quit.
            If you encounter the warning, ``UserWarning: FigureCanvasAgg is non-interactive, and thus cannot be shown``, install PyQt6 using
            ``sage -pip install PyQt6`` (this warning can be displayed when calling ``dyn_show.start()``, hence the ``...`` in the example).
        """
        # initialize the matplotlib figure
        if show_halfedges == "auto":
            show_halfedges = self.nEdges <= 10

        if not isinstance(show_halfedges, bool):
            raise ValueError("Invalid value for show_halfedges")
        size = 7

        self.fig = plt.figure(figsize=(size, size + 1), layout="constrained")
        gs = matplotlib.gridspec.GridSpec(2, 2, figure=self.fig, width_ratios=[
                                          1, 4], height_ratios=[1, size])

        self.ax = self.fig.add_subplot(gs[1, :])
        txt_ax = self.fig.add_subplot(gs[0, :])
        # slider_ax = self.fig.add_subplot(gs[0, 1])

        pos_centered = self.center_pos()
        pos_centered_dict = {i: pos_centered[i] for i in range(self.nVertices)}

        # we need # type: ignore because the node_size argument is typed as int only, even though both
        # the documentation and the implementation make it clear that you can also use a list of int
        self.nodes_plt = nx.draw_networkx_nodes(self.G, pos_centered_dict, ax=self.ax, nodelist=list(range(self.nVertices)),
                                                # type: ignore
                                                node_color="red", node_size=[min(300, 1000 / self.nVertices**.5)] * self.real_n_vertices + [0] * (self.nVertices - self.real_n_vertices))
        self.edges_plt = nx.draw_networkx_edges(
            self.G, pos_centered_dict, ax=self.ax, arrows=False)

        # if show_halfedges:
        #     self.labels_head = nx.draw_networkx_edge_labels(
        #         self.G, pos_centered, ax=self.ax, rotate=False, edge_labels=self.edge_labels_head, label_pos=0.7)
        #     self.labels_tail = nx.draw_networkx_edge_labels(
        #         self.G, pos_centered, ax=self.ax, rotate=False, edge_labels=self.edge_labels_tail, label_pos=0.3)
        #     self.labels_middle = nx.draw_networkx_edge_labels(
        #         self.G, pos_centered, ax=self.ax, rotate=False, edge_labels=self.edge_labels_middle, label_pos=0.5)

        if show_halfedges:
            self.edge_texts = []
            for (d, prop) in ((self.edge_labels_head, 0.7),
                              (self.edge_labels_tail, 0.3), (self.edge_labels_middle, 0.5)):
                for (pair, txt) in d.items():
                    x = pos_centered[pair[0]][0] * prop + \
                        pos_centered[pair[1]][0] * (1 - prop)
                    y = pos_centered[pair[0]][1] * prop + \
                        pos_centered[pair[1]][1] * (1 - prop)
                    self.edge_texts.append((self.ax.text(x, y, txt, ha="center", va="center", bbox={
                                           "facecolor": "white", "edgecolor": "white"}), pair[0], pair[1], prop))

        self.show_halfedges = show_halfedges

        self.ax.axis("off")
        self.ax.set_xlim(left=0, right=1)
        self.ax.set_ylim(bottom=0, top=1)

        txt_ax.axis("off")
        txt_ax.set_xlim(left=0, right=1)
        txt_ax.set_ylim(bottom=0, top=1)

        # slider_ax.axis("off")

        self.text = txt_ax.text(0, 0.5, "Frame: 0")

        # self.slider = matplotlib.widgets.Slider(ax = slider_ax, valmin = -3, valmax = 1, valinit = 0, label = "delta_t")
        # self.slider.on_changed(self.update_slider_val)
        # self.slider.valtext.set_text("1.000")
        # self.force_delta_t = False

        # self.fig.tight_layout()

        # self.fig.canvas.mpl_connect('button_press_event', self.onClick)
        self.fig.canvas.mpl_connect('key_press_event', self.onKey)

        self.doFrame = False
        self.frame = 0

        if frame_by_frame:
            self.anim_running = False
            # self.anim.event_source.stop()
        else:
            self.anim_running = True

        self.reset_forces()

        if self.time_profile:
            self.forces_times = {}
            for force in self.forces_to_compute:
                self.forces_times[force.__name__] = 0.0

            self.update_forces_time = 0.0
            self.check_pos_correct_time = 0.0

        if self.nVertices > 30:
            self.useLogSpring = True
        else:
            self.useLogSpring = False

        self.validFrames = 0
        self.currentMaxDispl = self.maxDispl
        self.prev_pos_centered = None

        self.anim = FuncAnimation(
            self.fig, self.update_fig, cache_frame_data=False, blit=True, interval=1)

        self.ax.callbacks.connect(
            'xlim_changed', lambda event: setattr(self.anim, "_blit_cache", {}))
        self.ax.callbacks.connect(
            'ylim_changed', lambda event: setattr(self.anim, "_blit_cache", {}))

        if plt_show:
            plt.show()

    def reset_forces(self) -> None:
        r"""
        Clear and fill the attribute forces_to_compute with each force method, provided the corresponding ``use..Force`` attribute is set to True.
        """
        self.forces_to_compute = []

        if self.useVVForce:
            self.forces_to_compute.append(self.computeRepulsionVVForces)
        if self.useSpringForce:
            self.forces_to_compute.append(self.computeSpringForces)
        if self.useWeakSpringForce:
            self.forces_to_compute.append(self.computeWeakSpringForces)
        if self.useTorsionForce:
            self.forces_to_compute.append(self.computeTorsionForces)
        if self.useEEForce:
            self.forces_to_compute.append(self.computeRepulsionEEForces)

    def center_pos(self) -> list[tuple[float, float]]:
        r"""
        Center ``self.pos`` to new positions in [0,1]*[0,1] (with a .05 margin on each side), return the result as a list of pairs (x,y).

        EXAMPLES::

            sage: from sage.graphs.maps.dynamic_planar_map_show import Vector2D
            sage: lm = LabelledMap(Permutation([(1,),(2,)]),Permutation([(1,2)]))
            sage: d = DynamicPlanarMapShow(lm)
            sage: d.pos = [Vector2D(1,1), Vector2D(3,2)]
            sage: d.center_pos()
            [(0.05, 0.05), (0.9500000000000001, 0.9500000000000001)]
        """

        xs, ys = [v.x for v in self.pos], [v.y for v in self.pos]

        minx, maxx, miny, maxy = min(xs), max(xs), min(ys), max(ys)

        for i in range(self.nVertices):
            if minx < maxx:
                xs[i] = (xs[i] - minx) / (maxx - minx) * 0.9 + 0.05
            else:
                xs[i] = 0.5

            if miny < maxy:
                ys[i] = (ys[i] - miny) / (maxy - miny) * 0.9 + 0.05
            else:
                ys[i] = 0.5

        self.max_extent = max(maxx - minx, maxy - miny)

        return [(xs[i], ys[i]) for i in range(self.nVertices)]

    def print_timers(self) -> None:
        r"""
        Print the average time, in milliseconds, needed to compute each force and the correctness check.

        NOTE:
            Debug function, only used to estimate which methods should be optimized first. Called when the P key is hit.
        """
        if self.time_profile and self.frame > 0:
            print("Average correctness check time:", int(
                self.check_pos_correct_time * 1000 / self.frame))
            print("Average forces computing time:", int(
                self.update_forces_time * 1000 / self.frame))

            for force in self.forces_to_compute:
                print('\t', force.__name__, '\t', int(
                    self.forces_times[force.__name__] * 1000 / self.frame))
            print()

    def onKey(self, event: matplotlib.backend_bases.KeyEvent) -> None:
        r"""
        Callback function called by matplotlib on key presses.
        """
        if event.key == ' ':
            self.doFrame = True

        if event.key == "enter":
            if self.anim_running:
                # self.anim.event_source.stop()
                self.anim_running = False
            else:
                # self.anim.event_source.start()
                self.anim_running = True

        if event.key == "p":
            self.print_timers()

        if event.key == "escape":
            plt.close()

    def computeRepulsionVVForces(self) -> None:
        r"""
        Compute the repulsion force between each pair of vertices.
        """
        for i in range(self.nVertices):
            for j in range(i + 1, self.nVertices):
                r = self.pos[j] - self.pos[i]
                if i < self.real_n_vertices and j < self.real_n_vertices:
                    coef = self.repulsionVVRRCoef
                elif i >= self.real_n_vertices and j >= self.real_n_vertices:
                    coef = self.repulsionVVWWCoef
                else:
                    coef = self.repulsionVVRWCoef

                # force applied to i (going away from j)
                force = coef * r.normalized() / r.normSq()
                self.forces[i] -= force
                self.forces[j] += force

    def computeSpringForces(self) -> None:
        r"""
        Compute the spring force along each edge.
        """
        for (i, j) in self.edges:
            r = self.pos[j] - self.pos[i]

            if i < self.real_n_vertices and j < self.real_n_vertices:
                length = self.springLength
            else:
                length = self.springLength / self.break_down_num

            if self.useLogSpring:
                force = self.springCoef * \
                    math.log(r.norm() / length) * r.normalized()
            else:
                force = self.springCoef * (r.norm() - length) * r.normalized()

            self.forces[i] += force
            self.forces[j] -= force

    def computeWeakSpringForces(self) -> None:
        r"""
        Compute the weak spring force pulling each vertex towards the origin.
        """
        for i in range(self.nVertices):
            r = self.pos[i]
            force = -self.weakSpringCoef * r.norm() * r.normalized()

            self.forces[i] += force

    def computeTorsionForces(self) -> None:
        r"""
        Compute the torsion force between each pair of consecutive edges.
        """
        for i in range(self.nVertices):
            if len(self.embedding[i]) > 1:
                for k in range(len(self.embedding[i])):
                    # simulates an angular spring force between u1 and u2
                    u1 = self.embedding[i][k]
                    u2 = self.embedding[i][(k + 1) % len(self.embedding[i])]

                    v1, v2 = self.pos[u1] - \
                        self.pos[i], self.pos[u2] - self.pos[i]

                    angle = abs(v1.angle_towards(v2))
                    goal = 2.0 * math.pi / len(self.embedding[i])

                    # energy = self.torsionCoef * \
                    #     math.tanh(
                    #         angle / self.torsionScale) ** self.torsionPower

                    # scale = -2.0 * self.torsionPower * 10 * energy / self.torsionScale / math.sinh(2.0 * angle / self.torsionScale)

                    scale = self.torsionCoef * \
                        abs(angle - goal) ** self.torsionPower
                    if angle < goal:
                        scale = -scale

                    force1 = v1.rotate90() * (scale / v1.normSq())
                    force2 = v2.rotate90() * (scale / v2.normSq())

                    self.forces[u1] += force1
                    self.forces[i] += force2 - force1
                    self.forces[u2] -= force2

    def computeNodeEdgeForce(self, node: int, e1: int, e2: int) -> None:
        r"""
        Compute the repulsion force between the vertex ``node`` and the edges ``e1`` and ``e2``.
        """
        if node in (e1, e2):
            return

        distdiff = (self.pos[node] - self.pos[e1]).norm() + (self.pos[node] -
                                                             self.pos[e2]).norm() - (self.pos[e2] - self.pos[e1]).norm()

        energy = (distdiff / self.springLength) ** (-self.repulsionVertexEdgePower)

        scale = self.repulsionVertexEdgeCoef * energy / distdiff

        e1ton = (self.pos[node] - self.pos[e1]).normalized()
        e2ton = (self.pos[node] - self.pos[e2]).normalized()
        e1toe2 = (self.pos[e2] - self.pos[e1]).normalized()

        self.forces[node] += (e1ton + e2ton) * scale
        self.forces[e1] -= (e1ton - e1toe2) * scale
        self.forces[e2] -= (e2ton + e1toe2) * scale

    def computeRepulsionEEForces(self) -> None:
        r"""
        Compute the repulsion force between each pair of edges.
        """
        for face in self.faces:
            for i in range(
                    len(face)):      # iterate over every pair (node, edge)
                for j in range(len(face)):
                    self.computeNodeEdgeForce(
                        face[i], face[j], face[(j + 1) % len(face)])

    def check_pos_correct(self, pos: list[Vector2D]) -> bool:
        r"""
        Returns True if the given vertices positions are correct, i.e. do not induce edge crossings and respect the map embedding.
        """
        for i in range(self.nVertices):
            if len(self.embedding) > 2:
                angle_sum = 0
                for k in range(len(self.embedding[i])):
                    # simulates an angular spring force between u1 and u2
                    u1 = self.embedding[i][k]
                    u2 = self.embedding[i][(k + 1) % len(self.embedding[i])]

                    # embedding is given clockwise
                    angle_sum += 2.0 * math.pi - \
                        (pos[u1] - pos[i]).angle_towards(pos[u2] - pos[i])

                # if the embedding is correct, angle_sum = 2pi
                # otherwise, it's way bigger (>= 4pi?)

                if angle_sum >= 2.1 * math.pi:  # account for float inaccuracies
                    # print ("bad!! bad angle (sum", angle_sum, "!)")
                    # print (pos[i])
                    # for k in self.embedding[i]:
                    #    print (pos[k])
                    return False

        # check that there is no edge crossing
        # it suffices to check that for each face, none of the edges cross each
        # other

        for face in self.faces:
            segments = []
            for i in range(len(face)):
                segments.append(
                    tuple(sorted((pos[face[i]], pos[face[(i + 1) % len(face)]]))))

            # remove all identical segments
            segments = list(set(segments))

            if check_polygon_intersecting(segments):
                # print ("bad!! intersecting")
                # print (segments)
                return False

        return True

    def update_forces(self) -> None:
        r"""
        Update self.forces with each force computed.
        """
        self.forces = [Vector2D() for i in range(self.nVertices)]

        for force in self.forces_to_compute:
            if self.time_profile:
                begin = time.perf_counter()
            force()
            if self.time_profile:
                end = time.perf_counter()
                self.forces_times[force.__name__] += end - begin  # type: ignore
                # error was '"begin" is possibly unbound', which cannot be the case here
        """
        self.computeRepulsionVVForces()
        self.computeSpringForces()
        self.computeWeakSpringForces()
        self.computeTorsionForces()
        self.computeRepulsionEEForces()
        """

    def tick(self) -> None:
        r"""
        Apply a single tick (compute forces, update positions, reduce delta_t if needed and repeat until the new disposition is valid).
        """
        if self.use_exponential_model:
            base_delta_t = self.default_delta_t * \
                math.exp(-self.frame / (self.exp_decrease * self.nEdges))
            iters = 1

        else:
            # after this number of iterations, we should be closer to the
            # optimum; start going slower to refine the positions
            full_force_frames = int(self.nVertices * 1.5)
            # after this number of iterations, we start doing several frame at
            # once with decreasing delta_t to reach the optimum
            several_iter_frames = int(self.nVertices * 3)

            # delta_t will be proportional to 1 / frame ^ power when this
            # number is reached

            if self.frame < full_force_frames:
                base_delta_t = self.fast_delta_t
            else:
                base_delta_t = self.slow_delta_t * full_force_frames**self.decrease_power / self.frame**self.decrease_power

            if self.frame < several_iter_frames:
                iters = 1
            else:
                # we do several iterations with geometrically decreasing
                # delta_t to reach the optimum without flickering
                iters = self.max_iter_tick

        delta_t = base_delta_t
        delta_t /= (self.current_penalty + 1)

        self.current_delta_t = delta_t

        while iters > 0:
            # self.current_penalty = 0
            # currentMaxDispl = self.maxDispl / math.sqrt(self.current_penalty + 1)
            # currentMaxDispl = self.maxDispl / (2 ** self.current_penalty)
            # delta_t = base_delta_t / (2 ** self.current_penalty)

            if self.time_profile:
                begin = time.perf_counter()
            self.update_forces()  # LONG ? Et si oui, quelle force ?
            if self.time_profile:
                end = time.perf_counter()
                self.update_forces_time += float(end - begin)  # type: ignore
                # error was '"begin" is possibly unbound', which cannot be the case here

            # maxForce = max(map(Vector2D.normSq, self.forces))
            # if maxForce > 0:
            #    delta_t = min(delta_t, 2.0 / maxForce)

            newpos = [Vector2D() for i in range(self.nVertices)]

            for i in range(self.nVertices):
                displ = self.forces[i] * delta_t
                if displ.normSq() > self.currentMaxDispl ** 2.0:
                    displ = displ.normalized() * self.currentMaxDispl

                newpos[i] = self.pos[i] + displ

                # considering speed (ie. assuming nodes have inertia) seems to make the convergence slower
                # self.pos[i] += self.speed[i] * self.maxDelta + 0.5 * self.forces[i] * self.maxDelta * self.maxDelta / self.vertexInertia
                # self.speed[i] /= 2              # friction forces (we need energy loss somewhere)
                # self.speed[i] += self.forces[i] * self.maxDelta / self.vertexInertia

            if self.time_profile:
                begin = time.perf_counter()

            if not self.is_planar or self.check_pos_correct(newpos):  # LONG ?
                self.pos = newpos
                iters -= 1
                delta_t /= 3.0

                self.frame += 1
                self.validFrames += 1
                if self.current_penalty > 0:
                    self.current_penalty -= 1
            else:
                # we're going too fast... reduce delta_t & maxDispl
                self.validFrames = 0
                self.current_penalty = (
                    self.current_penalty + self.wrong_move_penalty) * 2

                delta_t /= 2.0
                self.currentMaxDispl /= 2.0

                self.validFrames = 0

            if self.validFrames > self.validFramesThreshold:
                self.currentMaxDispl = 0.2 * self.maxDispl + 0.8 * self.currentMaxDispl

            if self.time_profile:
                end = time.perf_counter()
                self.check_pos_correct_time += float(end - begin)  # type: ignore
                # error was '"begin" is possibly unbound', which cannot be the case here

    def update_fig(self, frame) -> tuple[matplotlib.artist.Artist, ...]:
        r"""
        Update the matplotlib figure with the new positions, and call ``self.tick()`` if the animation is currently running.
        """
        if frame > 5 and (self.anim_running or self.doFrame):
            self.tick()
            self.doFrame = False
            # self.pos[0].x += (random.random() - 0.5) / 5
        else:
            self.current_delta_t = 0

        self.text.set_text("Frame: " + str(self.frame) + ("; done" if self.done else ""))

        # plt.axis("on")
        # plt.cla()

        # self.pos[0] += Vector2D((random.random()-.5), (random.random()-.5))

        # nx.draw(
        #     self.G,
        #     {i: (float(self.pos[i-1].x), float(self.pos[i-1].y)) for i in range(1, self.nVertices + 1)},
        #     ax=self.ax,
        #     labels={
        #         i: str(i)
        #         for i in range(1, self.nVertices + 1)
        #     },
        #     node_size=[300] * self.nVertices,
        #     nodelist=list(range(1, self.nVertices + 1)),
        #     arrows=False,
        #     with_labels=True,
        #     node_color="red",
        # )

        pos_centered = self.center_pos()

        if self.prev_pos_centered:
            dist = max((pos_centered[i][0] - self.prev_pos_centered[i][0]) ** 2.0 + (
                pos_centered[i][1] - self.prev_pos_centered[i][1]) ** 2.0 for i in range(self.nVertices))

            if dist < self.convergence_limit ** 2.0 and self.anim_running and not self.done and self.frame > 10:
                print("Convergence found")
                # self.anim.event_source.stop()
                self.anim_running = False
                self.done = True

        self.prev_pos_centered = pos_centered

        edges_pos = [(pos_centered[i], pos_centered[j])
                     for (i, j) in self.edges]

        self.nodes_plt.set_offsets(pos_centered)
        self.edges_plt.set_segments(edges_pos)

        ret = [self.nodes_plt, self.edges_plt, self.text]

        if self.show_halfedges:
            # labels_head = nx.draw_networkx_edge_labels(self.G, pos_centered, ax=self.ax, rotate=False, edge_labels=self.edge_labels_head, label_pos=0.7)
            # labels_tail = nx.draw_networkx_edge_labels(self.G, pos_centered, ax=self.ax, rotate=False, edge_labels=self.edge_labels_tail, label_pos=0.3)
            # labels_middle = nx.draw_networkx_edge_labels(self.G, pos_centered, ax=self.ax, rotate=False, edge_labels=self.edge_labels_middle, label_pos=0.5)

            # for (source, dest) in [(labels_head, self.labels_head), (labels_tail, self.labels_tail), (labels_middle, self.labels_middle)]:
            #    for (edge, txt) in source.items():
            #        dest[edge].set_position(txt.get_position())

            #        ret.append(dest[edge])

            # for (d, prop) in ((self.labels_head, 0.7),
            #                   (self.labels_tail, 0.3), (self.labels_middle, 0.5)):
            #     for (edge, txt) in d.items():
            #         txt.set_x(pos_centered[edge[0]][0] * prop +
            #                   pos_centered[edge[1]][0] * (1 - prop))
            #         txt.set_y(pos_centered[edge[0]][1] * prop +
            #                   pos_centered[edge[1]][1] * (1 - prop))

            #         ret.append(txt)

            for (txt, i, j, prop) in self.edge_texts:
                txt.set_x(pos_centered[i][0] * prop +
                          pos_centered[j][0] * (1 - prop))
                txt.set_y(pos_centered[i][1] * prop +
                          pos_centered[j][1] * (1 - prop))
                ret.append(txt)

        return tuple(ret)
