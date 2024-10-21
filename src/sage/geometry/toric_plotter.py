# sage_setup: distribution = sagemath-polyhedra
r"""
Toric plotter

This module provides a helper class :class:`ToricPlotter` for producing plots
of objects related to toric geometry. Default plotting objects can be adjusted
using :func:`options` and reset using :func:`reset_options`.

AUTHORS:

- Andrey Novoseltsev (2010-10-03): initial version, using some code bits by
  Volker Braun.

EXAMPLES:

In most cases, this module is used indirectly, e.g. ::

    sage: fan = toric_varieties.dP6().fan()                                             # needs palp sage.graphs
    sage: fan.plot()                                                                    # needs palp sage.graphs sage.plot
    Graphics object consisting of 31 graphics primitives

You may change default plotting options as follows::

    sage: toric_plotter.options("show_rays")
    True
    sage: toric_plotter.options(show_rays=False)
    sage: toric_plotter.options("show_rays")
    False
    sage: fan.plot()                                                                    # needs palp sage.graphs sage.plot
    Graphics object consisting of 19 graphics primitives
    sage: toric_plotter.reset_options()
    sage: toric_plotter.options("show_rays")
    True
    sage: fan.plot()                                                                    # needs palp sage.graphs sage.plot
    Graphics object consisting of 31 graphics primitives
"""


# ****************************************************************************
#       Copyright (C) 2010 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2010 Andrey Novoseltsev <novoselt@gmail.com>
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import copy
from math import pi

from sage.arith.misc import integer_ceil as ceil, integer_floor as floor
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.misc.lazy_import import lazy_import
from sage.modules.free_module_element import vector
lazy_import("sage.plot.all", ["Color", "Graphics",
                              "arrow", "disk", "line", "point",
                              "polygon", "rainbow", "text"])
lazy_import("sage.plot.plot3d.all", "text3d")
from sage.rings.real_double import RDF
from sage.structure.sage_object import SageObject


# These options are used to initialize/reset plotting options.
# Most of them are set to "None" and "real default values" are computed
# automatically based on the plotted object and parameters actually provided by
# the user.
_default_options = dict()
_default_options["mode"] = "round" # Can be also "box" and "generators"
_default_options["show_lattice"] = None # Default is "True for small plots"
_default_options["show_rays"] = True
_default_options["show_generators"] = True
_default_options["show_walls"] = True

_default_options["generator_color"] = "blue"
_default_options["label_color"] = "black"
_default_options["point_color"] = "black"
_default_options["ray_color"] = "purple"
_default_options["wall_color"] = "rainbow"
_default_options["wall_alpha"] = 0.4

_default_options["point_size"] = None
_default_options["ray_thickness"] = 3
_default_options["generator_thickness"] = None
_default_options["font_size"] = 14

_default_options["ray_label"] = "u"
_default_options["wall_label"] = r"\sigma"

# If none of these are given, the default will be 2.5
_default_options["radius"] = None
_default_options["xmin"] = None
_default_options["xmax"] = None
_default_options["ymin"] = None
_default_options["ymax"] = None
_default_options["zmin"] = None
_default_options["zmax"] = None

_default_options["lattice_filter"] = None

_default_options["wall_zorder"] = -5
_default_options["ray_zorder"] = -4
_default_options["generator_zorder"] = -3
_default_options["point_zorder"] = -2
_default_options["label_zorder"] = -1

# These options are actually used as "current defaults" in plotting functions.
_options = copy(_default_options)


class ToricPlotter(SageObject):
    r"""
    Create a toric plotter.

    INPUT:

    - ``all_options`` -- a :class:`dictionary <dict>`, containing any of the
      options related to toric objects (see :func:`options`) and any other
      options that will be passed to lower level plotting functions

    - ``dimension`` -- integer (1, 2, or 3); dimension of toric objects to
      be plotted

    - ``generators`` -- (optional) a list of ray generators; see examples for
      a detailed explanation of this argument

    OUTPUT: a toric plotter

    EXAMPLES:

    In most cases there is no need to create and use :class:`ToricPlotter`
    directly. Instead, use plotting method of the object which you want to
    plot, e.g. ::

        sage: fan = toric_varieties.dP6().fan()                                         # needs palp sage.graphs
        sage: fan.plot()                                                                # needs palp sage.graphs sage.plot
        Graphics object consisting of 31 graphics primitives
        sage: print(fan.plot())                                                         # needs palp sage.graphs sage.plot
        Graphics object consisting of 31 graphics primitives

    If you do want to create your own plotting function for some toric
    structure, the anticipated usage of toric plotters is the following:

    - collect all necessary options in a dictionary;

    - pass these options and ``dimension`` to :class:`ToricPlotter`;

    - call :meth:`include_points` on ray generators and any other points that
      you want to be present on the plot (it will try to set appropriate
      cut-off bounds);

    - call :meth:`adjust_options` to choose "nice" default values for all
      options that were not set yet and ensure consistency of rectangular and
      spherical cut-off bounds;

    - call :meth:`set_rays` on ray generators to scale them to the cut-off
      bounds of the plot;

    - call appropriate ``plot_*`` functions to actually construct the plot.

    For example, the plot from the previous example can be obtained as
    follows::

        sage: # needs palp sage.graphs sage.plot
        sage: from sage.geometry.toric_plotter import ToricPlotter
        sage: options = dict()  # use default for everything
        sage: tp = ToricPlotter(options, fan.lattice().degree())
        sage: tp.include_points(fan.rays())
        sage: tp.adjust_options()
        sage: tp.set_rays(fan.rays())
        sage: result = tp.plot_lattice()
        sage: result += tp.plot_rays()
        sage: result += tp.plot_generators()
        sage: result += tp.plot_walls(fan(2))
        sage: result
        Graphics object consisting of 31 graphics primitives

    In most situations it is only necessary to include generators of rays, in
    this case they can be passed to the constructor as an optional argument.
    In the example above, the toric plotter can be completely set up using ::

        sage: tp = ToricPlotter(options, fan.lattice().degree(), fan.rays())            # needs palp sage.graphs sage.plot

    All options are exposed as attributes of toric plotters and can be modified
    after constructions, however you will have to manually call
    :meth:`adjust_options` and :meth:`set_rays` again if you decide to change
    the plotting mode and/or cut-off bounds. Otherwise plots maybe invalid.
    """

    def __init__(self, all_options, dimension, generators=None):
        r"""
        See :class:`ToricPlotter` for documentation.

        TESTS::

            sage: from sage.geometry.toric_plotter import ToricPlotter
            sage: tp = ToricPlotter(dict(), 2)
            sage: TestSuite(tp).run()
        """
        super().__init__()
        sd = self.__dict__
        extra_options = dict()
        self.extra_options = extra_options
        toric_options = options()
        for option, value in all_options.items():
            if option in toric_options:
                sd[option] = value
            else:
                extra_options[option] = value
        for option, value in toric_options.items():
            if option not in sd:
                sd[option] = value
        if dimension not in [1, 2, 3]:
            raise ValueError("toric objects can be plotted only for "
                             "dimensions 1, 2, and 3, not %s!" % dimension)
        self.dimension = dimension
        self.origin = vector(RDF, max(dimension, 2)) # 1-d is plotted in 2-d
        if self.mode not in ["box", "generators", "round"]:
            raise ValueError("unrecognized plotting mode: %s!" % self.mode)
        # If radius was explicitly set by the user, it sets other bounds too.
        # If we don't take it into account here, they will be replaced by
        # automatically computed values.
        if sd["radius"] is not None:
            for key in ["xmin", "ymin", "zmin"]:
                if sd[key] is None:
                    sd[key] = - sd["radius"]
            for key in ["xmax", "ymax", "zmax"]:
                if sd[key] is None:
                    sd[key] = sd["radius"]
        # We also set some of the "extra_options" if they were not given.
        if "axes" not in extra_options:
            extra_options["axes"] = False
        if "frame" not in extra_options:
            extra_options["frame"] = False
        if "aspect_ratio" not in extra_options:
            extra_options["aspect_ratio"] = 1
        if generators is not None:
            # Completely prepare the plotter
            self.include_points(generators)
            self.adjust_options()
            self.set_rays(generators)

    def __eq__(self, other):
        r"""
        Check if ``self`` is equal to ``other``.

        INPUT:

        - ``other`` -- anything

        OUTPUT: ``True`` if ``self`` is equal to ``other``, ``False`` otherwise

        TESTS::

            sage: from sage.geometry.toric_plotter import ToricPlotter
            sage: ToricPlotter(dict(), 2) == ToricPlotter(dict(), 2)
            True
            sage: ToricPlotter(dict(), 2) == 0
            False
        """
        # Just to make TestSuite happy...
        return type(self) is type(other) and self.__dict__ == other.__dict__

    def adjust_options(self):
        r"""
        Adjust plotting options.

        This function determines appropriate default values for those options,
        that were not specified by the user, based on the other options. See
        :class:`ToricPlotter` for a detailed example.

        OUTPUT: none

        TESTS::

            sage: from sage.geometry.toric_plotter import ToricPlotter
            sage: tp = ToricPlotter(dict(), 2)
            sage: print(tp.show_lattice)
            None
            sage: tp.adjust_options()
            sage: print(tp.show_lattice)
            True
        """
        # Unfortunately, some of the defaults are dimension specific for no
        # good reason but to remedy 2-d/3-d plotting inconsistencies in Sage.
        d = self.dimension
        if d <= 2:
            if self.point_size is None:
                self.point_size = 50
        elif d == 3:
            if self.point_size is None:
                self.point_size = 10
        if self.generator_thickness is None:
            self.generator_thickness = self.ray_thickness
        sd = self.__dict__
        bounds = ["radius", "xmin", "xmax", "ymin", "ymax", "zmin", "zmax"]
        bounds = [abs(sd[bound]) for bound in bounds if sd[bound] is not None]
        r = RDF(max(bounds + [0.5]) if bounds else 2.5)
        self.radius = r
        round = self.mode == "round"
        for key in ["xmin", "ymin", "zmin"]:
            if round or sd[key] is None:
                sd[key] = - r
            sd[key] = min(sd[key], - 0.5)
            sd[key] = RDF(sd[key])
        for key in ["xmax", "ymax", "zmax"]:
            if round or sd[key] is None:
                sd[key] = r
            sd[key] = max(sd[key], 0.5)
            sd[key] = RDF(sd[key])
        if self.show_lattice is None:
            self.show_lattice = (r <= 5) if d <= 2 else r <= 3

    def include_points(self, points, force=False):
        r"""
        Try to include ``points`` into the bounding box of ``self``.

        INPUT:

        - ``points`` -- list of points

        - ``force`` -- boolean (default: ``False``); by default, only bounds
          that were not set before will be chosen to include ``points``. Use
          ``force=True`` if you don't mind increasing existing bounding box.

        OUTPUT: none

        EXAMPLES::

            sage: from sage.geometry.toric_plotter import ToricPlotter
            sage: tp = ToricPlotter(dict(), 2)
            sage: print(tp.radius)
            None
            sage: tp.include_points([(3, 4)])
            sage: print(tp.radius)
            5.5...
            sage: tp.include_points([(5, 12)])
            sage: print(tp.radius)
            5.5...
            sage: tp.include_points([(5, 12)], force=True)
            sage: print(tp.radius)
            13.5...
        """
        if not points:
            return
        points = [vector(RDF, pt) for pt in points]
        sd = self.__dict__

        def update(bound, new_value, points):
            if force or sd[bound] is None:
                new_value = eval(new_value)
                if sd[bound] is None:
                    sd[bound] = new_value
                elif abs(sd[bound]) < abs(new_value):
                    sd[bound] = new_value

        update("radius", "max(pt.norm() for pt in points) + 0.5", points)
        try:
            update("xmin", "min(pt[0] for pt in points) - 0.5", points)
            update("xmax", "max(pt[0] for pt in points) + 0.5", points)
            update("ymin", "min(pt[1] for pt in points) - 0.5", points)
            update("ymax", "max(pt[1] for pt in points) + 0.5", points)
            update("zmin", "min(pt[2] for pt in points) - 0.5", points)
            update("zmax", "max(pt[2] for pt in points) + 0.5", points)
        except IndexError:  # 1 or 2 dimensional case
            pass

    def plot_generators(self):
        r"""
        Plot ray generators.

        Ray generators must be specified during construction or using
        :meth:`set_rays` before calling this method.

        OUTPUT: a plot

        EXAMPLES::

            sage: from sage.geometry.toric_plotter import ToricPlotter
            sage: tp = ToricPlotter(dict(), 2, [(3,4)])
            sage: tp.plot_generators()                                                  # needs sage.plot
            Graphics object consisting of 1 graphics primitive
        """
        generators = self.generators
        result = Graphics()
        if not generators or not self.show_generators:
            return result
        colors = color_list(self.generator_color, len(generators))
        d = self.dimension
        extra_options = self.extra_options
        origin = self.origin
        thickness = self.generator_thickness
        zorder = self.generator_zorder
        for generator, ray, color in zip(generators, self.rays, colors):
            if ray.dot_product(ray) < generator.dot_product(generator):
                result += line([origin, ray],
                               color=color, thickness=thickness,
                               zorder=zorder, **extra_options)
            else:
                # This should not be the case, but as of 4.6 plotting
                # functions are inconsistent and arrows behave very
                # different compared to lines.
                if d <= 2:
                    result += arrow(origin, generator,
                                    color=color, width=thickness,
                                    arrowsize=thickness + 1,
                                    zorder=zorder, **extra_options)
                else:
                    result += line([origin, generator], arrow_head=True,
                                   color=color, thickness=thickness,
                                   zorder=zorder, **extra_options)
        return result

    def plot_labels(self, labels, positions):
        r"""
        Plot ``labels`` at specified ``positions``.

        INPUT:

        - ``labels`` -- string or list of strings

        - ``positions`` -- list of points

        OUTPUT: a plot

        EXAMPLES::

            sage: from sage.geometry.toric_plotter import ToricPlotter
            sage: tp = ToricPlotter(dict(), 2)
            sage: tp.plot_labels("u", [(1.5,0)])                                        # needs sage.plot
            Graphics object consisting of 1 graphics primitive
        """
        result = Graphics()
        color = self.label_color
        extra_options = self.extra_options
        zorder = self.label_zorder
        font_size = self.font_size
        twod = self.dimension <= 2
        labels = label_list(labels, len(positions), twod)
        for label, position in zip(labels, positions):
            if label is None:
                continue
            if twod:
                result += text(label, position,
                               color=color, fontsize=font_size,
                               zorder=zorder, **extra_options)
            else:
                result += text3d(label, position, color=color, **extra_options)
        return result

    def plot_lattice(self):
        r"""
        Plot the lattice (i.e. its points in the cut-off bounds of ``self``).

        OUTPUT: a plot

        EXAMPLES::

            sage: from sage.geometry.toric_plotter import ToricPlotter
            sage: tp = ToricPlotter(dict(), 2)
            sage: tp.adjust_options()
            sage: tp.plot_lattice()                                                     # needs sage.plot
            Graphics object consisting of 1 graphics primitive
        """
        if not self.show_lattice:
            # Plot the origin anyway, otherwise rays/generators may look ugly.
            return self.plot_points([self.origin])
        d = self.dimension
        if d == 1:
            points = ((x, 0)
                      for x in range(ceil(self.xmin), floor(self.xmax) + 1))
        elif d == 2:
            points = ((x, y)
                      for x in range(ceil(self.xmin), floor(self.xmax) + 1)
                      for y in range(ceil(self.ymin), floor(self.ymax) + 1))
        elif d == 3:
            points = ((x, y, z)
                      for x in range(ceil(self.xmin), floor(self.xmax) + 1)
                      for y in range(ceil(self.ymin), floor(self.ymax) + 1)
                      for z in range(ceil(self.zmin), floor(self.zmax) + 1))
        if self.mode == "round":
            r = 1.01 * self.radius # To make sure integer values work OK.
            points = (pt for pt in points if vector(pt).dot_product(vector(pt)) <= r)
        f = self.lattice_filter
        if f is not None:
            points = (pt for pt in points if f(pt))
        return self.plot_points(tuple(points))

    def plot_points(self, points):
        r"""
        Plot given ``points``.

        INPUT:

        - ``points`` -- list of points

        OUTPUT: a plot

        EXAMPLES::

            sage: from sage.geometry.toric_plotter import ToricPlotter
            sage: tp = ToricPlotter(dict(), 2)
            sage: tp.adjust_options()
            sage: tp.plot_points([(1,0), (0,1)])                                        # needs sage.plot
            Graphics object consisting of 1 graphics primitive
        """
        return point(points, color=self.point_color, size=self.point_size,
                     zorder=self.point_zorder, **self.extra_options)

    def plot_ray_labels(self):
        r"""
        Plot ray labels.

        Usually ray labels are plotted together with rays, but in some cases it
        is desirable to output them separately.

        Ray generators must be specified during construction or using
        :meth:`set_rays` before calling this method.

        OUTPUT: a plot

        EXAMPLES::

            sage: from sage.geometry.toric_plotter import ToricPlotter
            sage: tp = ToricPlotter(dict(), 2, [(3,4)])
            sage: tp.plot_ray_labels()                                                  # needs sage.plot
            Graphics object consisting of 1 graphics primitive
        """
        return self.plot_labels(self.ray_label,
                                [1.1 * ray for ray in self.rays])

    def plot_rays(self):
        r"""
        Plot rays and their labels.

        Ray generators must be specified during construction or using
        :meth:`set_rays` before calling this method.

        OUTPUT: a plot

        EXAMPLES::

            sage: from sage.geometry.toric_plotter import ToricPlotter
            sage: tp = ToricPlotter(dict(), 2, [(3,4)])
            sage: tp.plot_rays()                                                        # needs sage.plot
            Graphics object consisting of 2 graphics primitives
        """
        result = Graphics()
        rays = self.rays
        if not rays or not self.show_rays:
            return result
        extra_options = self.extra_options
        origin = self.origin
        colors = color_list(self.ray_color, len(rays))
        thickness = self.ray_thickness
        zorder = self.ray_zorder
        for end, color in zip(rays, colors):
            result += line([origin, end],
                           color=color, thickness=thickness,
                           zorder=zorder, **extra_options)
        result += self.plot_ray_labels()
        return result

    def plot_walls(self, walls):
        r"""
        Plot ``walls``, i.e. 2-d cones, and their labels.

        Ray generators must be specified during construction or using
        :meth:`set_rays` before calling this method and these specified
        ray generators will be used in conjunction with
        :meth:`~sage.geometry.cone.ConvexRationalPolyhedralCone.ambient_ray_indices`
        of ``walls``.

        INPUT:

        - ``walls`` -- list of 2-d cones

        OUTPUT: a plot

        EXAMPLES::

            sage: quadrant = Cone([(1,0), (0,1)])
            sage: from sage.geometry.toric_plotter import ToricPlotter
            sage: tp = ToricPlotter(dict(), 2, quadrant.rays())
            sage: tp.plot_walls([quadrant])                                             # needs sage.plot
            Graphics object consisting of 2 graphics primitives

        Let's also check that the truncating polyhedron is functioning
        correctly::

            sage: tp = ToricPlotter({"mode": "box"}, 2, quadrant.rays())
            sage: tp.plot_walls([quadrant])                                             # needs sage.plot
            Graphics object consisting of 2 graphics primitives
        """
        result = Graphics()
        if not walls or not self.show_walls:
            return result
        rays = self.rays
        extra_options = self.extra_options
        mode = self.mode
        alpha = self.wall_alpha
        colors = color_list(self.wall_color, len(walls))
        zorder = self.wall_zorder
        if mode == "box":
            if self.dimension <= 2:
                ieqs = [(self.xmax, -1, 0), (- self.xmin, 1, 0),
                        (self.ymax, 0, -1), (- self.ymin, 0, 1)]
            else:
                ieqs = [(self.xmax, -1, 0, 0), (- self.xmin, 1, 0, 0),
                        (self.ymax, 0, -1, 0), (- self.ymin, 0, 1, 0),
                        (self.zmax, 0, 0, -1), (- self.zmin, 0, 0, 1)]
            box = Polyhedron(ieqs=ieqs, base_ring=RDF)
            for wall, color in zip(walls, colors):
                result += box.intersection(wall.polyhedron()).render_solid(
                    alpha=alpha, color=color, zorder=zorder, **extra_options)
        elif mode == "generators":
            origin = self.origin
            for wall, color in zip(walls, colors):
                vertices = [rays[i] for i in wall.ambient_ray_indices()]
                vertices.append(origin)
                result += Polyhedron(vertices=vertices, base_ring=RDF).render_solid(
                    alpha=alpha, color=color, zorder=zorder, **extra_options)
        label_sectors = []
        round = mode == "round"
        for wall, color in zip(walls, colors):
            S = wall.linear_subspace()
            lsd = S.dimension()
            if lsd == 0:    # Strictly convex wall
                r1, r2 = (rays[i] for i in wall.ambient_ray_indices())
            elif lsd == 1:  # wall is a half-plane
                for i, ray in zip(wall.ambient_ray_indices(), wall.rays()):
                    if ray in S:
                        r1 = rays[i]
                    else:
                        r2 = rays[i]
                if round:
                    # Plot one "extra" sector
                    result += sector(- r1, r2,
                      alpha=alpha, color=color, zorder=zorder, **extra_options)
            else:           # wall is a plane
                r1, r2 = S.basis()
                r1 = vector(RDF, r1)
                r1 = r1 / r1.norm() * self.radius
                r2 = vector(RDF, r2)
                r2 = r2 / r2.norm() * self.radius
                if round:
                    # Plot three "extra" sectors
                    result += sector(r1, - r2,
                      alpha=alpha, color=color, zorder=zorder, **extra_options)
                    result += sector(- r1, r2,
                      alpha=alpha, color=color, zorder=zorder, **extra_options)
                    result += sector(- r1, - r2,
                      alpha=alpha, color=color, zorder=zorder, **extra_options)
            label_sectors.append([r1, r2])
            if round:
                result += sector(r1, r2,
                    alpha=alpha, color=color, zorder=zorder, **extra_options)
        result += self.plot_labels(self.wall_label,
                    [sum(label_sector) / 3 for label_sector in label_sectors])
        return result

    def set_rays(self, generators):
        r"""
        Set up rays and their ``generators`` to be used by plotting functions.

        As an alternative to using this method, you can pass ``generators`` to
        :class:`ToricPlotter` constructor.

        INPUT:

        - ``generators`` -- list of primitive nonzero ray generators

        OUTPUT: none

        EXAMPLES::

            sage: from sage.geometry.toric_plotter import ToricPlotter
            sage: tp = ToricPlotter(dict(), 2)
            sage: tp.adjust_options()
            sage: tp.plot_rays()                                                        # needs sage.plot
            Traceback (most recent call last):
            ...
            AttributeError: 'ToricPlotter' object has no attribute 'rays'...
            sage: tp.set_rays([(0,1)])
            sage: tp.plot_rays()                                                        # needs sage.plot
            Graphics object consisting of 2 graphics primitives
        """
        d = self.dimension
        if d == 1:
            generators = [vector(RDF, 2, (gen[0], 0)) for gen in generators]
        else:
            generators = [vector(RDF, d, gen) for gen in generators]
        self.generators = generators
        if self.mode == "box":
            rays = []
            bounds = [self.__dict__[bound]
                for bound in ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"]]
            bounds = bounds[:2 * d]
            for gen in generators:
                factors = []
                for i, gen_i in enumerate(gen):
                    factors.append(gen_i / bounds[2 * i])
                    factors.append(gen_i / bounds[2 * i + 1])
                rays.append(gen / max(factors))
        elif self.mode == "generators":
            rays = generators
        elif self.mode == "round":
            r = self.radius
            rays = [r * gen / gen.norm() for gen in generators]
        self.rays = rays


def _unrecognized_option(option):
    r"""
    Raise an exception about wrong ``option``.

    INPUT:

    - ``option`` -- string

    OUTPUT: none, a :exc:`KeyError` exception is raised

    TESTS::

        sage: from sage.geometry.toric_plotter import _unrecognized_option
        sage: _unrecognized_option("nontoric")
        Traceback (most recent call last):
        ...
        KeyError: "unrecognized toric plot option: 'nontoric'!
        Type 'toric_plotter.options?' to see available options."
    """
    raise KeyError("unrecognized toric plot option: '%s'! " % option
                + "Type 'toric_plotter.options?' to see available options.")


def color_list(color, n):
    r"""
    Normalize a list of ``n`` colors.

    INPUT:

    - ``color`` -- anything specifying a :class:`Color`, a list of such
      specifications, or the string "rainbow";

    - ``n`` -- integer

    OUTPUT: list of ``n`` colors

    If ``color`` specified a single color, it is repeated ``n`` times. If it
    was a list of ``n`` colors, it is returned without changes. If it was
    "rainbow", the rainbow of ``n`` colors is returned.

    EXAMPLES::

        sage: # needs sage.plot
        sage: from sage.geometry.toric_plotter import color_list
        sage: color_list("grey", 1)
        [RGB color (0.5019607843137255, 0.5019607843137255, 0.5019607843137255)]
        sage: len(color_list("grey", 3))
        3
        sage: L = color_list("rainbow", 3)
        sage: L
        [RGB color (1.0, 0.0, 0.0),
         RGB color (0.0, 1.0, 0.0),
         RGB color (0.0, 0.0, 1.0)]
        sage: color_list(L, 3)
        [RGB color (1.0, 0.0, 0.0),
         RGB color (0.0, 1.0, 0.0),
         RGB color (0.0, 0.0, 1.0)]
        sage: color_list(L, 4)
        Traceback (most recent call last):
        ...
        ValueError: expected 4 colors, got 3!
    """
    try:
        color = Color(color)
        return [color] * n
    except (ValueError, TypeError):
        if isinstance(color, (list, tuple)):
            if len(color) != n:
                raise ValueError("expected %d colors, got %d!"
                                 % (n, len(color)))
            return color
        if color == "rainbow":
            return [Color(c) for c in rainbow(n, "rgbtuple")]
    raise TypeError("cannot interpret %s as a color!" % color)


def label_list(label, n, math_mode, index_set=None):
    r"""
    Normalize a list of ``n`` labels.

    INPUT:

    - ``label`` -- ``None``, a string, or a list of string

    - ``n`` -- integer

    - ``math_mode`` -- boolean; if ``True``, will produce LaTeX expressions
      for labels

    - ``index_set`` -- list of integers (default: ``range(n)``) that will be
      used as subscripts for labels

    OUTPUT: list of ``n`` labels

    If ``label`` was a list of ``n`` entries, it is returned without changes.
    If ``label`` is ``None``, a list of ``n`` ``None``'s is returned. If
    ``label`` is a string, a list of strings of the form ``$label_{i}$`` is
    returned, where `i` ranges over ``index_set``. (If ``math_mode=False``, the
    form "label_i" is used instead.) If ``n=1``, there is no subscript added,
    unless ``index_set`` was specified explicitly.

    EXAMPLES::

        sage: from sage.geometry.toric_plotter import label_list
        sage: label_list("u", 3, False)
        ['u_0', 'u_1', 'u_2']
        sage: label_list("u", 3, True)
        ['$u_{0}$', '$u_{1}$', '$u_{2}$']
        sage: label_list("u", 1, True)
        ['$u$']
    """
    if label is None:
        return [None] * n
    if isinstance(label, (list, tuple)):
        if len(label) != n:
            raise ValueError("expected %d labels, got %d!" % (n, len(label)))
        return label
    if index_set is None:
        if n == 1:
            return ["$%s$" % label.strip("$")] if math_mode else [label]
        index_set = range(n)
    if math_mode:
        label = label.strip("$")
        return list("$%s_{%d}$" % (label, i) for i in index_set)
    else:
        return list("%s_%d" % (label, i) for i in index_set)


def options(option=None, **kwds):
    r"""
    Get or set options for plots of toric geometry objects.

    .. NOTE::

        This function provides access to global default options. Any of these
        options can be overridden by passing them directly to plotting
        functions. See also :func:`reset_options`.

    INPUT:

    - None;

    OR:

    - ``option`` -- string, name of the option whose value you wish to get;

    OR:

    - keyword arguments specifying new values for one or more options.

    OUTPUT:

    - if there was no input, the dictionary of current options for toric plots;

    - if ``option`` argument was given, the current value of ``option``;

    - if other keyword arguments were given, none.

    **Name Conventions**

    To clearly distinguish parts of toric plots, in options and methods we use
    the following name conventions:

    Generator
        A primitive integral vector generating a 1-dimensional cone, plotted as
        an arrow from the origin (or a line, if the head of the arrow is beyond
        cut-off bounds for the plot).

    Ray
        A 1-dimensional cone, plotted as a line from the origin to the cut-off
        bounds for  the plot.

    Wall
        A 2-dimensional cone, plotted as a region between rays (in the above
        sense). Its exact shape depends on the plotting mode (see below).

    Chamber
        A 3-dimensional cone, plotting is not implemented yet.

    **Plotting Modes**

    A plotting mode mostly determines the shape of the cut-off region (which is
    always relevant for toric plots except for trivial objects consisting of
    the origin only). The following options are available:

    Box
        The cut-off region is a box with edges parallel to coordinate axes.

    Generators
        The cut-off region is determined by primitive integral generators of
        rays. Note that this notion is well-defined only for rays and walls,
        in particular you should plot the lattice on your own
        (:meth:`~ToricPlotter.plot_lattice` will use box mode which is likely
        to be unsuitable). While this method may not be suitable for general
        fans, it is quite natural for fans of :class:`CPR-Fano toric varieties.
        <sage.schemes.toric.fano_variety.CPRFanoToricVariety_field`

    Round
        The cut-off regions is a sphere centered at the origin.

    **Available Options**

    Default values for the following options can be set using this function:


    - ``mode`` -- "box", "generators", or "round", see above for descriptions;

    - ``show_lattice`` -- boolean, whether to show lattice points in the
      cut-off region or not;

    - ``show_rays`` -- boolean, whether to show rays or not;

    - ``show_generators`` -- boolean, whether to show rays or not;

    - ``show_walls`` -- boolean, whether to show rays or not;

    - ``generator_color`` -- a color for generators;

    - ``label_color`` -- a color for labels;

    - ``point_color`` -- a color for lattice points;

    - ``ray_color`` -- a color for rays, a list of colors (one for each ray),
      or the string "rainbow";

    - ``wall_color`` -- a color for walls, a list of colors (one for each
      wall), or the string "rainbow";

    - ``wall_alpha`` -- a number between 0 and 1, the alpha-value for walls
      (determining their transparency);

    - ``point_size`` -- integer; the size of lattice points

    - ``ray_thickness`` -- integer; the thickness of rays

    - ``generator_thickness`` -- integer; the thickness of generators

    - ``font_size`` -- integer; the size of font used for labels

    - ``ray_label`` -- string or list of strings used for ray labels; use
      ``None`` to hide labels

    - ``wall_label`` -- string or list of strings used for wall labels; use
      ``None`` to hide labels

    - ``radius`` -- a positive number, the radius of the cut-off region for
      "round" mode

    - ``xmin``, ``xmax``, ``ymin``, ``ymax``, ``zmin``, ``zmax`` -- numbers
      determining the cut-off region for "box" mode. Note that you cannot
      exclude the origin - if you try to do so, bounds will be automatically
      expanded to include it.

    - ``lattice_filter`` -- a callable, taking as an argument a lattice point
      and returning ``True`` if this point should be included on the plot
      (useful, e.g. for plotting sublattices)

    - ``wall_zorder``, ``ray_zorder``, ``generator_zorder``, ``point_zorder``,
      ``label_zorder`` -- integers, z-orders for different classes of objects.
      By default all values are negative, so that you can add other graphic
      objects on top of a toric plot. You may need to adjust these parameters
      if you want to put a toric plot on top of something else or if you want
      to overlap several toric plots.

    You can see the current default value of any options by typing, e.g. ::

        sage: toric_plotter.options("show_rays")
        True

    If the default value is ``None``, it means that the actual default is
    determined later based on the known options. Note, that not all options can
    be determined in such a way, so you should not set options to ``None``
    unless it was its original state. (You can always revert to this "original
    state" using :meth:`reset_options`.)

    EXAMPLES:

    The following line will make all subsequent toric plotting commands to draw
    "rainbows" from walls::

        sage: toric_plotter.options(wall_color='rainbow')

    If you prefer a less colorful output (e.g. if you need black-and-white
    illustrations for a paper), you can use something like this::

        sage: toric_plotter.options(wall_color='grey')
    """
    global _options
    if option is None and not kwds:
        return copy(_options)
    elif option is not None and not kwds:
        try:
            return _options[option]
        except KeyError:
            _unrecognized_option(option)
    elif option is None and kwds:
        for option in kwds:
            try:
                _options[option] = kwds[option]
            except KeyError:
                _unrecognized_option(option)
    else:
        raise ValueError("you cannot specify 'option' and other arguments at "
                         "the same time!")


def reset_options():
    r"""
    Reset options for plots of toric geometry objects.

    OUTPUT: none

    EXAMPLES::

        sage: toric_plotter.options("show_rays")
        True
        sage: toric_plotter.options(show_rays=False)
        sage: toric_plotter.options("show_rays")
        False

    Now all toric plots will not show rays, unless explicitly requested. If you
    want to go back to "default defaults", use this method::

        sage: toric_plotter.reset_options()
        sage: toric_plotter.options("show_rays")
        True
    """
    global _options
    _options = copy(_default_options)


def sector(ray1, ray2, **extra_options):
    r"""
    Plot a sector between ``ray1`` and ``ray2`` centered at the origin.

    .. NOTE::

        This function was intended for plotting strictly convex cones, so it
        plots the smaller sector between ``ray1`` and ``ray2`` and, therefore,
        they cannot be opposite. If you do want to use this function for bigger
        regions, split them into several parts.

    .. NOTE::

        As of version 4.6 Sage does not have a graphic primitive for sectors in
        3-dimensional space, so this function will actually approximate them
        using polygons (the number of vertices used depends on the angle
        between rays).

    INPUT:

    - ``ray1``, ``ray2`` -- rays in 2- or 3-dimensional space of the same
      length

    - ``extra_options`` -- dictionary of options that should be passed to
      lower level plotting functions

    OUTPUT: a plot

    EXAMPLES::

        sage: from sage.geometry.toric_plotter import sector
        sage: sector((1,0), (0,1))                                                      # needs sage.symbolic
        Graphics object consisting of 1 graphics primitive
        sage: sector((3,2,1), (1,2,3))                                                  # needs sage.plot
        Graphics3d Object
    """
    from sage.functions.all import arccos, arctan2

    ray1 = vector(RDF, ray1)
    ray2 = vector(RDF, ray2)
    r = ray1.norm()
    if len(ray1) == 2:
        # Plot an honest sector
        phi1 = arctan2(ray1[1], ray1[0])
        phi2 = arctan2(ray2[1], ray2[0])
        if phi1 > phi2:
            phi1, phi2 = phi2, phi1
        if phi2 - phi1 > pi:
            phi1, phi2 = phi2, phi1 + 2 * pi
        return disk((0,0), r, (phi1, phi2), **extra_options)
    else:
        # Plot a polygon, 30 vertices per radian.
        vertices_per_radian = 30
        n = ceil(arccos(ray1 * ray2 / r**2) * vertices_per_radian)
        dr = (ray2 - ray1) / n
        points = (ray1 + i * dr for i in range(n + 1))
        points = [r / pt.norm() * pt for pt in points]
        points.append(vector(RDF, 3))
        return polygon(points, **extra_options)
