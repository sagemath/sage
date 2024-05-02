
from .plot3d            import plot3d, cylindrical_plot3d, spherical_plot3d, Spherical, SphericalElevation, Cylindrical
from .parametric_plot3d import parametric_plot3d
from .plot_field3d      import plot_vector_field3d

# We lazy_import the following modules since they import numpy which slows down sage startup
from sage.misc.lazy_import import lazy_import
lazy_import("sage.plot.plot3d.implicit_plot3d",["implicit_plot3d"])

from .list_plot3d       import list_plot3d
from .revolution_plot3d import revolution_plot3d

from .platonic          import tetrahedron, cube, octahedron, dodecahedron, icosahedron

from .shapes2           import sphere, line3d, polygon3d, polygons3d, point3d, text3d, bezier3d

from .shapes            import arrow3d

#from shapes import Box, ColorCube, Cone, Cylinder, LineSegment, Arrow, Sphere, Torus, Text as Text3D
#from parametric_surface import ParametricSurface, MoebiusStrip
#from plot3d import plot3d, axes as axes3d
