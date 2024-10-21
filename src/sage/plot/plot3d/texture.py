# sage_setup: distribution = sagemath-plot
r"""
Texture support

This module provides texture/material support for 3D Graphics
objects and plotting.  This is a very rough common interface for
Tachyon, x3d, and obj (mtl).   See :class:`Texture
<sage.plot.plot3d.texture.Texture>` for full details about options and use.

Initially, we have no textures set::

    sage: sage.plot.plot3d.base.Graphics3d().texture_set()
    set()

However, one can access these textures in the following manner::

    sage: G = tetrahedron(color='red') + tetrahedron(color='yellow') + tetrahedron(color='red', opacity=0.5)
    sage: [t for t in G.texture_set() if t.color == colors.red] # we should have two red textures
    [Texture(texture..., red, ff0000), Texture(texture..., red, ff0000)]
    sage: [t for t in G.texture_set() if t.color == colors.yellow] # ...and one yellow
    [Texture(texture..., yellow, ffff00)]

And the Texture objects keep track of all their data::

    sage: T = tetrahedron(color='red', opacity=0.5)
    sage: t = T.get_texture()
    sage: t.opacity
    0.5
    sage: T # should be translucent
    Graphics3d Object

AUTHOR:

- Robert Bradshaw (2007-07-07) Initial version.
"""
from textwrap import dedent

from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from sage.misc.fast_methods import WithEqualityById
from sage.structure.sage_object import SageObject

from sage.plot.colors import colors, Color


uniq_c = 0


def _new_global_texture_id():
    """
    Generate a new unique id for a texture.

    EXAMPLES::

        sage: sage.plot.plot3d.texture._new_global_texture_id()
        'texture...'
        sage: sage.plot.plot3d.texture._new_global_texture_id()
        'texture...'
    """
    global uniq_c
    uniq_c += 1
    return "texture%s" % uniq_c


def is_Texture(x):
    r"""
    Deprecated. Use ``isinstance(x, Texture)`` instead.

    EXAMPLES::

        sage: from sage.plot.plot3d.texture import is_Texture, Texture
        sage: t = Texture(0.5)
        sage: is_Texture(t)
        doctest:...: DeprecationWarning: Please use isinstance(x, Texture)
        See https://github.com/sagemath/sage/issues/27593 for details.
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(27593, "Please use isinstance(x, Texture)")
    return isinstance(x, Texture)


def parse_color(info, base=None):
    r"""
    Parse the color.

    It transforms a valid color string into a color object and
    a color object into an RBG tuple of length 3. Otherwise,
    it multiplies the info by the base color.

    INPUT:

    - ``info`` -- color, valid color str or number
    - ``base`` -- tuple of length 3 (default: ``None``)

    OUTPUT: a tuple or color

    EXAMPLES:

    From a color::

        sage: from sage.plot.plot3d.texture import parse_color
        sage: c = Color('red')
        sage: parse_color(c)
        (1.0, 0.0, 0.0)

    From a valid color str::

        sage: parse_color('red')
        RGB color (1.0, 0.0, 0.0)
        sage: parse_color('#ff0000')
        RGB color (1.0, 0.0, 0.0)

    From a non valid color str::

        sage: parse_color('redd')
        Traceback (most recent call last):
        ...
        ValueError: unknown color 'redd'

    From an info and a base::

        sage: opacity = 10
        sage: parse_color(opacity, base=(.2,.3,.4))
        (2.0, 3.0, 4.0)
    """
    if isinstance(info, Color):
        return info.rgb()
    elif isinstance(info, str):
        try:
            return Color(info)
        except KeyError:
            raise ValueError("unknown color '%s'" % info)
    else:
        r, g, b = base
        # We don't want to lose the data when we split it into its respective components.
        if not r:
            r = 1e-5
        if not g:
            g = 1e-5
        if not b:
            b = 1e-5
        return (float(info * r), float(info * g), float(info * b))


class Texture(WithEqualityById, SageObject, metaclass=ClasscallMetaclass):
    r"""
    Class representing a texture.

    See documentation of :meth:`Texture.__classcall__
    <sage.plot.plot3d.texture.Texture.__classcall__>` for more details and
    examples.

    EXAMPLES:

    We create a translucent texture::

        sage: from sage.plot.plot3d.texture import Texture
        sage: t = Texture(opacity=0.6)
        sage: t
        Texture(texture..., 6666ff)
        sage: t.opacity
        0.6
        sage: t.jmol_str('obj')
        'color obj translucent 0.4 [102,102,255]'
        sage: t.mtl_str()
        'newmtl texture...\nKa 0.2 0.2 0.5\nKd 0.4 0.4 1.0\nKs 0.0 0.0 0.0\nillum 1\nNs 1.0\nd 0.6'
        sage: t.x3d_str()
        "<Appearance><Material diffuseColor='0.4 0.4 1.0' shininess='1.0' specularColor='0.0 0.0 0.0'/></Appearance>"

    TESTS::

        sage: Texture(opacity=1/3).opacity
        0.3333333333333333

        sage: hash(Texture()) # random
        42
    """
    @staticmethod
    def __classcall__(cls, id=None, **kwds):
        r"""
        Construct a new texture by id.

        INPUT:

        - ``id`` -- a texture (default: ``None``), a dictionary, a color, a
          string, a tuple, ``None`` or any other type acting as an ID. If ``id`` is
          ``None`` and keyword ``texture`` is empty, then it returns a unique texture object.
        - ``texture`` -- a texture
        - ``color`` -- tuple or string (default: (.4, .4, 1))
        - ``opacity`` -- number between 0 and 1 (default: 1)
        - ``ambient`` -- number (default: 0.5)
        - ``diffuse`` -- number (default: 1)
        - ``specular`` -- number (default: 0)
        - ``shininess`` -- number (default: 1)
        - ``name`` -- string (default: ``None``)
        - ``**kwds`` -- other valid keywords

        OUTPUT: a texture object

        EXAMPLES:

        Texture from integer ``id``::

            sage: from sage.plot.plot3d.texture import Texture
            sage: Texture(17)
            Texture(17, 6666ff)

        Texture from rational ``id``::

            sage: Texture(3/4)
            Texture(3/4, 6666ff)

        Texture from a dict::

            sage: Texture({'color':'orange','opacity':0.5})
            Texture(texture..., orange, ffa500)

        Texture from a color::

            sage: c = Color('red')
            sage: Texture(c)
            Texture(texture..., ff0000)

        Texture from a valid string color::

            sage: Texture('red')
            Texture(texture..., red, ff0000)

        Texture from a non valid string color::

            sage: Texture('redd')
            Texture(redd, 6666ff)

        Texture from a tuple::

            sage: Texture((.2,.3,.4))
            Texture(texture..., 334c66)

        Now accepting negative arguments, reduced modulo 1::

            sage: Texture((-3/8, 1/2, 3/8))
            Texture(texture..., 9f7f5f)

        Textures using other keywords::

            sage: Texture(specular=0.4)
            Texture(texture..., 6666ff)
            sage: Texture(diffuse=0.4)
            Texture(texture..., 6666ff)
            sage: Texture(shininess=0.3)
            Texture(texture..., 6666ff)
            sage: Texture(ambient=0.7)
            Texture(texture..., 6666ff)
        """
        if isinstance(id, Texture):
            return id
        if 'texture' in kwds:
            t = kwds['texture']
            if isinstance(t, Texture):
                return t
            else:
                raise TypeError("texture keyword must be a texture object")
        if isinstance(id, dict):
            kwds = id
            if 'rgbcolor' in kwds:
                kwds['color'] = kwds['rgbcolor']
            id = None
        elif isinstance(id, Color):
            kwds['color'] = id.rgb()
            id = None
        elif isinstance(id, str) and id in colors:
            kwds['color'] = id
            id = None
        elif isinstance(id, tuple):
            kwds['color'] = id
            id = None
        if id is None:
            id = _new_global_texture_id()
        return typecall(cls, id, **kwds)

    def __init__(self, id, color=(.4, .4, 1), opacity=1, ambient=0.5,
                 diffuse=1, specular=0, shininess=1, name=None, **kwds):
        r"""
        Construction of a texture.

        See documentation of :meth:`Texture.__classcall__
        <sage.plot.plot3d.texture.Texture.__classcall__>` for more details and
        examples.

        EXAMPLES::

            sage: from sage.plot.plot3d.texture import Texture
            sage: Texture(3, opacity=0.6)
            Texture(3, 6666ff)
        """
        self.id = id
        if name is None and isinstance(color, str):
            name = color
        self.name = name

        if not isinstance(color, tuple):
            color = parse_color(color)
        else:
            if len(color) == 4:
                opacity = color[3]
            color = tuple(float(1) if c == 1 else float(c) % 1
                          for c in color[0: 3])

        self.color = color
        self.opacity = float(opacity)
        self.shininess = float(shininess)

        if not isinstance(ambient, tuple):
            ambient = parse_color(ambient, color)
        self.ambient = ambient

        if not isinstance(diffuse, tuple):
            diffuse = parse_color(diffuse, color)
        self.diffuse = diffuse

        if not isinstance(specular, tuple):
            specular = parse_color(specular, color)
        self.specular = specular

    def _repr_(self):
        """
        Return a string representation of the Texture object.

        EXAMPLES::

            sage: from sage.plot.plot3d.texture import Texture
            sage: Texture('yellow')               #indirect doctest
            Texture(texture..., yellow, ffff00)
            sage: Texture((1,1,0), opacity=.5)
            Texture(texture..., ffff00)
        """
        if self.name is not None:
            return f"Texture({self.id}, {self.name}, {self.hex_rgb()})"
        else:
            return f"Texture({self.id}, {self.hex_rgb()})"

    def hex_rgb(self):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.texture import Texture
            sage: Texture('red').hex_rgb()
            'ff0000'
            sage: Texture((1, .5, 0)).hex_rgb()
            'ff7f00'
        """
        return "{:02x}{:02x}{:02x}".format(*tuple(int(255 * s) for s in self.color))

    def tachyon_str(self):
        r"""
        Convert Texture object to string suitable for Tachyon ray tracer.

        EXAMPLES::

            sage: from sage.plot.plot3d.texture import Texture
            sage: t = Texture(opacity=0.6)
            sage: t.tachyon_str()
            'Texdef texture...\n  Ambient 0.3333333333333333 Diffuse 0.6666666666666666 Specular 0.0 Opacity 0.6\n   Color 0.4 0.4 1.0\n   TexFunc 0'
        """

        total_ambient = sum(self.ambient)
        total_diffuse = sum(self.diffuse)
        total_specular = sum(self.specular)

        total_color = total_ambient + total_diffuse + total_specular

        if total_color == 0:
            total_color = 1

        ambient = total_ambient / total_color
        diffuse = total_diffuse / total_color
        specular = total_specular / total_color

        return dedent("""\
            Texdef {id}
              Ambient {ambient!r} Diffuse {diffuse!r} Specular {specular!r} Opacity {opacity!r}
              Color {color[0]!r} {color[1]!r} {color[2]!r}
              TexFunc 0""").format(id=self.id, ambient=ambient,
                                   diffuse=diffuse, specular=specular,
                                   opacity=self.opacity, color=self.color)

    def x3d_str(self):
        r"""
        Convert Texture object to string suitable for x3d.

        EXAMPLES::

            sage: from sage.plot.plot3d.texture import Texture
            sage: t = Texture(opacity=0.6)
            sage: t.x3d_str()
            "<Appearance><Material diffuseColor='0.4 0.4 1.0' shininess='1.0' specularColor='0.0 0.0 0.0'/></Appearance>"
        """
        return (
            "<Appearance>"
            "<Material diffuseColor='{color[0]!r} {color[1]!r} {color[2]!r}' "
            "shininess='{shininess!r}' "
            "specularColor='{specular!r} {specular!r} {specular!r}'/>"
            "</Appearance>").format(color=self.color, shininess=self.shininess,
                                    specular=self.specular[0])

    def mtl_str(self):
        r"""
        Convert Texture object to string suitable for mtl output.

        EXAMPLES::

            sage: from sage.plot.plot3d.texture import Texture
            sage: t = Texture(opacity=0.6)
            sage: t.mtl_str()
            'newmtl texture...\nKa 0.2 0.2 0.5\nKd 0.4 0.4 1.0\nKs 0.0 0.0 0.0\nillum 1\nNs 1.0\nd 0.6'
        """
        return dedent("""\
            newmtl {id}
            Ka {ambient[0]!r} {ambient[1]!r} {ambient[2]!r}
            Kd {diffuse[0]!r} {diffuse[1]!r} {diffuse[2]!r}
            Ks {specular[0]!r} {specular[1]!r} {specular[2]!r}
            illum {illumination}
            Ns {shininess!r}
            d {opacity!r}"""
        ).format(id=self.id, ambient=self.ambient, diffuse=self.diffuse,
                 specular=self.specular,
                 illumination=(2 if sum(self.specular) > 0 else 1),
                 shininess=self.shininess, opacity=self.opacity)

    def jmol_str(self, obj):
        r"""
        Convert Texture object to string suitable for Jmol applet.

        INPUT:

        - ``obj`` -- str

        EXAMPLES::

            sage: from sage.plot.plot3d.texture import Texture
            sage: t = Texture(opacity=0.6)
            sage: t.jmol_str('obj')
            'color obj translucent 0.4 [102,102,255]'

        ::

            sage: sum([dodecahedron(center=[2.5*x, 0, 0], color=(1, 0, 0, x/10)) for x in range(11)]).show(aspect_ratio=[1,1,1], frame=False, zoom=2)
        """
        translucent = "translucent %s" % float(1 - self.opacity) if self.opacity < 1 else ""
        return "color {} {} [{},{},{}]".format(obj, translucent,
                                           int(255 * self.color[0]),
                                           int(255 * self.color[1]),
                                           int(255 * self.color[2]))
