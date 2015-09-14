r"""
Texture Support

This module provides texture/material support for 3D Graphics
objects and plotting.  This is a very rough common interface for
Tachyon, x3d, and obj (mtl).   See
:meth:`Texture <sage.plot.plot3d.texture.Texture>` and
:class:`Texture_class <sage.plot.plot3d.texture.Texture_class>`
for full details about options and use.

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
from sage.misc.fast_methods import WithEqualityById
from sage.structure.sage_object import SageObject

from sage.plot.colors import Color


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


from sage.plot.colors import colors

def is_Texture(x):
    r"""
    Return whether ``x`` is an instance of ``Texture_class``.

    EXAMPLES::

        sage: from sage.plot.plot3d.texture import is_Texture, Texture
        sage: t = Texture(0.5)
        sage: is_Texture(t)
        True

    ::

        sage: is_Texture(4)
        False
    """
    return isinstance(x, Texture_class)

def Texture(id=None, **kwds):
    r"""
    Return a texture.

    INPUT:

    - ``id`` - a texture (optional, default: None), a dict, a color, a
      str, a tuple, None or any other type acting as an ID. If ``id`` is
      None and keyword ``texture`` is empty, then it returns a unique texture object.
    - ``texture`` - a texture
    - ``color`` - tuple or str, (optional, default: (.4, .4, 1))
    - ``opacity`` - number between 0 and 1 (optional, default: 1)
    - ``ambient`` - number (optional, default: 0.5)
    - ``diffuse`` - number (optional, default: 1)
    - ``specular`` - number (optional, default: 0)
    - ``shininess`` - number (optional, default: 1)
    - ``name`` - str (optional, default: None)
    - ``**kwds`` - other valid keywords

    OUTPUT:

    A texture object.

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
    if isinstance(id, Texture_class):
        return id
    if 'texture' in kwds:
        t = kwds['texture']
        if is_Texture(t):
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
    return Texture_class(id, **kwds)

def parse_color(info, base=None):
    r"""
    Parses the color.

    It transforms a valid color string into a color object and
    a color object into an RBG tuple of length 3. Otherwise,
    it multiplies the info by the base color.

    INPUT:

    - ``info`` - color, valid color str or number
    - ``base`` - tuple of length 3 (optional, default: None)

    OUTPUT:

    A tuple or color.

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
            raise ValueError("unknown color '%s'"%info)
    else:
        r, g, b = base
        # We don't want to lose the data when we split it into its respective components.
        if not r: r = 1e-5
        if not g: g = 1e-5
        if not b: b = 1e-5
        return (float(info*r), float(info*g), float(info*b))


class Texture_class(WithEqualityById, SageObject):
    r"""
    Construction of a texture.

    See documentation of :meth:`Texture <sage.plot.plot3d.texture.Texture>`
    for more details and examples.

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
    def __init__(self, id, color=(.4, .4, 1), opacity=1, ambient=0.5, diffuse=1, specular=0, shininess=1, name=None, **kwds):
        r"""
        Construction of a texture.

        See documentation of :meth:`Texture <sage.plot.plot3d.texture.Texture>`
        for more details and examples.

        EXAMPLES::

            sage: from sage.plot.plot3d.texture import Texture_class
            sage: Texture_class(3, opacity=0.6)
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
            color = (float(color[0]), float(color[1]), float(color[2]))

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
        Gives string representation of the Texture object.

        EXAMPLES::

            sage: from sage.plot.plot3d.texture import Texture
            sage: Texture('yellow')               #indirect doctest
            Texture(texture..., yellow, ffff00)
            sage: Texture((1,1,0), opacity=.5)
            Texture(texture..., ffff00)
        """
        if self.name is not None:
            return "Texture(%s, %s, %s)" % (self.id, self.name, self.hex_rgb())
        else:
            return "Texture(%s, %s)" % (self.id, self.hex_rgb())

    def hex_rgb(self):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.texture import Texture
            sage: Texture('red').hex_rgb()
            'ff0000'
            sage: Texture((1, .5, 0)).hex_rgb()
            'ff7f00'
        """
        return "%02x%02x%02x" % tuple(int(255*s) for s in self.color)

    def tachyon_str(self):
        r"""
        Converts Texture object to string suitable for Tachyon ray tracer.

        EXAMPLES::

            sage: from sage.plot.plot3d.texture import Texture
            sage: t = Texture(opacity=0.6)
            sage: t.tachyon_str()
            'Texdef texture...\n  Ambient 0.333333333333 Diffuse 0.666666666667 Specular 0.0 Opacity 0.6\n   Color 0.4 0.4 1.0\n   TexFunc 0'
        """
        total_color = float(sum(self.ambient) + sum(self.diffuse) + sum(self.specular))
        if total_color == 0:
            total_color = 1
        return "Texdef %s\n" % self.id + \
         "  Ambient %s Diffuse %s Specular %s Opacity %s\n" % \
                (sum(self.ambient)/total_color,
                 sum(self.diffuse)/total_color,
                 sum(self.specular)/total_color,
                 self.opacity) + \
        "   Color %s %s %s\n" % (self.color[0], self.color[1], self.color[2]) + \
        "   TexFunc 0"

    def x3d_str(self):
        r"""
        Converts Texture object to string suitable for x3d.

        EXAMPLES::

            sage: from sage.plot.plot3d.texture import Texture
            sage: t = Texture(opacity=0.6)
            sage: t.x3d_str()
            "<Appearance><Material diffuseColor='0.4 0.4 1.0' shininess='1.0' specularColor='0.0 0.0 0.0'/></Appearance>"
        """
        return "<Appearance><Material diffuseColor='%s %s %s' shininess='%s' specularColor='%s %s %s'/></Appearance>" % \
                (self.color[0], self.color[1], self.color[2], self.shininess, self.specular[0], self.specular[0], self.specular[0])

    def mtl_str(self):
        r"""
        Converts Texture object to string suitable for mtl output.

        EXAMPLES::

            sage: from sage.plot.plot3d.texture import Texture
            sage: t = Texture(opacity=0.6)
            sage: t.mtl_str()
            'newmtl texture...\nKa 0.2 0.2 0.5\nKd 0.4 0.4 1.0\nKs 0.0 0.0 0.0\nillum 1\nNs 1.0\nd 0.6'
        """
        return "\n".join(["newmtl %s" % self.id,
                   "Ka %s %s %s" % self.ambient,
                   "Kd %s %s %s" % self.diffuse,
                   "Ks %s %s %s" % self.specular,
                   "illum %s" % (2 if sum(self.specular) > 0 else 1),
                   "Ns %s" % self.shininess,
                   "d %s" % self.opacity, ])

    def jmol_str(self, obj):
        r"""
        Converts Texture object to string suitable for Jmol applet.

        INPUT:

        - ``obj`` - str

        EXAMPLES::

            sage: from sage.plot.plot3d.texture import Texture
            sage: t = Texture(opacity=0.6)
            sage: t.jmol_str('obj')
            'color obj translucent 0.4 [102,102,255]'

        ::

            sage: sum([dodecahedron(center=[2.5*x, 0, 0], color=(1, 0, 0, x/10)) for x in range(11)]).show(aspect_ratio=[1,1,1], frame=False, zoom=2)
        """
        translucent = "translucent %s" % float(1-self.opacity) if self.opacity < 1 else ""
        return "color %s %s [%s,%s,%s]" % (obj, translucent,
                int(255*self.color[0]), int(255*self.color[1]), int(255*self.color[2]))

