# sage_setup: distribution = sagemath-plot
# sage.doctest: needs sage.plot
r"""
Colors

This module defines a :class:`Color` object and helper functions (see,
e.g., :func:`hue`, :func:`rainbow`), as well as a set of
:data:`colors` and :data:`colormaps` to use with
:class:`Graphics` objects in Sage.

For a list of pre-defined colors in Sage, evaluate::

    sage: sorted(colors)
    ['aliceblue', 'antiquewhite', 'aqua', 'aquamarine', 'automatic', ...]

Apart from 'automatic' which just an alias for 'lightblue', this list
comprises the "official" W3C CSS3_ / SVG_ colors.

.. _CSS3: http://www.w3.org/TR/css3-color/#svg-color
.. _SVG: http://www.w3.org/TR/SVG/types.html#ColorKeywords

For a list of color maps in Sage, evaluate::

    sage: sorted(colormaps)
    ['Accent', ...]

These are imported from matplotlib's colormaps_ collection.

.. _colormaps: https://matplotlib.org/stable/gallery/color/colormap_reference.html
"""

# ****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import math
from collections.abc import MutableMapping
from colorsys import hsv_to_rgb, hls_to_rgb, rgb_to_hsv, rgb_to_hls


colors_dict = {
    'automatic'              : '#add8e6',     # 173, 216, 230
    'aliceblue'              : '#f0f8ff',     # 240, 248, 255
    'antiquewhite'           : '#faebd7',     # 250, 235, 215
    'aqua'                   : '#00ffff',     #   0, 255, 255
    'aquamarine'             : '#7fffd4',     # 127, 255, 212
    'azure'                  : '#f0ffff',     # 240, 255, 255
    'beige'                  : '#f5f5dc',     # 245, 245, 220
    'bisque'                 : '#ffe4c4',     # 255, 228, 196
    'black'                  : '#000000',     #   0,   0,   0
    'blanchedalmond'         : '#ffebcd',     # 255, 235, 205
    'blue'                   : '#0000ff',     #   0,   0, 255
    'blueviolet'             : '#8a2be2',     # 138,  43, 226
    'brown'                  : '#a52a2a',     # 165,  42,  42
    'burlywood'              : '#deb887',     # 222, 184, 135
    'cadetblue'              : '#5f9ea0',     #  95, 158, 160
    'chartreuse'             : '#7fff00',     # 127, 255,   0
    'chocolate'              : '#d2691e',     # 210, 105,  30
    'coral'                  : '#ff7f50',     # 255, 127,  80
    'cornflowerblue'         : '#6495ed',     # 100, 149, 237
    'cornsilk'               : '#fff8dc',     # 255, 248, 220
    'crimson'                : '#dc143c',     # 220,  20,  60
    'cyan'                   : '#00ffff',     #   0, 255, 255
    'darkblue'               : '#00008b',     #   0,   0, 139
    'darkcyan'               : '#008b8b',     #   0, 139, 139
    'darkgoldenrod'          : '#b8860b',     # 184, 134,  11
    'darkgray'               : '#a9a9a9',     # 169, 169, 169
    'darkgreen'              : '#006400',     #   0, 100,   0
    'darkgrey'               : '#a9a9a9',     # 169, 169, 169
    'darkkhaki'              : '#bdb76b',     # 189, 183, 107
    'darkmagenta'            : '#8b008b',     # 139,   0, 139
    'darkolivegreen'         : '#556b2f',     #  85, 107,  47
    'darkorange'             : '#ff8c00',     # 255, 140,   0
    'darkorchid'             : '#9932cc',     # 153,  50, 204
    'darkred'                : '#8b0000',     # 139,   0,   0
    'darksalmon'             : '#e9967a',     # 233, 150, 122
    'darkseagreen'           : '#8fbc8f',     # 143, 188, 143
    'darkslateblue'          : '#483d8b',     #  72,  61, 139
    'darkslategray'          : '#2f4f4f',     #  47,  79,  79
    'darkslategrey'          : '#2f4f4f',     #  47,  79,  79
    'darkturquoise'          : '#00ced1',     #   0, 206, 209
    'darkviolet'             : '#9400d3',     # 148,   0, 211
    'deeppink'               : '#ff1493',     # 255,  20, 147
    'deepskyblue'            : '#00bfff',     #   0, 191, 255
    'dimgray'                : '#696969',     # 105, 105, 105
    'dimgrey'                : '#696969',     # 105, 105, 105
    'dodgerblue'             : '#1e90ff',     #  30, 144, 255
    'firebrick'              : '#b22222',     # 178,  34,  34
    'floralwhite'            : '#fffaf0',     # 255, 250, 240
    'forestgreen'            : '#228b22',     #  34, 139,  34
    'fuchsia'                : '#ff00ff',     # 255,   0, 255
    'gainsboro'              : '#dcdcdc',     # 220, 220, 220
    'ghostwhite'             : '#f8f8ff',     # 248, 248, 255
    'gold'                   : '#ffd700',     # 255, 215,   0
    'goldenrod'              : '#daa520',     # 218, 165,  32
    'gray'                   : '#808080',     # 128, 128, 128
    'green'                  : '#008000',     #   0, 128,   0
    'greenyellow'            : '#adff2f',     # 173, 255,  47
    'grey'                   : '#808080',     # 128, 128, 128
    'honeydew'               : '#f0fff0',     # 240, 255, 240
    'hotpink'                : '#ff69b4',     # 255, 105, 180
    'indianred'              : '#cd5c5c',     # 205,  92,  92
    'indigo'                 : '#4b0082',     #  75,   0, 130
    'ivory'                  : '#fffff0',     # 255, 255, 240
    'khaki'                  : '#f0e68c',     # 240, 230, 140
    'lavender'               : '#e6e6fa',     # 230, 230, 250
    'lavenderblush'          : '#fff0f5',     # 255, 240, 245
    'lawngreen'              : '#7cfc00',     # 124, 252,   0
    'lemonchiffon'           : '#fffacd',     # 255, 250, 205
    'lightblue'              : '#add8e6',     # 173, 216, 230
    'lightcoral'             : '#f08080',     # 240, 128, 128
    'lightcyan'              : '#e0ffff',     # 224, 255, 255
    'lightgoldenrodyellow'   : '#fafad2',     # 250, 250, 210
    'lightgray'              : '#d3d3d3',     # 211, 211, 211
    'lightgreen'             : '#90ee90',     # 144, 238, 144
    'lightgrey'              : '#d3d3d3',     # 211, 211, 211
    'lightpink'              : '#ffb6c1',     # 255, 182, 193
    'lightsalmon'            : '#ffa07a',     # 255, 160, 122
    'lightseagreen'          : '#20b2aa',     #  32, 178, 170
    'lightskyblue'           : '#87cefa',     # 135, 206, 250
    'lightslategray'         : '#778899',     # 119, 136, 153
    'lightslategrey'         : '#778899',     # 119, 136, 153
    'lightsteelblue'         : '#b0c4de',     # 176, 196, 222
    'lightyellow'            : '#ffffe0',     # 255, 255, 224
    'lime'                   : '#00ff00',     #   0, 255,   0
    'limegreen'              : '#32cd32',     #  50, 205,  50
    'linen'                  : '#faf0e6',     # 250, 240, 230
    'magenta'                : '#ff00ff',     # 255,   0, 255
    'maroon'                 : '#800000',     # 128,   0,   0
    'mediumaquamarine'       : '#66cdaa',     # 102, 205, 170
    'mediumblue'             : '#0000cd',     #   0,   0, 205
    'mediumorchid'           : '#ba55d3',     # 186,  85, 211
    'mediumpurple'           : '#9370db',     # 147, 112, 219
    'mediumseagreen'         : '#3cb371',     #  60, 179, 113
    'mediumslateblue'        : '#7b68ee',     # 123, 104, 238
    'mediumspringgreen'      : '#00fa9a',     #   0, 250, 154
    'mediumturquoise'        : '#48d1cc',     #  72, 209, 204
    'mediumvioletred'        : '#c71585',     # 199,  21, 133
    'midnightblue'           : '#191970',     #  25,  25, 112
    'mintcream'              : '#f5fffa',     # 245, 255, 250
    'mistyrose'              : '#ffe4e1',     # 255, 228, 225
    'moccasin'               : '#ffe4b5',     # 255, 228, 181
    'navajowhite'            : '#ffdead',     # 255, 222, 173
    'navy'                   : '#000080',     #   0,   0, 128
    'oldlace'                : '#fdf5e6',     # 253, 245, 230
    'olive'                  : '#808000',     # 128, 128,   0
    'olivedrab'              : '#6b8e23',     # 107, 142,  35
    'orange'                 : '#ffa500',     # 255, 165,   0
    'orangered'              : '#ff4500',     # 255,  69,   0
    'orchid'                 : '#da70d6',     # 218, 112, 214
    'palegoldenrod'          : '#eee8aa',     # 238, 232, 170
    'palegreen'              : '#98fb98',     # 152, 251, 152
    'paleturquoise'          : '#afeeee',     # 175, 238, 238
    'palevioletred'          : '#db7093',     # 219, 112, 147
    'papayawhip'             : '#ffefd5',     # 255, 239, 213
    'peachpuff'              : '#ffdab9',     # 255, 218, 185
    'peru'                   : '#cd853f',     # 205, 133,  63
    'pink'                   : '#ffc0cb',     # 255, 192, 203
    'plum'                   : '#dda0dd',     # 221, 160, 221
    'powderblue'             : '#b0e0e6',     # 176, 224, 230
    'purple'                 : '#800080',     # 128,   0, 128
    'red'                    : '#ff0000',     # 255,   0,   0
    'rosybrown'              : '#bc8f8f',     # 188, 143, 143
    'royalblue'              : '#4169e1',     #  65, 105, 225
    'saddlebrown'            : '#8b4513',     # 139,  69,  19
    'salmon'                 : '#fa8072',     # 250, 128, 114
    'sandybrown'             : '#f4a460',     # 244, 164,  96
    'seagreen'               : '#2e8b57',     #  46, 139,  87
    'seashell'               : '#fff5ee',     # 255, 245, 238
    'sienna'                 : '#a0522d',     # 160,  82,  45
    'silver'                 : '#c0c0c0',     # 192, 192, 192
    'skyblue'                : '#87ceeb',     # 135, 206, 235
    'slateblue'              : '#6a5acd',     # 106,  90, 205
    'slategray'              : '#708090',     # 112, 128, 144
    'slategrey'              : '#708090',     # 112, 128, 144
    'snow'                   : '#fffafa',     # 255, 250, 250
    'springgreen'            : '#00ff7f',     #   0, 255, 127
    'steelblue'              : '#4682b4',     #  70, 130, 180
    'tan'                    : '#d2b48c',     # 210, 180, 140
    'teal'                   : '#008080',     #   0, 128, 128
    'thistle'                : '#d8bfd8',     # 216, 191, 216
    'tomato'                 : '#ff6347',     # 255,  99,  71
    'turquoise'              : '#40e0d0',     #  64, 224, 208
    'violet'                 : '#ee82ee',     # 238, 130, 238
    'wheat'                  : '#f5deb3',     # 245, 222, 179
    'white'                  : '#ffffff',     # 255, 255, 255
    'whitesmoke'             : '#f5f5f5',     # 245, 245, 245
    'yellow'                 : '#ffff00',     # 255, 255,   0
    'yellowgreen'            : '#9acd32'      # 154, 205,  50
}


def mod_one(x):
    """
    Reduce a number modulo 1.

    INPUT:

    - ``x`` -- an instance of Integer, int, RealNumber, etc.; the
      number to reduce

    OUTPUT: float

    EXAMPLES::

        sage: from sage.plot.colors import mod_one
        sage: mod_one(1)
        1.0
        sage: mod_one(7.0)
        0.0
        sage: mod_one(-11/7)
        0.4285714285714286
        sage: mod_one(pi) + mod_one(-pi)                                                # needs sage.symbolic
        1.0
    """
    x = float(x)
    if x != 1:
        x = math.modf(x)[0]
        if x < 0:
            x += 1
    return x


def html_to_float(c):
    """
    Convert a HTML hex color to a Red-Green-Blue (RGB) tuple.

    INPUT:

    - ``c`` -- string; a valid HTML hex color

    OUTPUT: a RGB 3-tuple of floats in the interval [0.0, 1.0]

    EXAMPLES::

        sage: from sage.plot.colors import html_to_float
        sage: html_to_float('#fff')
        (1.0, 1.0, 1.0)
        sage: html_to_float('#abcdef')
        (0.6705882352941176, 0.803921568627451, 0.9372549019607843)
        sage: html_to_float('#123xyz')
        Traceback (most recent call last):
        ...
        ValueError: invalid literal for int() with base 16: '3x'
    """
    if not c or c[0] != '#':
        raise ValueError("'%s' must be a valid HTML hex color (e.g., '#f07' or '#d6e7da')" % c)
    h = c[1:]
    if len(h) == 3:
        h = f'{h[0]}{h[0]}{h[1]}{h[1]}{h[2]}{h[2]}'
    elif len(h) != 6:
        raise ValueError("color hex string (= '%s') must have length 3 or 6" % h)
    return tuple([int(h[i:i + 2], base=16) / 255 for i in [0, 2, 4]])


def rgbcolor(c, space='rgb'):
    """
    Convert a color (string, tuple, list, or :class:`Color`) to a
    mod-one reduced (see :func:`mod_one`) valid Red-Green-Blue (RGB)
    tuple.  The returned tuple is also a valid matplotlib RGB color.

    INPUT:

    - ``c`` -- a :class:`Color` instance, string (name or HTML hex),
      3-tuple, or 3-list; the color to convert

    - ``space`` -- string (default: ``'rgb'``); the color space
      coordinate system (other choices are ``'hsl'``, ``'hls'``, and ``'hsv'``)
      in which to interpret a 3-tuple or 3-list

    OUTPUT: a RGB 3-tuple of floats in the interval [0.0, 1.0]

    EXAMPLES::

        sage: from sage.plot.colors import rgbcolor
        sage: rgbcolor(Color(0.25, 0.4, 0.9))
        (0.25, 0.4, 0.9)
        sage: rgbcolor('purple')
        (0.5019607843137255, 0.0, 0.5019607843137255)
        sage: rgbcolor('#fa0')
        (1.0, 0.6666666666666666, 0.0)
        sage: rgbcolor('#ffffff')
        (1.0, 1.0, 1.0)
        sage: rgbcolor((1,1/2,1/3))
        (1.0, 0.5, 0.3333333333333333)
        sage: rgbcolor([1,1/2,1/3])
        (1.0, 0.5, 0.3333333333333333)
        sage: rgbcolor((1,1,1), space='hsv')
        (1.0, 0.0, 0.0)
        sage: rgbcolor((0.5,0.75,1), space='hls')
        (0.5, 0.9999999999999999, 1.0)
        sage: rgbcolor((0.5,1.0,0.75), space='hsl')
        (0.5, 0.9999999999999999, 1.0)
        sage: rgbcolor([1,2,255])   # WARNING -- numbers are reduced mod 1!!
        (1.0, 0.0, 0.0)
        sage: rgbcolor('#abcd')
        Traceback (most recent call last):
        ...
        ValueError: color hex string (= 'abcd') must have length 3 or 6
        sage: rgbcolor('fff')
        Traceback (most recent call last):
        ...
        ValueError: unknown color 'fff'
        sage: rgbcolor(1)
        Traceback (most recent call last):
        ...
        TypeError: '1' must be a Color, list, tuple, or string
        sage: rgbcolor((0.2,0.8,1), space='grassmann')
        Traceback (most recent call last):
        ...
        ValueError: space must be one of 'rgb', 'hsv', 'hsl', 'hls'
        sage: rgbcolor([0.4, 0.1])
        Traceback (most recent call last):
        ...
        ValueError: color list or tuple '[0.400000000000000, 0.100000000000000]' must have 3 entries, one for each RGB, HSV, HLS, or HSL channel
    """
    if isinstance(c, Color):
        return c.rgb()

    if isinstance(c, str):
        if len(c) > 0 and c[0] == '#':
            # Assume an HTML-like color, e.g., #00ffff or #ab0.
            return html_to_float(c)
        else:
            try:
                return colors[c].rgb()
            except KeyError:
                raise ValueError("unknown color '%s'" % c)

    elif isinstance(c, (list, tuple)):
        if len(c) != 3:
            raise ValueError(f"color list or tuple '{c}' must have 3 entries, one for each RGB, HSV, HLS, or HSL channel")
        c = [mod_one(comp) for comp in c]
        if space == 'rgb':
            return tuple(c)
        elif space == 'hsv':
            return tuple(map(float, hsv_to_rgb(*c)))
        elif space == 'hls':
            return tuple(map(float, hls_to_rgb(*c)))
        elif space == 'hsl':
            return tuple(map(float, hls_to_rgb(c[0], c[2], c[1])))
        else:
            raise ValueError("space must be one of 'rgb', 'hsv', 'hsl', 'hls'")

    raise TypeError("'%s' must be a Color, list, tuple, or string" % c)


# For backward compatibility.
to_mpl_color = rgbcolor


class Color:
    def __init__(self, r='#0000ff', g=None, b=None, space='rgb'):
        """
        A Red-Green-Blue (RGB) color model color object.  For most
        consumer-grade devices (e.g., CRTs, LCDs, and printers), as
        well as internet applications, this is a point in the sRGB
        absolute color space.  The Hue-Saturation-Lightness (HSL),
        Hue-Lightness-Saturation (HLS), and Hue-Saturation-Value (HSV)
        spaces are useful alternate representations, or coordinate
        transformations, of this space.  Coordinates in all of these
        spaces are floating point values in the interval [0.0, 1.0].

        .. NOTE:: All instantiations of :class:`Color` are converted
                  to an internal RGB floating point 3-tuple.  This is
                  likely to degrade precision.

        INPUT:

        - ``r``, ``g``, ``b`` -- either a triple of floats between 0 and 1,
          OR ``r`` - a color name string or HTML color hex string

        - ``space`` -- string (default: ``'rgb'``); the coordinate system
          (other choices are ``'hsl'``, ``'hls'``, and ``'hsv'``) in which to
          interpret a triple of floats

        EXAMPLES::

            sage: Color('purple')
            RGB color (0.5019607843137255, 0.0, 0.5019607843137255)
            sage: Color('#8000ff')
            RGB color (0.5019607843137255, 0.0, 1.0)
            sage: Color(0.5,0,1)
            RGB color (0.5, 0.0, 1.0)
            sage: Color(0.5, 1.0, 1, space='hsv')
            RGB color (0.0, 1.0, 1.0)
            sage: Color(0.25, 0.5, 0.5, space='hls')
            RGB color (0.5000000000000001, 0.75, 0.25)
            sage: Color(1, 0, 1/3, space='hsl')
            RGB color (0.3333333333333333, 0.3333333333333333, 0.3333333333333333)
            sage: from sage.plot.colors import chocolate
            sage: Color(chocolate)
            RGB color (0.8235294117647058, 0.4117647058823529, 0.11764705882352941)
        """
        if g is None and b is None:
            self._rgb = rgbcolor(r)
        else:
            self._rgb = rgbcolor((r, g, b), space=space)

    def __repr__(self):
        """
        Return a string representation of this color.

        OUTPUT: string

        EXAMPLES::

            sage: Color('#8000ff').__repr__()
            'RGB color (0.5019607843137255, 0.0, 1.0)'
            sage: Color(1, 0.5, 1/16, space='hsl').__repr__()
            'RGB color (0.09375, 0.03125, 0.03125)'
        """
        return f"RGB color {self._rgb}"

    def __lt__(self, right):
        """
        Check whether a :class:`Color` object is less than some other
        object. This doesn't make sense, and so we conclude that it is
        not less than the other object.

        INPUT:

        - ``right`` -- an object

        OUTPUT: boolean; ``False``

        EXAMPLES::

            sage: Color('red') < Color('red')
            False
            sage: Color('blue') < Color('red')
            False
            sage: Color('red') < "xyzzy"
            False
        """
        return False

    def __le__(self, right):
        """
        Check whether a :class:`Color` object is less than or equal to
        some other object. It wouldn't make sense for it to be less than
        the other object, so we treat this the same as an equality
        check.

        INPUT:

        - ``right`` -- an object

        OUTPUT: boolean; ``False``

        EXAMPLES::

            sage: Color('red') <= Color('red')
            True
            sage: Color('blue') <= Color('red')
            False
            sage: Color('red') <= "xyzzy"
            False
        """
        return self == right

    def __eq__(self, right):
        """
        Compare two :class:`Color` objects to determine whether
        they refer to the same color.

        INPUT:

        - ``right`` -- a :class:`Color` instance

        OUTPUT: boolean; ``True`` if the two colors are the same, ``False``
        if different

        EXAMPLES::

            sage: Color('red') == Color((1,0,0))
            True
            sage: Color('blue') == Color((0,1,0))
            False
            sage: Color('blue') + Color((0,1,0)) == Color((0,0.5,0.5))
            True
            sage: Color(0.2,0.3,0.2) == False
            False
        """
        if isinstance(right, Color):
            return self._rgb == right._rgb
        return False

    def __ne__(self, right):
        """
        Compare two :class:`Color` objects to determine whether
        they refer to different colors.

        INPUT:

        - ``right`` -- a :class:`Color` instance

        OUTPUT:

        boolean; ``True`` if the two colors are different, ``False`` if they're
        the same.

        EXAMPLES::

            sage: Color('green') != Color('yellow')
            True
            sage: Color('red') != Color(1,0,0)
            False
            sage: Color('yellow') != Color(1,1,0)
            False
            sage: Color('blue') != 23
            True
        """
        return not (self == right)

    def __gt__(self, right):
        """
        Check whether a :class:`Color` object is greater than some other
        object. This doesn't make sense, and so we conclude that it is
        not greater than the other object.

        INPUT:

        - ``right`` -- an object

        OUTPUT: boolean; ``False``

        EXAMPLES::

            sage: Color('red') > Color('red')
            False
            sage: Color('blue') > Color('red')
            False
            sage: Color('red') > "xyzzy"
            False
        """
        return False

    def __ge__(self, right):
        """
        Check whether a :class:`Color` object is greater than or equal
        to some other object. It wouldn't make sense for it to be
        greater than the other object, so we treat this the same as an
        equality check.

        INPUT:

        - ``right`` -- an object

        OUTPUT: boolean; ``False``

        EXAMPLES::

            sage: Color('red') >= Color('red')
            True
            sage: Color('blue') >= Color('red')
            False
            sage: Color('red') >= "xyzzy"
            False
        """
        return self == right

    def __hash__(self):
        """
        Return the hash value of a color.
        Equal colors return equal hash values.

        OUTPUT: a hash value

        EXAMPLES::

            sage: hash(Color('red')) # random
            873150856
            sage: hash(Color('red')) == hash(Color((1,0,0)))
            True
        """
        return hash(self._rgb)

    def blend(self, color, fraction=0.5):
        """
        Return a color blended with the given ``color`` by a given
        ``fraction``.  The algorithm interpolates linearly between the
        colors' corresponding R, G, and B coordinates.

        INPUT:

        - ``color`` -- a :class:`Color` instance or float-convertible
          3-tuple/list; the color with which to blend this color

        - ``fraction`` -- a float-convertible number; the fraction of
          ``color`` to blend with this color

        OUTPUT: a **new** :class:`Color` instance

        EXAMPLES::

            sage: from sage.plot.colors import red, blue, lime
            sage: red.blend(blue)
            RGB color (0.5, 0.0, 0.5)
            sage: red.blend(blue, fraction=0.0)
            RGB color (1.0, 0.0, 0.0)
            sage: red.blend(blue, fraction=1.0)
            RGB color (0.0, 0.0, 1.0)
            sage: lime.blend((0.3, 0.5, 0.7))
            RGB color (0.15, 0.75, 0.35)
            sage: blue.blend(blue)
            RGB color (0.0, 0.0, 1.0)
            sage: red.blend(lime, fraction=0.3)
            RGB color (0.7, 0.3, 0.0)
            sage: blue.blend((0.0, 0.9, 0.2), fraction=0.2)
            RGB color (0.0, 0.18000000000000002, 0.8400000000000001)
            sage: red.blend(0.2)
            Traceback (most recent call last):
            ...
            TypeError: 0.200000000000000 must be a Color or float-convertible 3-tuple/list
        """
        fraction = float(fraction)
        if isinstance(color, Color):
            color = color._rgb
        if isinstance(color, (list, tuple)) and len(color) == 3:
            color = [float(_) for _ in color]
            return Color(rgbcolor([(1 - fraction) * a + fraction * b
                                   for a, b in zip(self._rgb, color)]))
        raise TypeError(f"{color} must be a Color or float-convertible 3-tuple/list")

    def __add__(self, right):
        """
        Return a color "added" on the right to another color, with
        :meth:`blend`.

        INPUT:

        - ``right`` -- a :class:`Color` instance or float-convertible
          3-tuple/list

        OUTPUT: a **new** :class:`Color` instance

        EXAMPLES::

            sage: from sage.plot.colors import red, blue, lime
            sage: red + blue + lime
            RGB color (0.25, 0.5, 0.25)
            sage: from sage.plot.colors import cyan, magenta, yellow
            sage: cyan + magenta + yellow
            RGB color (0.75, 0.75, 0.5)
            sage: c1 = Color(0.1, 0.5, 0.8); c2 = Color(0.2, 0.4, 0.7, space='hsv')
            sage: c1 + 0.1
            Traceback (most recent call last):
            ...
            TypeError: 0.100000000000000 must be a Color or float-convertible 3-tuple/list
            sage: c2 + [0.5, 0.2, 0.9]
            RGB color (0.572, 0.44999999999999996, 0.66)
            sage: c1.__add__(red).__add__((0.9, 0.2, 1/3))
            RGB color (0.7250000000000001, 0.225, 0.3666666666666667)
            sage: c1 + c2
            RGB color (0.37199999999999994, 0.6, 0.61)
        """
        return self.blend(right)

    def __radd__(self, left):
        """
        Return a color "added" on the left to another color, with
        :meth:`blend`.

        INPUT:

        - ``left`` -- a :class:`Color` instance or float-convertible
          3-tuple/list

        OUTPUT: a **new** :class:`Color` instance

        EXAMPLES::

            sage: from sage.plot.colors import olive, orchid
            sage: olive + orchid
            RGB color (0.6784313725490196, 0.47058823529411764, 0.4196078431372549)
            sage: d1 = Color(0.1, 0.5, 0.8, space='hls'); d2 = Color(0.2, 0.4, 0.7)
            sage: [0.5, 0.2, 0.9] + d2
            RGB color (0.35, 0.30000000000000004, 0.8)
            sage: 0.1 + d1
            Traceback (most recent call last):
            ...
            TypeError: 0.100000000000000 must be a Color or float-convertible 3-tuple/list
            sage: d2.__radd__(Color('brown')).__radd__((0.9, 0.2, 1/3))
            RGB color (0.661764705882353, 0.2411764705882353, 0.38284313725490193)
        """
        return self + left

    def __mul__(self, right):
        """
        Return a color whose RGB coordinates are this color's
        coordinates multiplied on the right by a scalar.

        INPUT:

        - ``right`` -- a float-convertible number

        OUTPUT: a **new** :class:`Color` instance

        EXAMPLES::

            sage: Color('yellow') * 0.5
            RGB color (0.5, 0.5, 0.0)
            sage: Color('yellow') * (9.0 / 8.0)   # reduced modulo 1.0
            RGB color (0.125, 0.125, 0.0)
            sage: from sage.plot.colors import cyan, grey, indianred
            sage: cyan * 0.3 + grey * 0.1 + indianred * 0.6
            RGB color (0.25372549019607843, 0.1957843137254902, 0.1957843137254902)
            sage: indianred.__mul__(42)
            RGB color (0.764705882352942, 0.1529411764705877, 0.1529411764705877)
        """
        right = float(right)
        return Color([x * right for x in self._rgb])

    def __rmul__(self, left):
        """
        Return a color whose RGB coordinates are this color's
        coordinates multiplied on the left by a scalar.

        INPUT:

        - ``left`` -- a float-convertible number

        OUTPUT: a **new** :class:`Color` instance

        EXAMPLES::

            sage: from sage.plot.colors import aqua, cornsilk, tomato
            sage: 0.3 * aqua
            RGB color (0.0, 0.3, 0.3)
            sage: Color('indianred').__rmul__(42)
            RGB color (0.764705882352942, 0.1529411764705877, 0.1529411764705877)
        """
        return self * left

    def __truediv__(self, right):
        """
        Return a color whose RGB coordinates are this color's
        coordinates divided by a scalar.

        INPUT:

        - ``right`` -- a float-convertible, nonzero number

        OUTPUT: a **new** instance of :class:`Color`

        EXAMPLES::

            sage: from sage.plot.colors import papayawhip, yellow
            sage: yellow / 4
            RGB color (0.25, 0.25, 0.0)
            sage: yellow.__truediv__(4)
            RGB color (0.25, 0.25, 0.0)
            sage: (papayawhip + Color(0.5, 0.5, 0.1) + yellow) / 3.0
            RGB color (0.29166666666666663, 0.286437908496732, 0.07794117647058824)
            sage: vector((papayawhip / 2).rgb()) == vector((papayawhip * 0.5).rgb())
            True
            sage: yellow.__truediv__(1/4)
            RGB color (0.0, 0.0, 0.0)

        TESTS::

            sage: Color('black') / 0.0
            Traceback (most recent call last):
            ...
            ZeroDivisionError: float division by zero

            sage: papayawhip / yellow
            Traceback (most recent call last):
            ...
            TypeError: float() argument must be a string or a... number...
        """
        return self * (1 / float(right))

    def __int__(self):
        """
        Return the integer representation of this colour.

        OUTPUT:

        The integer 256^2 r_int + 256 g_int + b_int, where r_int, g_int,
        and b_int are obtained from r, g, and b by converting from the
        real interval [0.0, 1.0] to the integer range 0, 1, ..., 255.

        EXAMPLES::

            sage: from sage.plot.colors import whitesmoke
            sage: int(whitesmoke)
            16119285
        """
        return float_to_integer(*self._rgb)

    def __iter__(self):
        """
        Return an iterator over the RGB coordinates of this color.

        OUTPUT: a tupleiterator

        EXAMPLES::

            sage: from sage.plot.colors import dodgerblue, maroon
            sage: r, g, b = dodgerblue
            sage: r
            0.11764705882352941
            sage: g
            0.5647058823529412
            sage: b
            1.0
            sage: vector(maroon) == vector(Color(maroon)) == vector(Color('maroon'))
            True
        """
        return iter(self._rgb)

    def __getitem__(self, i):
        """
        Return the Red (0th), Green (1st), or Blue (2nd) coordinate of this
        color via index access.

        INPUT:

        - ``i`` -- integer; the 0-based coordinate to retrieve

        OUTPUT: float

        EXAMPLES::

            sage: from sage.plot.colors import crimson, midnightblue
            sage: Color('#badfad')[0]
            0.7294117647058823
            sage: (crimson[0], crimson[1], crimson[2]) == crimson.rgb()
            True
            sage: midnightblue[2] == midnightblue[-1]
            True
            sage: midnightblue[3]
            Traceback (most recent call last):
            ...
            IndexError: tuple index out of range
        """
        return self._rgb[i]

    def rgb(self):
        """
        Return the underlying Red-Green-Blue (RGB) coordinates of this
        color.

        OUTPUT: a 3-tuple of floats

        EXAMPLES::

            sage: Color(0.3, 0.5, 0.7).rgb()
            (0.3, 0.5, 0.7)
            sage: Color('#8000ff').rgb()
            (0.5019607843137255, 0.0, 1.0)
            sage: from sage.plot.colors import orange
            sage: orange.rgb()
            (1.0, 0.6470588235294118, 0.0)
            sage: Color('magenta').rgb()
            (1.0, 0.0, 1.0)
            sage: Color(1, 0.7, 0.9, space='hsv').rgb()
            (0.9, 0.2700000000000001, 0.2700000000000001)
        """
        return self._rgb

    def hls(self):
        """
        Return the Hue-Lightness-Saturation (HLS) coordinates of this
        color.

        OUTPUT: a 3-tuple of floats

        EXAMPLES::

            sage: Color(0.3, 0.5, 0.7, space='hls').hls()
            (0.30000000000000004, 0.5, 0.7)
            sage: Color(0.3, 0.5, 0.7, space='hsl').hls() # abs tol 1e-15
            (0.30000000000000004, 0.7, 0.5000000000000001)
            sage: Color('#aabbcc').hls() # abs tol 1e-15
            (0.5833333333333334, 0.7333333333333334, 0.25000000000000017)
            sage: from sage.plot.colors import orchid
            sage: orchid.hls() # abs tol 1e-15
            (0.8396226415094339, 0.6470588235294117, 0.5888888888888889)
        """
        return tuple(map(float, rgb_to_hls(*self._rgb)))

    def hsl(self):
        """
        Return the Hue-Saturation-Lightness (HSL) coordinates of this
        color.

        OUTPUT: a 3-tuple of floats

        EXAMPLES::

            sage: Color(1,0,0).hsl()
            (0.0, 1.0, 0.5)
            sage: from sage.plot.colors import orchid
            sage: orchid.hsl() # abs tol 1e-15
            (0.8396226415094339, 0.5888888888888889, 0.6470588235294117)
            sage: Color('#aabbcc').hsl() # abs tol 1e-15
            (0.5833333333333334, 0.25000000000000017, 0.7333333333333334)
        """
        h, l, s = tuple(map(float, rgb_to_hls(*self._rgb)))
        return (h, s, l)

    def hsv(self):
        """
        Return the Hue-Saturation-Value (HSV) coordinates of this
        color.

        OUTPUT: a 3-tuple of floats

        EXAMPLES::

            sage: from sage.plot.colors import red
            sage: red.hsv()
            (0.0, 1.0, 1.0)
            sage: Color(1,1,1).hsv()
            (0.0, 0.0, 1.0)
            sage: Color('gray').hsv()
            (0.0, 0.0, 0.5019607843137255)
        """
        return tuple(map(float, rgb_to_hsv(*self._rgb)))

    def html_color(self):
        """
        Return a HTML hex representation for this color.

        OUTPUT: string of length 7

        EXAMPLES::

            sage: Color('yellow').html_color()
            '#ffff00'
            sage: Color('#fedcba').html_color()
            '#fedcba'
            sage: Color(0.0, 1.0, 0.0).html_color()
            '#00ff00'
            sage: from sage.plot.colors import honeydew
            sage: honeydew.html_color()
            '#f0fff0'
        """
        return float_to_html(*self._rgb)

    def lighter(self, fraction=1/3):
        """
        Return a lighter "shade" of this RGB color by
        :meth:`blend`-ing it with white.  This is **not** an inverse
        of :meth:`darker`.

        INPUT:

        - ``fraction`` -- a float (default: 1/3); blending fraction
          to apply

        OUTPUT: a **new** instance of :class:`Color`

        EXAMPLES::

            sage: from sage.plot.colors import khaki
            sage: khaki.lighter()
            RGB color (0.9607843137254903, 0.934640522875817, 0.6993464052287582)
            sage: Color('white').lighter().darker()
            RGB color (0.6666666666666667, 0.6666666666666667, 0.6666666666666667)
            sage: Color('#abcdef').lighter(1/4)
            RGB color (0.7529411764705882, 0.8529411764705883, 0.9529411764705882)
            sage: Color(1, 0, 8/9, space='hsv').lighter()
            RGB color (0.925925925925926, 0.925925925925926, 0.925925925925926)
        """
        return self.blend((1.0, 1.0, 1.0), fraction)

    def darker(self, fraction=1/3):
        """
        Return a darker "shade" of this RGB color by :meth:`blend`-ing
        it with black.  This is **not** an inverse of :meth:`lighter`.

        INPUT:

        - ``fraction`` -- a float (default: 1/3); blending fraction
          to apply

        OUTPUT: a new instance of :class:`Color`

        EXAMPLES::

            sage: from sage.plot.colors import black
            sage: vector(black.darker().rgb()) == vector(black.rgb())
            True
            sage: Color(0.4, 0.6, 0.8).darker(0.1)
            RGB color (0.36000000000000004, 0.54, 0.7200000000000001)
            sage: Color(.1,.2,.3,space='hsl').darker()
            RGB color (0.24000000000000002, 0.20800000000000002, 0.16)
        """
        return self.blend((0.0, 0.0, 0.0), fraction)


class ColorsDict(dict):
    """
    A dict-like collection of colors, accessible via key or attribute.
    For a list of color names, evaluate::

        sage: sorted(colors)
        ['aliceblue', 'antiquewhite', 'aqua', 'aquamarine', ...]
    """
    def __init__(self):
        """
        Construct a dict-like collection of colors.  The keys are the
        color names (i.e., strings) and the values are RGB 3-tuples of
        floats.

        EXAMPLES::

            sage: from sage.plot.colors import ColorsDict
            sage: cols = ColorsDict()
            sage: set([(type(c), type(cols[c])) for c in cols])
            {(<... 'str'>, <class 'sage.plot.colors.Color'>)}
            sage: sorted(cols)
            ['aliceblue', 'antiquewhite', 'aqua', 'aquamarine', ...]
            sage: len(cols)
            148
        """
        # Convert the colors_dict defined above to Color instances.
        for k in colors_dict:
            self[k] = Color(colors_dict[k])

    def __getattr__(self, name):
        """
        Get a color via attribute access.

        INPUT:

        - ``name`` -- string; the name of the color to return

        OUTPUT: a RGB 3-tuple of floats

        EXAMPLES::

            sage: from sage.plot.colors import ColorsDict, blue
            sage: cols = ColorsDict()
            sage: cols.blue
            RGB color (0.0, 0.0, 1.0)
            sage: cols['blue']
            RGB color (0.0, 0.0, 1.0)
            sage: blue
            RGB color (0.0, 0.0, 1.0)
            sage: cols.punk
            Traceback (most recent call last):
            ...
            AttributeError: 'ColorsDict' has no attribute or colormap punk...
        """
        try:
            return self[name]
        except KeyError:
            raise AttributeError("'{}' has no attribute or colormap {}".format(type(self).__name__, name))

    def __dir__(self):
        """
        Return an approximate list of attribute names, including the
        color names.

        OUTPUT: list of strings

        EXAMPLES::

            sage: from sage.plot.colors import ColorsDict
            sage: cols = ColorsDict()
            sage: 'green' in dir(cols)
            True
        """
        methods = ['__dir__', '__getattr__']
        return dir(super()) + methods + sorted(self)


colors = ColorsDict()


# Add convenient module-scope colors.  These are Color instances.
for c in colors:
    vars()[c] = colors[c]


def hue(h, s=1, v=1):
    r"""
    Convert a Hue-Saturation-Value (HSV) color tuple to a valid
    Red-Green-Blue (RGB) tuple.  All three inputs should lie in the
    interval [0.0, 1.0]; otherwise, they are reduced modulo 1 (see
    :func:`mod_one`).  In particular ``h=0`` and ``h=1`` yield red,
    with the intermediate hues orange, yellow, green, cyan, blue, and
    violet as ``h`` increases.

    This function makes it easy to sample a broad range of colors for
    graphics::

        sage: # needs sage.symbolic
        sage: p = Graphics()
        sage: for phi in xsrange(0, 2 * pi, 1 / pi):
        ....:     p += plot(sin(x + phi), (x, -7, 7), rgbcolor=hue(phi))
        sage: p
        Graphics object consisting of 20 graphics primitives

    INPUT:

    - ``h`` -- a number; the color's hue

    - ``s`` -- a number (default: 1); the color's saturation

    - ``v`` -- a number (default: 1); the color's value

    OUTPUT: a RGB 3-tuple of floats in the interval [0.0, 1.0]

    EXAMPLES::

        sage: hue(0.6)
        (0.0, 0.40000000000000036, 1.0)
        sage: from sage.plot.colors import royalblue
        sage: royalblue
        RGB color (0.2549019607843137, 0.4117647058823529, 0.8823529411764706)
        sage: hue(*royalblue.hsv())
        (0.2549019607843137, 0.4117647058823529, 0.8823529411764706)
        sage: hue(.5, .5, .5)
        (0.25, 0.5, 0.5)

    .. NOTE::

        The HSV to RGB coordinate transformation itself is
        given in the source code for the Python library's
        :mod:`colorsys` module::

            sage: from colorsys import hsv_to_rgb    # not tested
            sage: hsv_to_rgb??                       # not tested
    """
    return tuple(map(float, hsv_to_rgb(mod_one(h), mod_one(s), mod_one(v))))


def float_to_html(r, g, b):
    """
    Convert a Red-Green-Blue (RGB) color tuple to a HTML hex color.

    Each input value should be in the interval [0.0, 1.0]; otherwise,
    the values are first reduced modulo one (see :func:`mod_one`).

    INPUT:

    - ``r`` -- a real number; the RGB color's "red" intensity

    - ``g`` -- a real number; the RGB color's "green" intensity

    - ``b`` -- a real number; the RGB color's "blue" intensity

    OUTPUT: string of length 7, starting with ``'#'``

    EXAMPLES::

        sage: from sage.plot.colors import float_to_html
        sage: float_to_html(1.,1.,0.)
        '#ffff00'
        sage: float_to_html(.03,.06,.02)
        '#070f05'
        sage: float_to_html(*Color('brown').rgb())
        '#a52a2a'

    TESTS::

        sage: float_to_html((0.2, 0.6, 0.8))
        Traceback (most recent call last):
        ...
        TypeError: ...float_to_html() missing 2 required positional arguments: 'g' and 'b'
    """
    return "#%06x" % float_to_integer(r, g, b)


def float_to_integer(r, g, b):
    """
    Convert a Red-Green-Blue (RGB) color tuple to an integer.

    Each input value should be in the interval [0.0, 1.0]; otherwise,
    the values are first reduced modulo one (see :func:`mod_one`).

    INPUT:

    - ``r`` -- a real number; the RGB color's "red" intensity

    - ``g`` -- a real number; the RGB color's "green" intensity

    - ``b`` -- a real number; the RGB color's "blue" intensity

    OUTPUT:

    The integer 256^2 r_int + 256 g_int + b_int, where r_int, g_int, and
    b_int are obtained from r, g, and b by converting from the real
    interval [0.0, 1.0] to the integer range 0, 1, ..., 255.

    EXAMPLES::

        sage: from sage.plot.colors import float_to_integer
        sage: float_to_integer(1.,1.,0.)
        16776960
        sage: float_to_integer(.03,.06,.02)
        462597
        sage: float_to_integer(*Color('brown').rgb())
        10824234

    TESTS::

        sage: float_to_integer((0.2, 0.6, 0.8))
        Traceback (most recent call last):
        ...
        TypeError: ...float_to_integer() missing 2 required positional arguments: 'g' and 'b'
    """
    r, g, b = map(mod_one, (r, g, b))
    return int(r * 255) << 16 | int(g * 255) << 8 | int(b * 255)


def rainbow(n, format='hex'):
    """
    Return a list of colors sampled at equal intervals over the
    spectrum, from Hue-Saturation-Value (HSV) coordinates (0, 1, 1) to
    (1, 1, 1).  This range is red at the extremes, but it covers
    orange, yellow, green, cyan, blue, violet, and many other hues in
    between.  This function is particularly useful for representing
    vertex partitions on graphs.

    INPUT:

    - ``n`` -- a number; the length of the list

    - ``format`` -- string (default: ``'hex'``); the output format for
      each color in the list. The other choice is ``'rgbtuple'``.

    OUTPUT: a list of strings or RGB 3-tuples of floats in the interval
    [0.0, 1.0]

    EXAMPLES::

        sage: from sage.plot.colors import rainbow
        sage: rainbow(7)
        ['#ff0000', '#ffda00', '#48ff00', '#00ff91', '#0091ff', '#4800ff', '#ff00da']
        sage: rainbow(int(7))
        ['#ff0000', '#ffda00', '#48ff00', '#00ff91', '#0091ff', '#4800ff', '#ff00da']
        sage: rainbow(7, 'rgbtuple')
        [(1.0, 0.0, 0.0), (1.0, 0.8571428571428571, 0.0), (0.2857142857142858, 1.0, 0.0), (0.0, 1.0, 0.5714285714285712), (0.0, 0.5714285714285716, 1.0), (0.2857142857142856, 0.0, 1.0), (1.0, 0.0, 0.8571428571428577)]

    AUTHORS:

    - Robert L. Miller

    - Karl-Dieter Crisman (directly use :func:`hsv_to_rgb` for hues)
    """
    R = [tuple(map(float, hsv_to_rgb(i / n, 1, 1))) for i in range(n)]

    if format == 'rgbtuple':
        return R
    elif format == 'hex':
        for j in range(len(R)):
            R[j] = float_to_html(*R[j])
        return R


# If you change what this accepts, remember to change the documentation
# about cmap where it is used and to test these classes.
def get_cmap(cmap):
    r"""
    Return a color map (actually, a matplotlib :class:`Colormap`
    object), given its name or a [mixed] list/tuple of RGB list/tuples
    and color names.  For a list of map names, evaluate::

        sage: sorted(colormaps)
        ['Accent', ...]

    See :func:`rgbcolor` for valid list/tuple element formats.

    INPUT:

    - ``cmap`` -- string, list, tuple, or :class:`matplotlib.colors.Colormap`;
      a string must be a valid color map name

    OUTPUT: a :class:`matplotlib.colors.Colormap` instance

    EXAMPLES::

        sage: from sage.plot.colors import get_cmap
        sage: get_cmap('jet')
        <matplotlib.colors.LinearSegmentedColormap object at 0x...>
        sage: get_cmap([(0,0,0), (0.5,0.5,0.5), (1,1,1)])
        <matplotlib.colors.ListedColormap object at 0x...>
        sage: get_cmap(['green', 'lightblue', 'blue'])
        <matplotlib.colors.ListedColormap object at 0x...>
        sage: get_cmap(((0.5, 0.3, 0.2), [1.0, 0.0, 0.5], 'purple', Color(0.5,0.5,1, space='hsv')))
        <matplotlib.colors.ListedColormap object at 0x...>
        sage: get_cmap('jolies')
        Traceback (most recent call last):
        ...
        RuntimeError: Color map jolies not known (type "import matplotlib; list(matplotlib.colormaps.keys())" for valid names)
        sage: get_cmap('mpl')
        Traceback (most recent call last):
        ...
        RuntimeError: Color map mpl not known (type "import matplotlib; list(matplotlib.colormaps.keys())" for valid names)

    TESTS:

    Check that :issue:`33491` is fixed::

        sage: get_cmap('turbo')
        <matplotlib.colors.ListedColormap object at 0x...>
    """
    # matplotlib color maps
    from matplotlib.colors import ListedColormap, Colormap

    if isinstance(cmap, Colormap):
        return cmap

    elif isinstance(cmap, str):
        from matplotlib import colormaps
        try:
            return colormaps[cmap]
        except KeyError:
            raise RuntimeError('Color map %s not known (type "import matplotlib; list(matplotlib.colormaps.keys())" for valid names)' % cmap)

    elif isinstance(cmap, (list, tuple)):
        cmap = [rgbcolor(_) for _ in cmap]
        return ListedColormap(cmap)


def check_color_data(cfcm):
    """
    Make sure that the arguments are in order (coloring function, colormap).

    This will allow users to use both possible orders.

    EXAMPLES::

        sage: from sage.plot.colors import check_color_data
        sage: cf = lambda x,y : (x+y) % 1
        sage: cm = colormaps.autumn
        sage: check_color_data((cf, cm)) == (cf, cm)
        True
        sage: check_color_data((cm, cf)) == (cf, cm)
        True

    TESTS::

        sage: check_color_data(('a', 33))
        Traceback (most recent call last):
        ...
        ValueError: color data must be (color function, colormap)
    """
    cf, cm = cfcm
    from matplotlib.colors import Colormap
    if isinstance(cm, Colormap):
        return cf, cm
    elif isinstance(cf, Colormap):
        return cm, cf
    else:
        raise ValueError('color data must be (color function, colormap)')


class Colormaps(MutableMapping):
    """
    A dict-like collection of lazily-loaded matplotlib color maps.
    For a list of map names, evaluate::

        sage: sorted(colormaps)
        ['Accent', ...]
    """
    def __init__(self):
        """
        Construct an empty mutable collection of color maps.

        EXAMPLES::

            sage: from sage.plot.colors import Colormaps
            sage: maps = Colormaps()
            sage: len(maps.maps)
            0
        """
        self.maps = {}

    def load_maps(self):
        """
        If it's necessary, loads matplotlib's color maps and adds them
        to the collection.

        EXAMPLES::

            sage: from sage.plot.colors import Colormaps
            sage: maps = Colormaps()
            sage: len(maps.maps)
            0
            sage: maps.load_maps()
            sage: len(maps.maps)>60
            True
            sage: 'viridis' in maps
            True
        """
        if not self.maps:
            from matplotlib import colormaps
            self.maps.update(colormaps)

    def __dir__(self):
        """
        Return an approximate list of attribute names, including the
        color map names.

        OUTPUT: list of strings

        EXAMPLES::

            sage: from sage.plot.colors import Colormaps
            sage: maps = Colormaps()
            sage: 'Accent' in dir(maps)
            True
        """
        self.load_maps()
        methods = ['load_maps', '__dir__', '__len__', '__iter__',
                   '__contains__', '__getitem__', '__getattr__',
                   '__setitem__', '__delitem__']
        return dir(super()) + methods + sorted(self)

    def __len__(self):
        """
        Return the number of color maps.

        OUTPUT: integer

        EXAMPLES::

            sage: from sage.plot.colors import Colormaps
            sage: maps = Colormaps()
            sage: len(maps)>60
            True
        """
        self.load_maps()
        return len(self.maps)

    def __iter__(self):
        """
        Return an iterator over the color map collection.

        OUTPUT: a dictionary key iterator instance

        EXAMPLES::

            sage: from sage.plot.colors import Colormaps
            sage: maps = Colormaps()
            sage: count = 0
            sage: for m in maps: count += 1
            sage: count == len(maps)
            True
        """
        self.load_maps()
        return iter(self.maps)

    def __contains__(self, name):
        """
        Return whether a map is in the color maps collection.

        INPUT:

        - ``name`` -- string; the name of the map to query

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.plot.colors import Colormaps
            sage: maps = Colormaps()
            sage: 'summer' in maps
            True
            sage: 'not really a color map' in maps
            False
        """
        self.load_maps()
        return name in self.maps

    def __getitem__(self, name):
        """
        Get a color map from the collection via key access.

        INPUT:

        - ``name`` -- string; the name of the map return

        OUTPUT: an instance of :class:`matplotlib.colors.Colormap`

        EXAMPLES::

            sage: from sage.plot.colors import Colormaps
            sage: maps = Colormaps()
            sage: maps.get('Oranges')
            <matplotlib.colors.LinearSegmentedColormap object at ...>
            sage: maps['copper']
            <matplotlib.colors.LinearSegmentedColormap object at ...>
            sage: maps.get('not a color map')
            sage: maps['not a color map']
            Traceback (most recent call last):
            ...
            KeyError: "no colormap with name 'not a color map'"
        """
        self.load_maps()
        try:
            return self.maps[name]
        except KeyError:
            raise KeyError("no colormap with name '%s'" % name)

    def __getattr__(self, name):
        """
        Get a color map from the collection via attribute access.

        INPUT:

        - ``name`` -- string; the name of the map to return

        OUTPUT: an instance of :class:`matplotlib.colors.Colormap`

        EXAMPLES::

            sage: from sage.plot.colors import Colormaps
            sage: maps = Colormaps()
            sage: maps.pink
            <matplotlib.colors.LinearSegmentedColormap object at ...>
            sage: maps.punk
            Traceback (most recent call last):
            ...
            AttributeError: 'Colormaps' has no attribute or colormap punk...
            sage: maps['punk']
            Traceback (most recent call last):
            ...
            KeyError: "no colormap with name 'punk'"
            sage: maps['bone'] == maps.bone
            True
        """
        try:
            return self[name]
        except KeyError:
            raise AttributeError("'{}' has no attribute or colormap {}".format(type(self).__name__, name))

    def __repr__(self):
        """
        Return a string representation of the color map collection.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.plot.colors import Colormaps
            sage: maps = Colormaps()
            sage: maps
            {...}
            sage: type(repr(maps))
            <... 'str'>
        """
        self.load_maps()
        return repr(self.maps)

    def __setitem__(self, name, colormap):
        """
        Add a color map to the collection.

        INPUT:

        - ``name`` -- string; the name of the map to add

        - ``colormap`` -- an instance of
          :class:`matplotlib.colors.Colormap`; the color map to add

        EXAMPLES::

            sage: from sage.plot.colors import Colormaps, get_cmap
            sage: maps = Colormaps()
            sage: count = len(maps)
            sage: my_map = get_cmap(['chartreuse', '#007', (1.0, 0.0, 0.0)])
            sage: maps['my_map'] = my_map
            sage: 'my_map' in maps
            True
            sage: count + 1 == len(maps)
            True
        """
        self.load_maps()
        self.maps[name] = colormap

    def __delitem__(self, name):
        """
        Removes a color map from the collection.

        INPUT:

        - ``name`` -- string; the name of the map to remove

        EXAMPLES::

            sage: from sage.plot.colors import Colormaps
            sage: maps = Colormaps()
            sage: count = len(maps)
            sage: maps.popitem()  # random
            ('Spectral', <matplotlib.colors.LinearSegmentedColormap object at ...>)
            sage: count - 1 == len(maps)
            True
        """
        self.load_maps()
        del self.maps[name]


colormaps = Colormaps()
