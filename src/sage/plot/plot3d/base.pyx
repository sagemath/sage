r"""
Base classes for 3D graphics objects and plotting

The most important facts about these classes are
that you can simply add graphics objects
together (``G1+G2``, see :meth:`Graphics3d.__add__`),
and the :meth:`Graphics3d.show` method with its options for
choosing a viewer and setting
various parameters for displaying the graphics.

Most of the other methods of these classes are technical and
for special usage.

AUTHORS:

- Robert Bradshaw (2007-02): initial version

- Robert Bradshaw (2007-08): Cythonization, much optimization

- William Stein (2008)

- Paul Masson (2016): Three.js support

- Joshua Campbell (2020): Three.js animation support

- Günter Rote (2021): camera and light parameters for tachyon

.. TODO::

    finish integrating tachyon -- good default lights

    full documentation of three.js viewer parameters

    zoom by changing camera parameters instead of scaling objects
"""

# ****************************************************************************
#      Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cpython.list cimport *

import os
import zipfile

from functools import reduce
from io import StringIO
from random import randint

from sage.misc.temporary_file import tmp_filename
from sage.misc.fast_methods cimport hash_by_id
from sage.modules.free_module_element import vector
from sage.rings.real_double import RDF
from sage.plot.plot3d.texture import Texture
from sage.plot.plot3d.transform cimport Transformation, point_c, face_c
include "point_c.pxi"

from sage.interfaces.tachyon import tachyon_rt

from libc.math cimport INFINITY


default_texture = Texture()
pi = RDF.pi()


cdef class Graphics3d(SageObject):
    """
    This is the baseclass for all 3d graphics objects.

    .. automethod:: __add__
    .. automethod:: _rich_repr_
    """
    def __cinit__(self):
        """
        The Cython constructor.

        EXAMPLES::

            sage: gfx = sage.plot.plot3d.base.Graphics3d()
            sage: gfx._extra_kwds
            {}
        """
        self._extra_kwds = dict()

    def __hash__(self):
        r"""
        TESTS::

            sage: from sage.plot.plot3d.base import Graphics3d
            sage: hash(Graphics3d()) # random
            140658972348064
        """
        return hash_by_id(<void *> self)

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT: string

        EXAMPLES::

            sage: S = sphere((0, 0, 0), 1)
            sage: print(S)
            Graphics3d Object
        """
        return str(self)

    def _rich_repr_(self, display_manager, **kwds):
        """
        Rich Output Magic Method.

        See :mod:`sage.repl.rich_output` for details.

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: g = sphere()
            sage: g._rich_repr_(dm)  # OutputSceneThreejs container outside doctest mode
            OutputSceneJmol container
        """
        ### First, figure out the best graphics format
        types = display_manager.types
        can_view_jmol = (types.OutputSceneJmol in display_manager.supported_output())
        can_view_canvas3d  = (types.OutputSceneCanvas3d in display_manager.supported_output())
        can_view_wavefront = (types.OutputSceneWavefront in display_manager.supported_output())
        can_view_threejs = (types.OutputSceneThreejs in display_manager.supported_output())
        opts = self._process_viewing_options(kwds)
        viewer = opts.get('viewer', None)
        # make sure viewer is one of the supported options
        if viewer not in [None, 'jmol', 'tachyon', 'canvas3d', 'wavefront', 'threejs']:
            import warnings
            warnings.warn('viewer={} is not supported'.format(viewer))
            viewer = None
        # select suitable default
        if viewer is None:
            viewer = SHOW_DEFAULTS['viewer']
        # fall back to 2d image if necessary
        if viewer == 'canvas3d' and not can_view_canvas3d:
            viewer = 'jmol'
        if viewer == 'wavefront' and not can_view_wavefront:
            viewer = 'jmol'
        if viewer == 'threejs' and not can_view_threejs:
            viewer = 'jmol'
        if viewer == 'jmol' and not can_view_jmol:
            viewer = 'tachyon'
        # Second, return the corresponding graphics file
        if viewer == 'threejs':
            return self._rich_repr_threejs(**opts)
        elif viewer == 'jmol':
            return self._rich_repr_jmol(**opts)
        elif viewer == 'tachyon':
            preferred = (
                types.OutputImagePng,
                types.OutputImageJpg,
                types.OutputImageGif,
            )
            for output_container in preferred:
                if output_container in display_manager.supported_output():
                    return self._rich_repr_tachyon(output_container, **opts)
        elif viewer == 'canvas3d':
            return self._rich_repr_canvas3d(**opts)
        elif viewer == 'wavefront':
            return self._rich_repr_wavefront(**opts)
        else:
            assert False   # unreachable

    def _rich_repr_tachyon(self, output_container, **kwds):
        """
        Rich Representation using Tachyon.

        INPUT:

        - ``output_container`` -- the rich output container to contain
          the rendered image (can be png, gif, or jpg). Determines the
          type of the output.

        - ``**kwds`` -- optional keyword arguments are passed to the
          Tachyon raytracer

        OUTPUT:

        Instance of
        :class:`~sage.repl.rich_output.output_graphics.OutputImagePng`,
        :class:`~sage.repl.rich_output.output_graphics.OutputImageGif`, or
        :class:`~sage.repl.rich_output.output_graphics.OutputImageJpg`.

        EXAMPLES::

            sage: import sage.repl.rich_output.output_catalog as catalog
            sage: sphere()._rich_repr_tachyon(catalog.OutputImagePng)
            OutputImagePng container
            sage: import sage.repl.rich_output.output_catalog as catalog
            sage: sphere()._rich_repr_tachyon(catalog.OutputImageJpg)  # optional -- libjpeg
            OutputImageJpg container
        """
        filename = tmp_filename(ext='.png')
        opts = self._process_viewing_options(kwds)
        T = self._prepare_for_tachyon(
            opts['frame'], opts['axes'], opts['frame_aspect_ratio'],
            opts['aspect_ratio'],
            1 # opts['zoom']. Let zoom be handled by tachyon.
            # We don't want the perspective to change by zooming
        )

        tachyon_args = dict((key,val) for key,val in opts.items() if key in Graphics3d.tachyon_keywords)
        extra_opts = opts.get("extra_opts", "")
        if "shade" in opts:
            if opts["shade"] not in ["full", "medium", "low", "lowest"]:
                raise ValueError("shade must be set to 'full', 'medium', 'low' or 'lowest'")
            extra_opts += " -" + opts["shade"] + "shade"
        tachyon_rt(T.tachyon(**tachyon_args), filename, opts['verbosity'], extra_opts)
        from sage.repl.rich_output.buffer import OutputBuffer
        import sage.repl.rich_output.output_catalog as catalog
        import PIL.Image as Image
        if output_container is catalog.OutputImagePng:
            buf = OutputBuffer.from_file(filename)
        elif output_container is catalog.OutputImageGif:
            gif = tmp_filename(ext='.gif')
            Image.open(filename).save(gif)
            buf = OutputBuffer.from_file(gif)
        elif output_container is catalog.OutputImageJpg:
            jpg = tmp_filename(ext='.jpg')
            Image.open(filename).save(jpg)
            buf = OutputBuffer.from_file(jpg)
        else:
            raise ValueError('output_container not supported')
        return output_container(buf)

    def _rich_repr_jmol(self, **kwds):
        """
        Rich Representation as JMol scene.

        INPUT:

        - ``**kwds`` -- optional keyword arguments are passed to JMol

        OUTPUT:

        Instance of
        :class:`sage.repl.rich_output.output_graphics3d.OutputSceneJmol`.

        EXAMPLES::

            sage: sphere()._rich_repr_jmol()
            OutputSceneJmol container
        """
        from sage.misc.temporary_file import tmp_dir
        root_dir = os.path.abspath(tmp_dir())
        scene_zip     = os.path.join(root_dir, 'scene.spt.zip')
        preview_png   = os.path.join(root_dir, 'preview.png')
        opts = self._process_viewing_options(kwds)
        zoom = opts['zoom']
        T = self._prepare_for_jmol(
            opts['frame'],
            opts['axes'],
            opts['frame_aspect_ratio'],
            opts['aspect_ratio'],
            zoom,
        )
        T.export_jmol(scene_zip, **opts)
        from sage.interfaces.jmoldata import JmolData
        jdata = JmolData()
        if not jdata.is_jmol_available():
            # We can only use JMol to generate preview if a jvm is installed
            from sage.repl.rich_output.output_graphics import OutputImagePng
            tachyon = self._rich_repr_tachyon(OutputImagePng, **opts)
            tachyon.png.save_as(preview_png)
        else:
            # Java needs absolute paths
            scene_native = scene_zip

            script = '''set defaultdirectory "{0}"\nscript SCRIPT\n'''.format(scene_native)
            jdata.export_image(targetfile=preview_png, datafile=script,
                               image_type='PNG',
                               figsize=opts['figsize'][0])
        from sage.repl.rich_output.output_graphics3d import OutputSceneJmol
        from sage.repl.rich_output.buffer import OutputBuffer
        scene_zip     = OutputBuffer.from_file(scene_zip)
        preview_png   = OutputBuffer.from_file(preview_png)
        return OutputSceneJmol(scene_zip, preview_png)

    def _rich_repr_wavefront(self, **kwds):
        r"""
        Rich Representation as Wavefront (obj + mtl) Scene.

        INPUT:

        - ``**kwds`` -- optional keyword arguments are ignored

        OUTPUT:

        Instance of
        :class:`sage.repl.rich_output.output_graphics3d.OutputSceneWavefront`.

        EXAMPLES::

            sage: line = line3d([(0,0,0), (1,1,1)])
            sage: out = line._rich_repr_wavefront()
            sage: out
            OutputSceneWavefront container
            sage: out.obj.get_str()
            'mtllib ... 6 3 8 11\nf 8 7 12\nf 7 9 12\nf 9 10 12\nf 10 11 12\nf 11 8 12\n'
            sage: out.mtl.get_str()
            'newmtl texture...\nKd 0.4 0.4 1.0\nKs 0.0 0.0 0.0\nillum 1\nNs 1.0\nd 1.0\n'
        """
        from sage.repl.rich_output.output_graphics3d import OutputSceneWavefront
        from sage.repl.rich_output.buffer import OutputBuffer
        obj = OutputBuffer('mtllib scene.mtl\n' + self.obj())
        return OutputSceneWavefront(obj, self.mtl_str())

    def _rich_repr_canvas3d(self, **kwds):
        r"""
        Rich Representation as Canvas3D Scene.

        INPUT:

        - ``**kwds`` -- optional keyword arguments

        OUTPUT:

        Instance of
        :class:`sage.repl.rich_output.output_graphics3d.OutputSceneCanvas3d`.

        EXAMPLES::

            sage: out = sphere()._rich_repr_canvas3d()
            sage: out
            OutputSceneCanvas3d container
            sage: out.canvas3d.get_str()
            '[{"vertices":[{"x":0,"y":0,"z":-1},..., "color":"#6666ff", "opacity":1.0}]'
        """
        opts = self._process_viewing_options(kwds)
        aspect_ratio = opts['aspect_ratio'] # this necessarily has a value now
        frame_aspect_ratio = opts['frame_aspect_ratio']
        zoom = opts['zoom']
        frame = opts['frame']
        axes = opts['axes']
        T = self._prepare_for_tachyon(frame, axes, frame_aspect_ratio, aspect_ratio, zoom)
        data = flatten_list(T.json_repr(T.default_render_params()))
        canvas3d = '[' + ','.join(data) + ']'
        from sage.repl.rich_output.output_catalog import OutputSceneCanvas3d
        return OutputSceneCanvas3d(canvas3d)

    def _rich_repr_threejs(self, **kwds):
        r"""
        Rich Representation as Three.js Scene.

        INPUT:

        - ``**kwds`` -- optional keyword arguments

        OUTPUT:

        Instance of
        :class:`sage.repl.rich_output.output_graphics3d.OutputSceneThreejs`.

        EXAMPLES::

            sage: sphere(online=True)._rich_repr_threejs()
            OutputSceneThreejs container

        TESTS::

            sage: js = '// animation.js'
            sage: css = '/* animation.css */'
            sage: html = '<!-- animation.html -->'
            sage: d = dodecahedron()
            sage: i = icosahedron()
            sage: g1 = animate([d]).interactive()
            sage: g2 = animate([d, i]).interactive()

        Animation files are only included when at least 2 frames are present::

            sage: str = g1._rich_repr_threejs(online=True).html.get_str()
            sage: (js in str) or (css in str) or (html in str)
            False
            sage: str = g2._rich_repr_threejs(online=True).html.get_str()
            sage: (js in str) and (css in str) and (html in str)
            True

        Animation can be explicitly disabled by setting animate=False::

            sage: str = g2._rich_repr_threejs(online=True, animate=False).html.get_str()
            sage: (js in str) or (css in str) or (html in str)
            False

        Animation CSS and HTML are not included when animation_controls=False::

            sage: str = g2._rich_repr_threejs(online=True, animation_controls=False).html.get_str()
            sage: js in str
            True
            sage: (css in str) or (html in str)
            False

        "Fat" line scripts are only included when at least one line
        (or one surface with ``mesh=True``) exists with ``thickness > 1``::

            sage: fat = '// fat_lines.js'
            sage: L = line3d([(0, 0, 0), (1, 1, 1)], thickness=1)
            sage: str = L._rich_repr_threejs(online=True).html.get_str()
            sage: fat in str
            False
            sage: L = line3d([(0, 0, 0), (1, 1, 1)], thickness=10)
            sage: str = L._rich_repr_threejs(online=True).html.get_str()
            sage: fat in str
            True
            sage: d = dodecahedron(mesh=False, thickness=10)
            sage: str = d._rich_repr_threejs(online=True).html.get_str()
            sage: fat in str
            False
            sage: d = dodecahedron(mesh=True, thickness=1)
            sage: str = d._rich_repr_threejs(online=True).html.get_str()
            sage: fat in str
            False
            sage: d = dodecahedron(mesh=True, thickness=10)
            sage: str = d._rich_repr_threejs(online=True).html.get_str()
            sage: fat in str
            True

        If a page title is provided, it is stripped and HTML-escaped::

            sage: d = dodecahedron(page_title='\t"Page" & <Title>\n')
            sage: str = d._rich_repr_threejs(online=True).html.get_str()
            sage: '<title>&quot;Page&quot; &amp; &lt;Title&gt;</title>' in str
            True
        """
        options = self._process_viewing_options(kwds)
        options.setdefault('online', False)

        js_options = {} # options passed to Three.js template

        js_options['animate'] = options.get('animate', True)
        js_options['animationControls'] = options.get('animation_controls', True)
        js_options['aspectRatio'] = options.get('aspect_ratio', [1,1,1])
        js_options['autoScaling'] = options.get('auto_scaling', [False, False, False])
        js_options['autoPlay'] = options.get('auto_play', True)
        js_options['axes'] = options.get('axes', False)
        js_options['axesLabels'] = options.get('axes_labels', ['x','y','z'])
        js_options['axesLabelsStyle'] = options.get('axes_labels_style')
        js_options['decimals'] = options.get('decimals', 2)
        js_options['delay'] = options.get('delay', 20)
        js_options['frame'] = options.get('frame', True)
        js_options['loop'] = options.get('loop', True)
        js_options['projection'] = options.get('projection', 'perspective')
        js_options['theme'] = options.get('theme', 'light')
        js_options['viewpoint'] = options.get('viewpoint', False)

        if js_options['projection'] not in ['perspective', 'orthographic']:
            import warnings
            warnings.warn('projection={} is not supported; using perspective'.format(js_options['projection']))
            js_options['projection'] = 'perspective'

        if js_options['theme'] not in ['light', 'dark']:
            import warnings
            warnings.warn('theme={} is not supported; using light theme'.format(js_options['theme']))
            js_options['theme'] = 'light'

        # Normalization of options values for proper JSONing
        js_options['aspectRatio'] = [float(i) for i in js_options['aspectRatio']]
        js_options['decimals'] = int(js_options['decimals'])
        js_options['delay'] = int(js_options['delay'])

        if js_options['viewpoint']:
            if len(js_options['viewpoint']) != 2 or len(js_options['viewpoint'][0]) != 3:
                import warnings
                warnings.warn('viewpoint must be of the form [[x,y,z],angle]')
                js_options['viewpoint'] = False
            else:
                if type(js_options['viewpoint']) is tuple:
                    js_options['viewpoint'] = list(js_options['viewpoint'])
                if type(js_options['viewpoint'][0]) is tuple:
                    js_options['viewpoint'][0] = list(js_options['viewpoint'][0])
                js_options['viewpoint'][0] = [float(i) for i in js_options['viewpoint'][0]]
                js_options['viewpoint'][1] = float(js_options['viewpoint'][1])

        if not js_options['frame']:
            js_options['axesLabels'] = False
            js_options['axesLabelsStyle'] = None

        if js_options['axesLabelsStyle'] is not None:
            from sage.plot.plot3d.shapes import _validate_threejs_text_style
            style = js_options['axesLabelsStyle']
            if isinstance(style, dict):
                style = _validate_threejs_text_style(style)
                style = [style, style, style]
            elif isinstance(style, list) and len(style) == 3 and all([isinstance(s, dict) for s in style]):
                style = [_validate_threejs_text_style(s) for s in style]
            else:
                import warnings
                warnings.warn("axes_labels_style must be a dict or a list of 3 dicts")
                style = [dict(), dict(), dict()]
            js_options['axesLabelsStyle'] = style

        from sage.repl.rich_output import get_display_manager
        scripts = get_display_manager().threejs_scripts(options['online'])
        styles = ''
        extra_html = ''

        b = self.bounding_box()
        bounds = '[{{"x":{}, "y":{}, "z":{}}}, {{"x":{}, "y":{}, "z":{}}}]'.format(
                 b[0][0], b[0][1], b[0][2], b[1][0], b[1][1], b[1][2])

        from sage.plot.colors import Color
        lights = '[{{"x":-5, "y":3, "z":0, "color":"{}", "parent":"camera"}}]'.format(
                 Color(.5,.5,.5).html_color())
        ambient = '{{"color":"{}"}}'.format(Color(.5,.5,.5).html_color())

        import json

        reprs = {'point': [], 'line': [], 'text': [], 'surface': []}
        frame_count = 0
        fat_lines = False
        for kind, desc in self.threejs_repr(self.default_render_params()):
            reprs[kind].append(desc)
            keyframe = int(desc.get('keyframe', -1))
            frame_count = max(frame_count, keyframe + 1)
            if kind == 'line' or (kind == 'surface' and desc.get('showMeshGrid')):
                linewidth = float(desc.get('linewidth', 1))
                if linewidth > 1:
                    fat_lines = True
        reprs = {kind: json.dumps(descs) for kind, descs in reprs.items()}

        from sage.env import SAGE_EXTCODE
        with open(os.path.join(
                SAGE_EXTCODE, 'threejs', 'threejs_template.html')) as f:
            html = f.read()

        if fat_lines:
            with open(os.path.join(SAGE_EXTCODE, 'threejs', 'fat_lines.js')) as f:
                scripts += '<script>' + f.read() + '</script>'

        js_options['animate'] = js_options['animate'] and frame_count > 1
        if js_options['animate']:
            if js_options['animationControls']:
                with open(os.path.join(SAGE_EXTCODE, 'threejs', 'animation.css')) as f:
                    css = f.read()
                    css = css.replace('SAGE_FRAME_COUNT', str(frame_count))
                    styles += '<style>' + css + '</style>'
                with open(os.path.join(SAGE_EXTCODE, 'threejs', 'animation.html')) as f:
                    extra_html += f.read()
            with open(os.path.join(SAGE_EXTCODE, 'threejs', 'animation.js')) as f:
                extra_html += '<script>' + f.read() + '</script>'

        page_title = options.get('page_title')
        if page_title is None:
            page_title = ""
        else:
            from html import escape as html_escape
            page_title = html_escape(str(page_title).strip())

        html = html.replace('SAGE_TITLE', page_title)
        html = html.replace('SAGE_SCRIPTS', scripts)
        html = html.replace('SAGE_STYLES', styles)
        html = html.replace('SAGE_EXTRA_HTML', extra_html)
        html = html.replace('SAGE_OPTIONS', json.dumps(js_options))
        html = html.replace('SAGE_BOUNDS', bounds)
        html = html.replace('SAGE_LIGHTS', lights)
        html = html.replace('SAGE_AMBIENT', ambient)
        html = html.replace('SAGE_TEXTS', str(reprs['text']))
        html = html.replace('SAGE_POINTS', str(reprs['point']))
        html = html.replace('SAGE_LINES', str(reprs['line']))
        html = html.replace('SAGE_SURFACES', str(reprs['surface']))

        from sage.repl.rich_output.output_catalog import OutputSceneThreejs
        return OutputSceneThreejs(html)

    def __str__(self):
        """
        EXAMPLES::

            sage: S = sphere((0, 0, 0), 1)
            sage: str(S)
            'Graphics3d Object'
        """
        return "Graphics3d Object"

    def __add__(left, right):
        """
        Addition of objects adds them to the same scene.

        EXAMPLES::

            sage: A = sphere((0,0,0), 1, color='red')
            sage: B = dodecahedron((2, 0, 0), color='yellow')
            sage: A+B
            Graphics3d Object

        For convenience, we take 0 and ``None`` to be the additive identity::

            sage: A + 0 is A
            True
            sage: A + None is A, 0 + A is A, None + A is A
            (True, True, True)

        In particular, this allows us to use the sum() function without
        having to provide an empty starting object::

            sage: sum(point3d((cos(n), sin(n), n)) for n in [0..10, step=.1])
            Graphics3d Object

        A Graphics 3d object and a 2d object can also be added::

            sage: A = sphere((0, 0, 0), 1) + circle((0, 0), 1.5)
            sage: A.show(aspect_ratio=1)
        """
        if right == 0 or right is None:
            return left
        elif left == 0 or left is None:
            return right
        elif not isinstance(left, Graphics3d):
            left = left.plot3d()
        elif not isinstance(right, Graphics3d):
            right = right.plot3d()
        return Graphics3dGroup([left, right])

    def _set_extra_kwds(self, kwds):
        """
        Allow one to pass rendering arguments on as if they were set
        in the constructor.

        EXAMPLES::

            sage: S = sphere((0, 0, 0), 1)
            sage: S._set_extra_kwds({'aspect_ratio': [1, 2, 2]})
            sage: S
            Graphics3d Object
        """
        self._extra_kwds = kwds

    def aspect_ratio(self, v=None):
        """
        Set or get the preferred aspect ratio.

        INPUT:

        - ``v`` -- (default: ``None``) must be a list or tuple of length three,
          or the integer ``1``. If no arguments are provided then the
          default aspect ratio is returned.

        EXAMPLES::

            sage: D = dodecahedron()
            sage: D.aspect_ratio()
            [1.0, 1.0, 1.0]
            sage: D.aspect_ratio([1,2,3])
            sage: D.aspect_ratio()
            [1.0, 2.0, 3.0]
            sage: D.aspect_ratio(1)
            sage: D.aspect_ratio()
            [1.0, 1.0, 1.0]
        """
        if v is not None:
            if v == 1:
                v = (1, 1, 1)
            if not isinstance(v, (tuple, list)):
                raise TypeError("aspect_ratio must be a list or tuple of "
                                "length 3 or the integer 1")
            self._aspect_ratio = [float(z) for z in v]
        else:
            if self._aspect_ratio is None:
                self._aspect_ratio = [1.0, 1.0, 1.0]
            return self._aspect_ratio

    def frame_aspect_ratio(self, v=None):
        """
        Set or get the preferred frame aspect ratio.

        INPUT:

        - ``v`` -- (default: ``None``) must be a list or tuple of
          length three, or the integer ``1``. If no arguments are
          provided then the default frame aspect ratio is returned.

        EXAMPLES::

            sage: D = dodecahedron()
            sage: D.frame_aspect_ratio()
            [1.0, 1.0, 1.0]
            sage: D.frame_aspect_ratio([2,2,1])
            sage: D.frame_aspect_ratio()
            [2.0, 2.0, 1.0]
            sage: D.frame_aspect_ratio(1)
            sage: D.frame_aspect_ratio()
            [1.0, 1.0, 1.0]
        """
        if v is not None:
            if v == 1:
                v = (1, 1, 1)
            if not isinstance(v, (tuple, list)):
                raise TypeError("frame_aspect_ratio must be a list or tuple of "
                                "length 3 or the integer 1")
            self._frame_aspect_ratio = [float(z) for z in v]
        else:
            if self._frame_aspect_ratio is None:
                self._frame_aspect_ratio = [1.0, 1.0, 1.0]
            return self._frame_aspect_ratio

    def _determine_frame_aspect_ratio(self, aspect_ratio):
        a_min, a_max = self._safe_bounding_box()
        return [(a_max[i] - a_min[i])*aspect_ratio[i] for i in range(3)]

    def _safe_bounding_box(self):
        """
        Return a bounding box but where no side length is 0.

        This is used to avoid zero-division errors for pathological
        plots.

        EXAMPLES::

            sage: G = line3d([(0, 0, 0), (0, 0, 1)])
            sage: G.bounding_box()
            ((0.0, 0.0, 0.0), (0.0, 0.0, 1.0))
            sage: G._safe_bounding_box()
            ([-1.0, -1.0, 0.0], [1.0, 1.0, 1.0])
        """
        a_min, a_max = self.bounding_box()
        a_min = list(a_min)
        a_max = list(a_max)
        for i in range(3):
            if a_min[i] == a_max[i]:
                a_min[i] = a_min[i] - 1
                a_max[i] = a_max[i] + 1
        return a_min, a_max

    def bounding_box(self):
        """
        Return the lower and upper corners of a 3d bounding box.

        This is used for rendering, and the scene should fit entirely
        within this box.

        Specifically, the first point returned has x, y, and z
        coordinates that are the respective minimum over all points
        in the graphics, and the second point is the maximum.

        The default return value is simply the box containing the origin.

        EXAMPLES::

            sage: sphere((1,1,1), 2).bounding_box()
            ((-1.0, -1.0, -1.0), (3.0, 3.0, 3.0))
            sage: G = line3d([(1, 2, 3), (-1,-2,-3)])
            sage: G.bounding_box()
            ((-1.0, -2.0, -3.0), (1.0, 2.0, 3.0))
        """
        return ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))

    def transform(self, **kwds):
        """
        Apply a transformation, where the inputs are
        passed onto a TransformGroup object.

        Mostly for internal use; see the translate, scale, and rotate
        methods for more details.

        EXAMPLES::

            sage: sphere((0,0,0), 1).transform(trans=(1, 0, 0), scale=(2,3,4)).bounding_box()
            ((-1.0, -3.0, -4.0), (3.0, 3.0, 4.0))
        """
        return TransformGroup([self], **kwds)

    def translate(self, *x):
        """
        Return the object translated by the given vector (which can be
        given either as a 3-iterable or via positional arguments).

        EXAMPLES::

            sage: icosahedron() + sum(icosahedron(opacity=0.25).translate(2*n, 0, 0) for n in [1..4])
            Graphics3d Object
            sage: icosahedron() + sum(icosahedron(opacity=0.25).translate([-2*n, n, n^2]) for n in [1..4])
            Graphics3d Object

        TESTS::

            sage: G = sphere((0, 0, 0), 1)
            sage: G.bounding_box()
            ((-1.0, -1.0, -1.0), (1.0, 1.0, 1.0))
            sage: G.translate(0, 0, 1).bounding_box()
            ((-1.0, -1.0, 0.0), (1.0, 1.0, 2.0))
            sage: G.translate(-1, 5, 0).bounding_box()
            ((-2.0, 4.0, -1.0), (0.0, 6.0, 1.0))
        """
        if len(x) == 1:
            x = x[0]
        return self.transform(trans=x)

    def scale(self, *x):
        """
        Return the object scaled in the x, y, and z directions.

        EXAMPLES::

            sage: G = dodecahedron() + dodecahedron(opacity=.5).scale(2)
            sage: G.show(aspect_ratio=1)
            sage: G = icosahedron() + icosahedron(opacity=.5).scale([1, 1/2, 2])
            sage: G.show(aspect_ratio=1)

        TESTS::

            sage: G = sphere((0, 0, 0), 1)
            sage: G.scale(2)
            Graphics3d Object
            sage: G.scale(1, 2, 1/2).show(aspect_ratio=1)
            sage: G.scale(2).bounding_box()
            ((-2.0, -2.0, -2.0), (2.0, 2.0, 2.0))
        """
        if isinstance(x[0], (tuple, list)):
            x = x[0]
        return self.transform(scale=x)

    def rotate(self, v, theta):
        r"""
        Return the object rotated about the vector `v` by `\theta` radians.

        EXAMPLES::

            sage: from math import pi
            sage: from sage.plot.plot3d.shapes import Cone
            sage: v = (1,2,3)
            sage: G = arrow3d((0, 0, 0), v)
            sage: G += Cone(1/5, 1).translate((0, 0, 2))
            sage: C = Cone(1/5, 1, opacity=.25).translate((0, 0, 2))
            sage: G += sum(C.rotate(v, pi*t/4) for t in [1..7])
            sage: G.show(aspect_ratio=1)

            sage: from sage.plot.plot3d.shapes import Box
            sage: Box(1/3, 1/5, 1/7).rotate((1, 1, 1), pi/3).show(aspect_ratio=1)
        """
        vx, vy, vz = v
        return self.transform(rot=[vx, vy, vz, theta])

    def rotateX(self, theta):
        """
        Return the object rotated about the `x`-axis by the given angle.

        EXAMPLES::

            sage: from math import pi
            sage: from sage.plot.plot3d.shapes import Cone
            sage: G = Cone(1/5, 1) + Cone(1/5, 1, opacity=.25).rotateX(pi/2)
            sage: G.show(aspect_ratio=1)
        """
        return self.rotate((1, 0, 0), theta)

    def rotateY(self, theta):
        """
        Return the object rotated about the `y`-axis by the given angle.

        EXAMPLES::

            sage: from math import pi
            sage: from sage.plot.plot3d.shapes import Cone
            sage: G = Cone(1/5, 1) + Cone(1/5, 1, opacity=.25).rotateY(pi/3)
            sage: G.show(aspect_ratio=1)
        """
        return self.rotate((0, 1, 0), theta)

    def rotateZ(self, theta):
        """
        Return the object rotated about the `z`-axis by the given angle.

        EXAMPLES::

            sage: from math import pi
            sage: from sage.plot.plot3d.shapes import Box
            sage: G = Box(1/2, 1/3, 1/5) + Box(1/2, 1/3, 1/5, opacity=.25).rotateZ(pi/5)
            sage: G.show(aspect_ratio=1)
        """
        return self.rotate((0, 0, 1), theta)

    def viewpoint(self):
        """
        Return the viewpoint of this plot.

        Currently only a stub for x3d.

        EXAMPLES::

            sage: type(dodecahedron().viewpoint())
            <class 'sage.plot.plot3d.base.Viewpoint'>
        """
        # This should probably be reworked somehow.
        return Viewpoint(0,0,6)

    def default_render_params(self):
        """
        Return an instance of RenderParams suitable for plotting this object.

        EXAMPLES::

            sage: type(dodecahedron().default_render_params())
            <class 'sage.plot.plot3d.base.RenderParams'>
        """
        return RenderParams(ds=.075)

    def testing_render_params(self):
        """
        Return an instance of RenderParams suitable for testing this object.

        In particular, it opens up a temporary file as an auxiliary zip
        file for jmol.

        EXAMPLES::

            sage: type(dodecahedron().testing_render_params())
            <class 'sage.plot.plot3d.base.RenderParams'>
        """
        params = RenderParams(ds=.075)
        fn = tmp_filename(ext='.zip')
        params.output_archive = zipfile.ZipFile(fn, 'w', zipfile.ZIP_STORED, True)
        return params

    def x3d(self):
        """
        An x3d scene file (as a string) containing the this object.

        EXAMPLES::

            sage: print(sphere((1, 2, 3), 5).x3d())
            <X3D version='3.0' profile='Immersive' xmlns:xsd='http://www.w3.org/2001/XMLSchema-instance' xsd:noNamespaceSchemaLocation=' http://www.web3d.org/specifications/x3d-3.0.xsd '>
            <head>
            <meta name='title' content='sage3d'/>
            </head>
            <Scene>
            <Viewpoint position='0 0 6'/>
            <Transform translation='1 2 3'>
            <Shape><Sphere radius='5.0'/><Appearance><Material diffuseColor='0.4 0.4 1.0' shininess='1.0' specularColor='0.0 0.0 0.0'/></Appearance></Shape>
            </Transform>
            </Scene>
            </X3D>

            sage: G = icosahedron() + sphere((0,0,0), 0.5, color='red')
            sage: print(G.x3d())
            <X3D version='3.0' profile='Immersive' xmlns:xsd='http://www.w3.org/2001/XMLSchema-instance' xsd:noNamespaceSchemaLocation=' http://www.web3d.org/specifications/x3d-3.0.xsd '>
            <head>
            <meta name='title' content='sage3d'/>
            </head>
            <Scene>
            <Viewpoint position='0 0 6'/>
            <Shape>
            <IndexedFaceSet coordIndex='...'>
              <Coordinate point='...'/>
            </IndexedFaceSet>
            <Appearance><Material diffuseColor='0.4 0.4 1.0' shininess='1.0' specularColor='0.0 0.0 0.0'/></Appearance></Shape>
            <Transform translation='0 0 0'>
            <Shape><Sphere radius='0.5'/><Appearance><Material diffuseColor='1.0 0.0 0.0' shininess='1.0' specularColor='0.0 0.0 0.0'/></Appearance></Shape>
            </Transform>
            </Scene>
            </X3D>
        """
        return """
<X3D version='3.0' profile='Immersive' xmlns:xsd='http://www.w3.org/2001/XMLSchema-instance' xsd:noNamespaceSchemaLocation=' http://www.web3d.org/specifications/x3d-3.0.xsd '>
<head>
<meta name='title' content='sage3d'/>
</head>
<Scene>
%s
%s
</Scene>
</X3D>
""" % (self.viewpoint().x3d_str(), self.x3d_str())

    ################ TACHYON ################

    ####### insertion of camera parameters

    tachyon_keywords = (
      "antialiasing",
      # "aspectratio",
      "zoom",  # zoom was previously handled directly by scaling the scene.
               # This has now been disabled, and zoom is handled by tachyon.
      "raydepth", "figsize", "light_position",
      "camera_position","updir",
      # "look_at", # omit look_at. viewdir is sufficient for most purposes
      "viewdir")

    # The tachyon "aspectratio" parameter is outdated for normal users:
    # From the tachyon documentation:
    # "By using the aspect ratio parameter, one can produce images which look
    # correct on any screen.  Aspect ratio alters the relative width of the image,
    # while keeping plane the height of the image plane constant.  In general,
    # most workstation displays have an aspect ratio of 1.0."

    # Instead, the tachyon aspectratio is set to match nonsquare
    # drawing area in "figsize".

    # Parameters are mostly taken from tachyon.py,
    # but camera_center is renamed camera_position.
    # Apparently reST strips () from default parameters in the automatic documentation.
    # Thus, I replaced () by [] as default values.

    def tachyon(self,
        zoom=1.0,
        antialiasing=False,
        figsize=[5,5], # resolution = 100*figsize
        raydepth=8,
        camera_position=[2.3, 2.4, 2.0], # old default values
        updir=[0, 0, 1],
        # look_at=(0, 0, 0), # could be nice to have, but viewdir is good enough
        light_position=[4.0, 3.0, 2.0],
        viewdir=None,
        # projection='PERSPECTIVE', # future extension, allow different projection types
    ):
        """
        A tachyon input file (as a string) containing the this object.

        EXAMPLES::

            sage: print(sphere((1, 2, 3), 5, color='yellow').tachyon())
            <BLANKLINE>
            begin_scene
            resolution 500 500
            <BLANKLINE>
                     camera
                    ...
                  plane
                    center -592.870151560437 618.647114671761 -515.539262226467
                    normal -2.3 2.4 -2.0
                    TEXTURE
                        AMBIENT 1.0 DIFFUSE 0.0 SPECULAR 0.0 OPACITY 1.0
                        COLOR 1.0 1.0 1.0
                        TEXFUNC 0
            <BLANKLINE>
                Texdef texture...
              Ambient 0.3333333333333333 Diffuse 0.6666666666666666 Specular 0.0 Opacity 1.0
              Color 1.0 1.0 0.0
              TexFunc 0
            <BLANKLINE>
                Sphere center 1.0 -2.0 3.0 Rad 5.0 texture...
            <BLANKLINE>
            end_scene

            sage: G = icosahedron(color='red') + sphere((1,2,3), 0.5, color='yellow')
            sage: G.show(viewer='tachyon', frame=false)
            sage: print(G.tachyon())
            begin_scene
            ...
            Texdef texture...
              Ambient 0.3333333333333333 Diffuse 0.6666666666666666 Specular 0.0 Opacity 1.0
               Color 1.0 1.0 0.0
               TexFunc 0
            TRI V0 ...
            Sphere center 1.0 -2.0 3.0 Rad 0.5 texture...
            end_scene
        """
        render_params = self.default_render_params()
        # switch from LH to RH coords to be consistent with java rendition
        render_params.push_transform(Transformation(scale=[1,-1,1]))

        if len(camera_position)!=3:
            raise ValueError('Camera center must consist of three numbers')

        if viewdir is None:
            viewdir = [float(- camera_position[i]) for i in range(3)]
            if viewdir == [0.0,0.0,0.0]:
                viewdir = (1,0,0) # issue a Warning? "camera_position at origin"
        # switch from LH to RH coords to be consistent with java rendition
        viewdir = _flip_orientation(viewdir)
        updir = _flip_orientation(updir)
        camera_position = _flip_orientation(camera_position)
        light_position = _flip_orientation(light_position)

        return """
begin_scene
resolution {resolution_x:d} {resolution_y:d}

         camera
            zoom {zoom:f}
            aspectratio {aspectratio:f}
            antialiasing {antialiasing:d}
            raydepth {raydepth:d}
            center {camera_position}
            viewdir {viewdir}
            updir {updir}
         end_camera

      light center {light_position}
            rad 0.2
            color  1.0 1.0 1.0

      plane
        center {viewdir1000}
        normal {viewdir}
        TEXTURE
            AMBIENT 1.0 DIFFUSE 0.0 SPECULAR 0.0 OPACITY 1.0
            COLOR 1.0 1.0 1.0
            TEXFUNC 0

    {scene}

    {render_parameters}

end_scene""".format(
    # render_params.antialiasing, this only provided the default value of 8
    scene = "\n".join(sorted([t.tachyon_str() for t in self.texture_set()])),
    render_parameters =
             "\n".join(flatten_list(self.tachyon_repr(render_params))),
    viewdir1000=self._tostring(1000*vector(viewdir).normalized().n()),
    viewdir=self._tostring(viewdir),
    camera_position=self._tostring(camera_position),
    updir=self._tostring(updir),
    light_position=self._tostring(light_position),
    zoom=zoom,
    antialiasing=antialiasing,
    resolution_x=figsize[0]*100,
    resolution_y=figsize[1]*100,
    aspectratio=float(figsize[1])/float(figsize[0]),
    raydepth=raydepth,
                   )

    def obj(self):
        """
        An .obj scene file (as a string) containing the this object.

        A .mtl file of the same name must also be produced for
        coloring.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import ColorCube
            sage: print(ColorCube(1, ['red', 'yellow', 'blue']).obj())
            g obj_1
            usemtl ...
            v 1 1 1
            v -1 1 1
            v -1 -1 1
            v 1 -1 1
            f 1 2 3 4
            ...
            g obj_6
            usemtl ...
            v -1 -1 1
            v -1 1 1
            v -1 1 -1
            v -1 -1 -1
            f 21 22 23 24
        """
        return "\n".join(flatten_list([self.obj_repr(self.default_render_params()), ""]))

    @staticmethod
    def _tostring(s):
        r"""
        Convert vector information to a space-separated string.

        EXAMPLES::

            sage: sage.plot.plot3d.base.Graphics3d._tostring((1.0,1.2,-1.3))
            '1.00000000000000 1.20000000000000 -1.30000000000000'
        """
        return ' '.join(map(str,s))

    def export_jmol(self, filename='jmol_shape.jmol', force_reload=False,
                    zoom=1, spin=False, background=(1,1,1), stereo=False,
                    mesh=False, dots=False,
                    perspective_depth = True,
                    orientation = (-764,-346,-545,76.39), **ignored_kwds):
                    # orientation chosen to look same as tachyon
        """
        A jmol scene consists of a script which refers to external files.
        Fortunately, we are able to put all of them in a single zip archive,
        which is the output of this call.

        EXAMPLES::

            sage: out_file = tmp_filename(ext='.jmol')
            sage: G = sphere((1, 2, 3), 5) + cube() + sage.plot.plot3d.shapes.Text("hi")
            sage: G.export_jmol(out_file)
            sage: import zipfile
            sage: z = zipfile.ZipFile(out_file)
            sage: z.namelist()
            ['obj_...pmesh', 'SCRIPT']

            sage: print(z.read('SCRIPT').decode('ascii'))
            data "model list"
            2
            empty
            Xx 0 0 0
            Xx 5.5 5.5 5.5
            end "model list"; show data
            select *
            wireframe off; spacefill off
            set labelOffset 0 0
            background [255,255,255]
            spin OFF
            moveto 0 -764 -346 -545 76.39
            centerAt absolute {0 0 0}
            zoom 100
            frank OFF
            set perspectivedepth ON
            isosurface sphere_1  center {1.0 2.0 3.0} sphere 5.0
            color isosurface  [102,102,255]
            pmesh obj_... "obj_...pmesh"
            color pmesh  [102,102,255]
            select atomno = 1
            color atom  [102,102,255]
            label "hi"
            isosurface fullylit; pmesh o* fullylit; set antialiasdisplay on;

            sage: print(z.read(z.namelist()[0]).decode('ascii'))
            24
            0.5 0.5 0.5
            -0.5 0.5 0.5
            ...
            -0.5 -0.5 -0.5
            6
            5
            0
            1
            ...
        """
        render_params = self.default_render_params()
        render_params.mesh = mesh
        render_params.dots = dots
        render_params.output_file = filename
        render_params.force_reload = render_params.randomize_counter = force_reload
        render_params.output_archive = zipfile.ZipFile(filename, 'w', zipfile.ZIP_DEFLATED, True)
        # Render the data
        all = flatten_list([self.jmol_repr(render_params), ""])

        f = StringIO()

        if render_params.atom_list:
            # Load the atom model
            f.write('data "model list"\n')
            f.write('%s\nempty\n' % (len(render_params.atom_list) + 1))
            for atom in render_params.atom_list:
                f.write('Xx %s %s %s\n' % atom)
            f.write('Xx 5.5 5.5 5.5\n') # so the zoom fits the box
            f.write('end "model list"; show data\n')
            f.write('select *\n')
            f.write('wireframe off; spacefill off\n')
            f.write('set labelOffset 0 0\n')

        # Set the scene background color
        f.write('background [%s,%s,%s]\n' % tuple(int(a*255) for a in background))
        if spin:
            f.write('spin ON\n')
        else:
            f.write('spin OFF\n')
        if stereo:
            if stereo is True:
                stereo = "redblue"
            f.write('stereo %s\n' % stereo)
        if orientation:
            f.write('moveto 0 %s %s %s %s\n' % tuple(orientation))

        f.write('centerAt absolute {0 0 0}\n')
        f.write('zoom {0}\n'.format(zoom * 100))
        f.write('frank OFF\n')  # jmol logo

        if perspective_depth:
            f.write('set perspectivedepth ON\n')
        else:
            f.write('set perspectivedepth OFF\n')

        # Put the rest of the object in
        f.write("\n".join(all))
        # Make sure the lighting is correct
        f.write("isosurface fullylit; pmesh o* fullylit; set antialiasdisplay on;\n")

        render_params.output_archive.writestr('SCRIPT', f.getvalue())
        render_params.output_archive.close()

    def json_repr(self, render_params):
        """
        A (possibly nested) list of strings. Each entry is formatted
        as JSON, so that a JavaScript client could eval it and get an
        object. Each object has fields to encapsulate the faces and
        vertices of the object. This representation is intended to be
        consumed by the canvas3d viewer backend.

        EXAMPLES::

            sage: G = sage.plot.plot3d.base.Graphics3d()
            sage: G.json_repr(G.default_render_params())
            []
        """
        return []

    def jmol_repr(self, render_params):
        r"""
        A (possibly nested) list of strings which will be concatenated and
        used by jmol to render the object.

        (Nested lists of strings are used because otherwise all the
        intermediate concatenations can kill performance). This may
        refer to several remove files, which are stored in
        render_parames.output_archive.

        EXAMPLES::

            sage: G = sage.plot.plot3d.base.Graphics3d()
            sage: G.jmol_repr(G.default_render_params())
            []
            sage: G = sphere((1, 2, 3))
            sage: G.jmol_repr(G.default_render_params())
            [['isosurface sphere_1  center {1.0 2.0 3.0} sphere 1.0\ncolor isosurface  [102,102,255]']]
        """
        return []

    def tachyon_repr(self, render_params):
        r"""
        A (possibly nested) list of strings which will be concatenated and
        used by tachyon to render the object.

        (Nested lists of strings are used because otherwise all the
        intermediate concatenations can kill performance). This may
        include a reference to color information which is stored
        elsewhere.

        EXAMPLES::

            sage: G = sage.plot.plot3d.base.Graphics3d()
            sage: G.tachyon_repr(G.default_render_params())
            []
            sage: G = sphere((1, 2, 3))
            sage: G.tachyon_repr(G.default_render_params())
            ['Sphere center 1.0 2.0 3.0 Rad 1.0 texture...']
        """
        return []

    def obj_repr(self, render_params):
        """
        A (possibly nested) list of strings which will be concatenated and
        used to construct an .obj file of the object.

        (Nested lists of strings are used because otherwise all the
        intermediate concatenations can kill performance). This may
        include a reference to color information which is stored
        elsewhere.

        EXAMPLES::

            sage: G = sage.plot.plot3d.base.Graphics3d()
            sage: G.obj_repr(G.default_render_params())
            []
            sage: G = cube()
            sage: G.obj_repr(G.default_render_params())
            ['g obj_1',
             'usemtl ...',
             ['v 0.5 0.5 0.5',
              'v -0.5 0.5 0.5',
              'v -0.5 -0.5 0.5',
              'v 0.5 -0.5 0.5',
              'v 0.5 0.5 -0.5',
              'v -0.5 0.5 -0.5',
              'v 0.5 -0.5 -0.5',
              'v -0.5 -0.5 -0.5'],
             ['f 1 2 3 4',
              'f 1 5 6 2',
              'f 1 4 7 5',
              'f 6 5 7 8',
              'f 7 4 3 8',
              'f 3 2 6 8'],
             []]
        """
        return []

    def threejs_repr(self, render_params):
        """
        A flat list of ``(kind, desc)`` tuples where ``kind`` is one of:
        'point', 'line', 'text', or 'surface'; and where ``desc`` is a dictionary
        describing a point, line, text, or surface.

        EXAMPLES::

            sage: G = sage.plot.plot3d.base.Graphics3d()
            sage: G.threejs_repr(G.default_render_params())
            []
        """
        return []

    def texture_set(self):
        """
        Often the textures of a 3d file format are kept separate from the
        objects themselves. This function returns the set of textures used,
        so they can be defined in a preamble or separate file.

        EXAMPLES::

            sage: sage.plot.plot3d.base.Graphics3d().texture_set()
            set()

            sage: G = tetrahedron(color='red') + tetrahedron(color='yellow') + tetrahedron(color='red', opacity=0.5)
            sage: [t for t in G.texture_set() if t.color == colors.red] # we should have two red textures
            [Texture(texture..., red, ff0000), Texture(texture..., red, ff0000)]
            sage: [t for t in G.texture_set() if t.color == colors.yellow] # ...and one yellow
            [Texture(texture..., yellow, ffff00)]
        """
        return set()

    def mtl_str(self):
        """
        Return the contents of a .mtl file, to be used to provide coloring
        information for an .obj file.

        EXAMPLES::

            sage: G = tetrahedron(color='red') + tetrahedron(color='yellow', opacity=0.5)
            sage: print(G.mtl_str())
            newmtl ...
            Ka 0.5 5e-06 5e-06
            Kd 1.0 1e-05 1e-05
            Ks 0.0 0.0 0.0
            illum 1
            Ns 1.0
            d 1.0
            newmtl ...
            Ka 0.5 0.5 5e-06
            Kd 1.0 1.0 1e-05
            Ks 0.0 0.0 0.0
            illum 1
            Ns 1.0
            d 0.5
        """
        return "\n\n".join(sorted([t.mtl_str() for t in self.texture_set()])) + "\n"

    def flatten(self):
        """
        Try to reduce the depth of the scene tree by consolidating groups
        and transformations.

        The generic Graphics3d object cannot be made flatter.

        EXAMPLES::

            sage: G = sage.plot.plot3d.base.Graphics3d()
            sage: G.flatten() is G
            True
        """
        return self

    def _rescale_for_frame_aspect_ratio_and_zoom(self, b, frame_aspect_ratio, zoom):
        if frame_aspect_ratio is None:
            return (b*zoom,b*zoom,b*zoom), (-b*zoom,-b*zoom,-b*zoom)
        box = [b*w for w in frame_aspect_ratio]
        # Now take the maximum length in box and rescale to b.
        s = b / max(box)
        box_max = tuple([s*w*zoom for w in box])
        box_min = tuple([-w*zoom for w in box_max])
        return box_min, box_max

    def _prepare_for_jmol(self, frame, axes, frame_aspect_ratio, aspect_ratio, zoom):
        s = 3
        box_min, box_max = self._rescale_for_frame_aspect_ratio_and_zoom(s, frame_aspect_ratio, zoom)
        a_min, a_max = self._box_for_aspect_ratio(aspect_ratio, box_min, box_max)
        return self._transform_to_bounding_box(box_min, box_max, a_min, a_max, frame=frame,
                                               axes=axes, thickness=1,
                                               labels = True)   # jmol labels are implemented

    def _prepare_for_tachyon(self, frame, axes, frame_aspect_ratio, aspect_ratio, zoom):
        box_min, box_max = self._rescale_for_frame_aspect_ratio_and_zoom(1.0, frame_aspect_ratio, zoom)
        a_min, a_max = self._box_for_aspect_ratio(aspect_ratio, box_min, box_max)
        return self._transform_to_bounding_box(box_min, box_max, a_min, a_max,
                                               frame=frame, axes=axes, thickness=.75,
                                               labels = False)  # no tachyon text implemented yet

    def _box_for_aspect_ratio(self, aspect_ratio, box_min, box_max):
        # 1. Find a box around self so that when self gets rescaled into the
        # box defined by box_min, box_max, it has the right aspect ratio
        a_min, a_max = self._safe_bounding_box()

        if aspect_ratio == "automatic" or aspect_ratio == [1.0]*3:
            return a_min, a_max

        longest_side = 0
        longest_length = a_max[0] - a_min[0]
        shortest_side = 0
        shortest_length = a_max[0] - a_min[0]

        for i in range(3):
            s = a_max[i] - a_min[i]
            if s > longest_length:
                longest_length = s
                longest_side = i
            if s < shortest_length:
                shortest_length = s
                shortest_side = i

        # 2. Rescale aspect_ratio so the shortest side is 1.
        r = float(aspect_ratio[shortest_side])
        aspect_ratio = [a/r for a in aspect_ratio]

        # 3. Extend the bounding box of self by rescaling so the sides
        # have the same ratio as aspect_ratio, and without changing
        # the longest side.
        long_box_side = box_max[longest_side] - box_min[longest_side]
        sc = [1.0,1.0,1.0]
        for i in range(3):
            # compute the length we want:
            new_length = longest_length / aspect_ratio[i]
            # change the side length by a_min and a_max so
            # that a_max[i] - a_min[i] = new_length

            # We have to take into account the ratio of the
            # sides after transforming to the bounding box.
            z = long_box_side / (box_max[i] - box_min[i])
            w = new_length / ((a_max[i] - a_min[i]) * z)
            sc[i] = w

        w = min(sc)
        sc = [z/w for z in sc]
        for i in range(3):
            a_min[i] *= sc[i]
            a_max[i] *= sc[i]

        return a_min, a_max

    def _transform_to_bounding_box(self, xyz_min, xyz_max, a_min, a_max, frame, axes, thickness, labels):

        a_min_orig = a_min
        a_max_orig = a_max

        # Rescale in each direction
        scale = [float(xyz_max[i] - xyz_min[i]) / (a_max[i] - a_min[i]) for i in range(3)]
        X = self.scale(scale)
        a_min = [scale[i]*a_min[i] for i in range(3)]
        a_max = [scale[i]*a_max[i] for i in range(3)]

        # Translate so lower left corner of original bounding box
        # is in the right spot
        T = [xyz_min[i] - a_min[i] for i in range(3)]
        X = X.translate(T)
        if frame:
            from sage.plot.plot3d.shapes2 import frame3d, frame_labels
            F = frame3d(xyz_min, xyz_max, opacity=0.5, color=(0,0,0), thickness=thickness)
            if labels:
                F += frame_labels(xyz_min, xyz_max, a_min_orig, a_max_orig)

            X += F

        if axes:
            # draw axes
            from sage.plot.plot3d.shapes import arrow3d
            A = (arrow3d((min(0,a_min[0]),0, 0), (max(0,a_max[0]), 0,0),
                             thickness, color='blue'),
                 arrow3d((0,min(0,a_min[1]), 0), (0, max(0,a_max[1]), 0),
                             thickness, color='blue'),
                 arrow3d((0, 0, min(0,a_min[2])), (0, 0, max(0,a_max[2])),
                             thickness, color='blue'))
            X += sum(A).translate([-z for z in T])

        return X

    def _process_viewing_options(self, kwds):
        """
        Process viewing options (the keywords passed to show()) and return a new
        dictionary. Defaults will be filled in for missing options and taken from
        self._extra_kwds as well. Options that have the value "automatic" will be
        automatically determined. Finally, the provided dictionary is modified
        to remove all of the keys that were used -- so that the unused keywords
        can be used elsewhere.
        """
        opts = {}
        opts.update(SHOW_DEFAULTS)
        opts.update(self._extra_kwds)
        opts.update(kwds)

        # Remove all of the keys that are viewing options, since the remaining
        # kwds might be passed on.
        for key_to_remove in SHOW_DEFAULTS:
            kwds.pop(key_to_remove, None)

        # deal with any aspect_ratio instances passed from the default options to plot
        if opts['aspect_ratio'] == 'auto':
            opts['aspect_ratio'] = 'automatic'
        if opts['aspect_ratio'] != 'automatic':
            # We need this round about way to make sure that we do not
            # store the aspect ratio that was passed on to show() by the
            # user. We let the .aspect_ratio() method take care of the
            # validity of the arguments that was passed on to show()
            original_aspect_ratio = self.aspect_ratio()
            self.aspect_ratio(opts['aspect_ratio'])
            opts['aspect_ratio'] = self.aspect_ratio()
            self.aspect_ratio(original_aspect_ratio)

        if opts['frame_aspect_ratio'] == 'automatic':
            if opts['aspect_ratio'] != 'automatic':
                # Set the aspect_ratio of the frame to be the same as that
                # of the object we are rendering given the aspect_ratio
                # we'll use for it.
                opts['frame_aspect_ratio'] = \
                    self._determine_frame_aspect_ratio(opts['aspect_ratio'])
            else:
                opts['frame_aspect_ratio'] = self.frame_aspect_ratio()
        else:
            # We need this round about way to make sure that we do not
            # store the frame aspect ratio that was passed on to show() by
            # the user. We let the .frame_aspect_ratio() method take care
            # of the validity of the arguments that was passed on to show()
            original_aspect_ratio = self.frame_aspect_ratio()
            self.frame_aspect_ratio(opts['frame_aspect_ratio'])
            opts['frame_aspect_ratio'] = self.frame_aspect_ratio()
            self.frame_aspect_ratio(original_aspect_ratio)

        if opts['aspect_ratio'] == 'automatic':
            opts['aspect_ratio'] = self.aspect_ratio()

        if not isinstance(opts['figsize'], (list,tuple)):
            opts['figsize'] = [opts['figsize'], opts['figsize']]

        return opts

    def show(self, **kwds):
        r"""
        Display graphics immediately.

        This method attempts to display the graphics immediately,
        without waiting for the currently running code (if any) to
        return to the command line. Be careful, calling it from within
        a loop will potentially launch a large number of external
        viewer programs.

        INPUT:

        - ``viewer`` -- string (default: ``'threejs'``); how to view the plot.
          Admissible values are

           * ``'threejs'``: interactive web-based 3D viewer using JavaScript
             and a WebGL renderer

           * ``'jmol'``: interactive 3D viewer using Java

           * ``'tachyon'``: ray tracer generating a static PNG image;
             can produce high-resolution graphics, but does
             not show any text labels

           * ``'canvas3d'``: web-based 3D viewer using JavaScript
             and a canvas renderer (Sage notebook only)

        - ``verbosity`` -- display information about rendering
          the figure

        - ``figsize`` -- (default: 5) x or pair [x,y] for
          numbers, e.g., [5,5]; controls the size of the output figure.
          With 'tachyon', the resolution (in number of pixels) is 100 times
          ``figsize``. This is ignored for the jmol embedded renderer.

        - ``aspect_ratio`` -- (default: ``'automatic'``) aspect
          ratio of the coordinate system itself; give [1,1,1] or 1 to make
          spheres look round

        - ``frame_aspect_ratio`` -- (default: ``'automatic'``)
          aspect ratio of frame that contains the 3d scene

        - ``zoom`` -- (default: 1) how zoomed in

        - ``frame`` -- boolean (default: ``True``); if ``True``, draw a
          bounding frame with labels

        - ``axes`` -- boolean (default: ``False``); if ``True``, draw coordinate
          axes

        - ``camera_position`` -- (for tachyon) (default: (2.3, 2.4, 2.0))
          the viewpoint, with respect to the cube
          $[-1,1]\\times[-1,1]\\times[-1,1]$,
          into which the bounding box of the scene
          is scaled and centered.
          The default viewing direction is towards the origin.

        - ``viewdir`` -- (for tachyon) (default: ``None``) three coordinates
          specifying the viewing direction

        - ``updir`` -- (for tachyon) (default: (0,0,1)) the "upward"
          direction of the camera

        - ``light_position`` -- (for tachyon) (default: (4,3,2)) the position
          of the single light source in the scene (in addition to ambient light)

        - ``antialiasing`` -- (for tachyon) (default: ``False``)

        - ``raydepth`` -- (for tachyon) (default: 8)
          see the :class:`sage.plot.plot3d.tachyon.Tachyon` class

        - ``shade`` -- (for tachyon) string (default: ``'full'``);
          shading options. Admissible values are

           * ``'full'``: best quality rendering (and slowest).
             Sets tachyon command line flag ``-fullshade``.

           * ``'medium``: good quality rendering, but no shadows.
             Sets tachyon command line flag ``-mediumshade``.

           * ``'low'``: low quality rendering, preview (and fast).
             Sets tachyon command line flag ``-lowshade``.

           * ``'lowest'``: worst quality rendering, preview (and fastest).
             Sets tachyon command line flag ``-lowestshade``.

        - ``extra_opts`` -- (for tachyon) string (default: empty string);
          extra options that will be appended to the tachyon command line

        - ``**kwds`` -- other options, which make sense for particular
          rendering engines

        OUTPUT:

        This method does not return anything. Use :meth:`save` if you
        want to save the figure as an image file.

        .. WARNING::

            By default, the jmol and tachyon viewers perform
            some non-uniform scaling of the axes.

        If this is not desired, one can set ``aspect_ratio=1``::

            sage: p = plot3d(lambda u,v:(cos(u)-cos(v)), (-0.2,0.2),(-0.2,0.2))
            sage: p.show(viewer='threejs')
            sage: p.show(viewer='jmol')
            sage: p.show(viewer='jmol',aspect_ratio=1)
            sage: p.show(viewer='tachyon',camera_position=(4,0,0))
            sage: p.show(viewer='tachyon',camera_position=(2,2,0.3),aspect_ratio=1)

        CHANGING DEFAULTS: Defaults can be uniformly changed by importing a
        dictionary and changing it. For example, here we change the default
        so images display without a frame instead of with one::

            sage: from sage.plot.plot3d.base import SHOW_DEFAULTS
            sage: SHOW_DEFAULTS['frame'] = False

        This sphere will not have a frame around it::

            sage: sphere((0,0,0))
            Graphics3d Object

        We change the default back::

            sage: SHOW_DEFAULTS['frame'] = True

        Now this sphere is enclosed in a frame::

            sage: sphere((0,0,0))
            Graphics3d Object

        EXAMPLES: We illustrate use of the ``aspect_ratio`` option::

            sage: x, y = var('x,y')                                                     # needs sage.symbolic
            sage: p = plot3d(2*sin(x*y), (x, -pi, pi), (y, -pi, pi))                    # needs sage.symbolic
            sage: p.show(aspect_ratio=[1,1,1])                                          # needs sage.symbolic

        This looks flattened, but filled with the plot::

            sage: p.show(frame_aspect_ratio=[1,1,1/16])                                 # needs sage.symbolic

        This looks flattened, but the plot is square and smaller::

            sage: p.show(aspect_ratio=[1,1,1], frame_aspect_ratio=[1,1,1/8])            # needs sage.symbolic

        This example shows indirectly that the defaults
        from :func:`~sage.plot.plot.plot` are dealt with properly::

            sage: plot(vector([1,2,3]))
            Graphics3d Object

        We use the 'canvas3d' backend from inside the notebook to get a view of
        the plot rendered inline using HTML canvas::

            sage: p.show(viewer='canvas3d')                                             # needs sage.symbolic

        Sometimes shadows in Tachyon-produced images can lead to confusing
        plots. To remove them::

            sage: p.show(viewer='tachyon', shade='medium')                              # needs sage.symbolic

        One can also pass Tachyon command line flags directly. For example,
        the following line produces the same result as the previous one::

            sage: p.show(viewer='tachyon', extra_opts='-mediumshade')                   # needs sage.symbolic
        """
        from sage.repl.rich_output import get_display_manager
        dm = get_display_manager()
        dm.display_immediately(self, **kwds)

    def _save_image_png(self, filename, **kwds):
        r"""
        Save a PNG rendering.

        This private method is only for use by :meth:`save_image`.

        EXAMPLES::

            sage: s = sphere()
            sage: png = tmp_filename(ext='.png')
            sage: s._save_image_png(png)
            sage: with open(png, 'rb') as fobj:
            ....:     fobj.read().startswith(b'\x89PNG')
            True

            sage: s._save_image_png('/path/to/foo.bar')
            Traceback (most recent call last):
            ...
            AssertionError
        """
        assert filename.endswith('.png')
        opts = self._process_viewing_options(kwds)
        viewer = opts['viewer']
        if viewer == 'threejs':
            viewer = 'jmol'  # since threejs has no png dump
        if viewer == 'tachyon':
            from sage.repl.rich_output.output_catalog import OutputImagePng
            render = self._rich_repr_tachyon(OutputImagePng, **opts)
            render.png.save_as(filename)
        elif viewer == 'jmol':
            scene = self._rich_repr_jmol(**opts)
            scene.preview_png.save_as(filename)
        else:
            raise ValueError('cannot use viewer={0} to render image'.format(viewer))

    def save_image(self, filename, **kwds):
        r"""
        Save a 2-D image rendering.

        The image type is determined by the extension of the filename.
        For example, this could be ``.png``, ``.jpg``, ``.gif``,
        ``.pdf``, ``.svg``.

        INPUT:

        - ``filename`` -- string; the file name under which to save
          the image

        Any further keyword arguments are passed to the renderer.

        EXAMPLES::

            sage: G = sphere()
            sage: png = tmp_filename(ext='.png')
            sage: G.save_image(png)
            sage: with open(png, 'rb') as fobj:
            ....:     assert fobj.read().startswith(b'\x89PNG')

            sage: gif = tmp_filename(ext='.gif')
            sage: G.save_image(gif)
            sage: with open(gif, 'rb') as fobj:
            ....:     assert fobj.read().startswith(b'GIF')
        """
        ext = os.path.splitext(filename)[1].lower()
        if ext not in ['.bmp', '.png', '.gif', '.ppm', '.tiff', '.tif',
                       '.jpg', '.jpeg']:
            raise ValueError('unknown image file type: {0}'.format(ext))
        if ext == '.png':
            self._save_image_png(filename, **kwds)
        else:
            png = tmp_filename(ext='.png')
            self._save_image_png(png, **kwds)
            import PIL.Image as Image
            Image.open(png).save(filename)

    def save(self, filename, **kwds):
        """
        Save the graphic in a file.

        The file type depends on the file extension you give in the
        filename. This can be either:

        - an image file (of type: PNG, BMP, GIF, PPM, or TIFF) rendered
          using Jmol (default) or Tachyon,

        - a Sage object file (of type ``.sobj``) that you can load back later
          (a pickle),

        - an HTML file depicting the graphic using the Three.js viewer,

        - a data file (of type: X3D, STL, AMF, PLY) for export and use in
          other software.

        For data files, the support is only partial. For instance STL and
        AMF only works for triangulated surfaces. The prefered format is X3D.

        INPUT:

        - ``filename`` -- string; where to save the image or object

        - ``**kwds`` -- when specifying an image file to be rendered by Tachyon
          or Jmol, any of the viewing options accepted by :meth:`show` are valid as
          keyword arguments to this function and they will behave in the same
          way. Accepted keywords include: ``viewer``, ``verbosity``,
          ``figsize``, ``aspect_ratio``, ``frame_aspect_ratio``, ``zoom``,
          ``frame``, and ``axes``. Default values are provided.

        EXAMPLES::

            sage: f = tmp_filename(ext='.png')
            sage: G = sphere()
            sage: G.save(f)

        We demonstrate using keyword arguments to control the appearance of the
        output image::

            sage: G.save(f, zoom=2, figsize=[5, 10])

        Using Tachyon instead of the default viewer (Jmol) to create the
        image::

            sage: G.save(f, viewer='tachyon')

        Since Tachyon only outputs PNG images, PIL will be used to convert to
        alternate formats::

            sage: cube().save(tmp_filename(ext='.gif'), viewer='tachyon')

        Here is how to save in one of the data formats::

            sage: f = tmp_filename(ext='.x3d')
            sage: cube().save(f)

            sage: open(f).read().splitlines()[7]
            "<Shape><Box size='0.5 0.5 0.5'/><Appearance><Material diffuseColor='0.4 0.4 1.0' shininess='1.0' specularColor='0.0 0.0 0.0'/></Appearance></Shape>"

        Producing a Three.js-based HTML file::

            sage: f = tmp_filename(ext='.html')
            sage: G.save(f, frame=False, online=True)
        """
        ext = os.path.splitext(filename)[1].lower()
        if ext == '' or ext == '.sobj':
            SageObject.save(self, filename)
        elif ext in ['.bmp', '.png', '.gif', '.ppm', '.tiff', '.tif',
                     '.jpg', '.jpeg']:
            self.save_image(filename, **kwds)
        elif filename.endswith('.spt.zip'):
            scene = self._rich_repr_jmol(**kwds)
            scene.jmol.save(filename)
        elif ext == '.x3d':
            with open(filename, 'w') as outfile:
                outfile.write(self.x3d())
        elif ext == '.stl':
            with open(filename, 'wb') as outfile:
                outfile.write(self.stl_binary())
        elif ext == '.amf':
            # todo : zip the output file ?
            with open(filename, 'w') as outfile:
                outfile.write(self.amf_ascii_string())
        elif ext == '.ply':
            with open(filename, 'w') as outfile:
                outfile.write(self.ply_ascii_string())
        elif ext == '.html':
            self._rich_repr_threejs(**kwds).html.save_as(filename)
        else:
            raise ValueError('filetype {} not supported by save()'.format(ext))

    def stl_binary(self):
        """
        Return an STL (STereoLithography) binary representation of the surface.

        .. WARNING::

            This only works for surfaces, transforms and unions of surfaces,
            but not for general plot objects!

        OUTPUT: a binary string that represents the surface in the binary STL format

        See :wikipedia:`STL_(file_format)`

        .. SEEALSO:: :meth:`stl_ascii_string`

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: x,y,z = var('x,y,z')
            sage: a = implicit_plot3d(x^2+y^2+z^2-9,[x,-5,5],[y,-5,5],[z,-5,5])
            sage: astl = a.stl_binary()
            sage: print(astl[:40].decode('ascii'))
            STL binary file / made by SageMath / ###

            sage: p = polygon3d([[0,0,0], [1,2,3], [3,0,0]])
            sage: print(p.stl_binary()[:40].decode('ascii'))
            STL binary file / made by SageMath / ###

        This works when faces have more then 3 sides::

            sage: # needs sage.geometry.polyhedron sage.groups
            sage: P = polytopes.dodecahedron()
            sage: Q = P.plot().all[-1]
            sage: print(Q.stl_binary()[:40].decode('ascii'))
            STL binary file / made by SageMath / ###
        """
        import struct
        header = b'STL binary file / made by SageMath / '
        header += b'#' * (80 - len(header))
        # header = 80 bytes, arbitrary ascii characters
        data = self.stl_binary_repr(self.default_render_params())
        N_triangles = len(data)
        return b''.join([header, struct.pack('I', N_triangles)] + data)

    def stl_ascii_string(self, name='surface'):
        """
        Return an STL (STereoLithography) representation of the surface.

        .. WARNING::

            This only works for surfaces, not for general plot objects!

        INPUT:

        - ``name`` -- string (default: ``'surface'``); name of the surface

        OUTPUT: string that represents the surface in the STL format

        See :wikipedia:`STL_(file_format)`

        .. SEEALSO:: :meth:`stl_binary`

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: x,y,z = var('x,y,z')
            sage: a = implicit_plot3d(x^2+y^2+z^2-9,[x,-5,5],[y,-5,5],[z,-5,5])
            sage: astl = a.stl_ascii_string()
            sage: astl.splitlines()[:7]  # abs tol 1e-10
            ['solid surface',
            'facet normal 0.9733285267845754 -0.16222142113076257 -0.16222142113076257',
            '    outer loop',
            '        vertex 2.94871794872 -0.384615384615 -0.39358974359',
            '        vertex 2.95021367521 -0.384615384615 -0.384615384615',
            '        vertex 2.94871794872 -0.39358974359 -0.384615384615',
            '    endloop']

            sage: p = polygon3d([[0,0,0], [1,2,3], [3,0,0]])
            sage: print(p.stl_ascii_string(name='triangle'))
            solid triangle
            facet normal 0.0 0.8320502943378436 -0.5547001962252291
                outer loop
                    vertex 0.0 0.0 0.0
                    vertex 1.0 2.0 3.0
                    vertex 3.0 0.0 0.0
                endloop
            endfacet
            endsolid triangle

        Now works when faces have more then 3 sides::

            sage: # needs sage.geometry.polyhedron sage.groups
            sage: P = polytopes.dodecahedron()
            sage: Q = P.plot().all[-1]
            sage: print(Q.stl_ascii_string().splitlines()[:7])
            ['solid surface',
             'facet normal 0.0 0.5257311121191338 0.8506508083520399',
             '    outer loop',
             '        vertex -0.7639320225002102 0.7639320225002102 0.7639320225002102',
             '        vertex -0.4721359549995796 0.0 1.2360679774997898',
             '        vertex 0.4721359549995796 0.0 1.2360679774997898',
             '    endloop']
        """
        from sage.modules.free_module import FreeModule
        RR3 = FreeModule(RDF, 3)

        faces = self.face_list()
        if not faces:
            self.triangulate()
            faces = self.face_list()

        code = ("facet normal {!r} {!r} {!r}\n"
                "    outer loop\n"
                "        vertex {!r} {!r} {!r}\n"
                "        vertex {!r} {!r} {!r}\n"
                "        vertex {!r} {!r} {!r}\n"
                "    endloop\n"
                "endfacet\n")

        faces_iter = faces.__iter__()

        def chopped_faces_iter():
            for face in faces_iter:
                n = len(face)
                if n == 3:
                    yield face
                else:
                    # naive cut into triangles
                    v = face[-1]
                    for i in range(n - 2):
                        yield [v, face[i], face[i + 1]]

        string_list = ["solid {}\n".format(name)]
        for i, j, k in chopped_faces_iter():
            ij = RR3(j) - RR3(i)
            ik = RR3(k) - RR3(i)
            n = ij.cross_product(ik)
            n = n / n.norm()
            string_list += [code.format(n[0], n[1], n[2],
                                        i[0], i[1], i[2],
                                        j[0], j[1], j[2],
                                        k[0], k[1], k[2])]
        string_list += ["endsolid {}".format(name)]
        return "".join(string_list)

    def ply_ascii_string(self, name='surface'):
        """
        Return a PLY (Polygon File Format) representation of the surface.

        INPUT:

        - ``name`` -- string (default: ``'surface'``); name of the surface

        OUTPUT: string that represents the surface in the PLY format

        See :wikipedia:`PLY_(file_format)`.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: x,y,z = var('x,y,z')
            sage: a = implicit_plot3d(x^2+y^2+z^2-9,[x,-5,5],[y,-5,5],[z,-5,5])
            sage: astl = a.ply_ascii_string()
            sage: astl.splitlines()[:10]
            ['ply',
            'format ascii 1.0',
            'comment surface',
            'element vertex 15540',
            'property float x',
            'property float y',
            'property float z',
            'element face 5180',
            'property list uchar int vertex_indices',
            'end_header']

            sage: p = polygon3d([[0,0,0], [1,2,3], [3,0,0]])
            sage: print(p.ply_ascii_string(name='triangle'))
            ply
            format ascii 1.0
            comment triangle
            element vertex 3
            property float x
            property float y
            property float z
            element face 1
            property list uchar int vertex_indices
            end_header
            0.0 0.0 0.0
            1.0 2.0 3.0
            3.0 0.0 0.0
            3 0 1 2
        """
        faces = self.index_faces()
        if not faces:
            self.triangulate()
            faces = self.index_faces()

        string_list = ["ply\nformat ascii 1.0\ncomment {}\nelement vertex {}\nproperty float x\nproperty float y\nproperty float z\nelement face {}\nproperty list uchar int vertex_indices\nend_header\n".format(name, len(self.vertex_list()), len(faces))]

        vertex_template = '{} {} {}\n'
        for v in self.vertices():
            string_list += [vertex_template.format(*v)]

        for f in faces:
            string_list += [str(len(f))
                            + ''.join(' {}'.format(k) for k in f) + '\n']

        return "".join(string_list)

    def amf_ascii_string(self, name='surface'):
        """
        Return an AMF (Additive Manufacturing File Format) representation of
        the surface.

        .. WARNING::

            This only works for triangulated surfaces!

        INPUT:

        - ``name`` -- string (default: ``'surface'``); name of the surface

        OUTPUT: string that represents the surface in the AMF format

        See :wikipedia:`Additive_Manufacturing_File_Format`.

        .. TODO::

            This should rather be saved as a ZIP archive to save space.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: x,y,z = var('x,y,z')
            sage: a = implicit_plot3d(x^2+y^2+z^2-9,[x,-5,5],[y,-5,5],[z,-5,5])
            sage: a_amf = a.amf_ascii_string()
            sage: a_amf[:160]
            '<?xml version="1.0" encoding="utf-8"?><amf><object id="surface"><mesh><vertices><vertex><coordinates><x>2.948717948717948</x><y>-0.384615384615385</y><z>-0.3935'

            sage: p = polygon3d([[0,0,0], [1,2,3], [3,0,0]])
            sage: print(p.amf_ascii_string(name='triangle'))
            <?xml version="1.0" encoding="utf-8"?><amf><object id="triangle"><mesh><vertices><vertex><coordinates><x>0.0</x><y>0.0</y><z>0.0</z></coordinates></vertex><vertex><coordinates><x>1.0</x><y>2.0</y><z>3.0</z></coordinates></vertex><vertex><coordinates><x>3.0</x><y>0.0</y><z>0.0</z></coordinates></vertex></vertices><volume><triangle><v1>0</v1><v2>1</v2><v3>2</v3></triangle></volume></mesh></object></amf>
        """
        faces = self.index_faces()
        if not faces:
            self.triangulate()
            faces = self.index_faces()

        if len(faces[0]) > 3:
            raise ValueError('not made of triangles')

        string_list = ['<?xml version="1.0" encoding="utf-8"?><amf><object id="{}"><mesh>'.format(name)]

        string_list += ['<vertices>']
        vertex_template = '<vertex><coordinates><x>{!r}</x><y>{!r}</y><z>{!r}</z></coordinates></vertex>'
        for v in self.vertices():
            string_list += [vertex_template.format(*v)]
        string_list += ['</vertices><volume>']

        face_template = '<triangle><v1>{}</v1><v2>{}</v2><v3>{}</v3></triangle>'
        for i, j, k in faces:
            string_list += face_template.format(i, j, k)

        string_list += ['</volume></mesh></object></amf>']
        return "".join(string_list)

    def plot(self):
        """
        Draw a 3D plot of this graphics object, which just returns this
        object since this is already a 3D graphics object.
        Needed to support PLOT in docstrings, see :issue:`17498`

        EXAMPLES::

            sage: S = sphere((0,0,0), 2)
            sage: S.plot() is S
            True
        """
        return self

# if you add any default parameters you must update some code below
SHOW_DEFAULTS = {'viewer': 'threejs',
                 'verbosity': 0,
                 'figsize': 5,
                 'aspect_ratio': "automatic",
                 'frame_aspect_ratio': "automatic",
                 'zoom': 1,
                 'frame': True,
                 'axes': False}


class Graphics3dGroup(Graphics3d):
    """
    This class represents a collection of 3d objects. Usually they are formed
    implicitly by summing.
    """
    def __init__(self, all=(), rot=None, trans=None, scale=None, T=None):
        """
        EXAMPLES::

            sage: sage.plot.plot3d.base.Graphics3dGroup([icosahedron(), dodecahedron(opacity=.5)])
            Graphics3d Object
            sage: type(icosahedron() + dodecahedron(opacity=.5))
            <class 'sage.plot.plot3d.base.Graphics3dGroup'>
        """
        self.all = list(all)
        self.frame_aspect_ratio(optimal_aspect_ratios([a.frame_aspect_ratio() for a in all]))
        self.aspect_ratio(optimal_aspect_ratios([a.aspect_ratio() for a in all]))
        self._set_extra_kwds(optimal_extra_kwds([a._extra_kwds for a in all]))

    def __add__(self, other):
        """
        We override this here to make large sums more efficient.

        EXAMPLES::

            sage: G = sum(tetrahedron(opacity=1-t/11).translate(t, 0, 0) for t in range(10))
            sage: G
            Graphics3d Object
            sage: len(G.all)
            10

        We check that :issue:`17258` is solved::

            sage: g = point3d([0,-2,-2]); g += point3d([2,-2,-2])
            sage: len(g.all)
            2
            sage: h = g + arrow([0,-2,-2], [2,-2,-2])
            sage: len(g.all)
            2
            sage: g == h
            False
        """
        if type(self) is Graphics3dGroup and isinstance(other, Graphics3d):
            s_all = list(self.all)
            s_all.append(other)
            return Graphics3dGroup(s_all)
        else:
            return Graphics3d.__add__(self, other)

    def bounding_box(self):
        """
        Box that contains the bounding boxes of
        the objects.

        EXAMPLES::

            sage: A = sphere((0,0,0), 5)
            sage: B = sphere((1, 5, 10), 1)
            sage: A.bounding_box()
            ((-5.0, -5.0, -5.0), (5.0, 5.0, 5.0))
            sage: B.bounding_box()
            ((0.0, 4.0, 9.0), (2.0, 6.0, 11.0))
            sage: (A+B).bounding_box()
            ((-5.0, -5.0, -5.0), (5.0, 6.0, 11.0))
            sage: (A+B).show(aspect_ratio=1, frame=True)

            sage: sage.plot.plot3d.base.Graphics3dGroup([]).bounding_box()
            ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
        """
        if len(self.all) == 0:
            return Graphics3d.bounding_box(self)
        v = [obj.bounding_box() for obj in self.all]
        return min3([a[0] for a in v]), max3([a[1] for a in v])

    def transform(self, **kwds):
        """
        Transforming this entire group simply makes a transform group with
        the same contents.

        EXAMPLES::

            sage: G = dodecahedron(color='red', opacity=.5) + icosahedron(color='blue')
            sage: G
            Graphics3d Object
            sage: G.transform(scale=(2,1/2,1))
            Graphics3d Object
            sage: G.transform(trans=(1,1,3))
            Graphics3d Object
        """
        T = TransformGroup(self.all, **kwds)
        T._set_extra_kwds(self._extra_kwds)
        return T

    def set_texture(self, **kwds):
        """
        EXAMPLES::

            sage: G = dodecahedron(color='red', opacity=.5) + icosahedron((3, 0, 0), color='blue')
            sage: G
            Graphics3d Object
            sage: G.set_texture(color='yellow')
            sage: G
            Graphics3d Object
        """
        for g in self.all:
            g.set_texture(**kwds)

    def json_repr(self, render_params):
        """
        The JSON representation of a group is simply the concatenation of the
        representations of its objects.

        EXAMPLES::

            sage: G = sphere() + sphere((1, 2, 3))
            sage: G.json_repr(G.default_render_params())
            [[['{"vertices":...']], [['{"vertices":...']]]
        """
        return [g.json_repr(render_params) for g in self.all]

    def tachyon_repr(self, render_params):
        """
        The tachyon representation of a group is simply the concatenation of
        the representations of its objects.

        EXAMPLES::

            sage: G = sphere() + sphere((1,2,3))
            sage: G.tachyon_repr(G.default_render_params())
            [['Sphere center 0.0 0.0 0.0 Rad 1.0 texture...'],
             ['Sphere center 1.0 2.0 3.0 Rad 1.0 texture...']]
        """
        return [g.tachyon_repr(render_params) for g in self.all]

    def x3d_str(self):
        """
        The x3d representation of a group is simply the concatenation of
        the representation of its objects.

        EXAMPLES::

            sage: G = sphere() + sphere((1,2,3))
            sage: print(G.x3d_str())
            <Transform translation='0 0 0'>
            <Shape><Sphere radius='1.0'/><Appearance><Material diffuseColor='0.4 0.4 1.0' shininess='1.0' specularColor='0.0 0.0 0.0'/></Appearance></Shape>
            </Transform>
            <Transform translation='1 2 3'>
            <Shape><Sphere radius='1.0'/><Appearance><Material diffuseColor='0.4 0.4 1.0' shininess='1.0' specularColor='0.0 0.0 0.0'/></Appearance></Shape>
            </Transform>
        """
        return "\n".join(g.x3d_str() for g in self.all)

    def obj_repr(self, render_params):
        """
        The obj representation of a group is simply the concatenation of
        the representation of its objects.

        EXAMPLES::

            sage: G = tetrahedron() + tetrahedron().translate(10, 10, 10)
            sage: G.obj_repr(G.default_render_params())
            [['g obj_1',
              'usemtl ...',
              ['v 0 0 1',
               'v 0.942809 0 -0.333333',
               'v -0.471405 0.816497 -0.333333',
               'v -0.471405 -0.816497 -0.333333'],
              ['f 1 2 3', 'f 2 4 3', 'f 1 3 4', 'f 1 4 2'],
              []],
             [['g obj_2',
               'usemtl ...',
               ['v 10 10 11',
                'v 10.9428 10 9.66667',
                'v 9.5286 10.8165 9.66667',
                'v 9.5286 9.1835 9.66667'],
               ['f 5 6 7', 'f 6 8 7', 'f 5 7 8', 'f 5 8 6'],
               []]]]
        """
        return [g.obj_repr(render_params) for g in self.all]

    def jmol_repr(self, render_params):
        r"""
        The jmol representation of a group is simply the concatenation of
        the representation of its objects.

        EXAMPLES::

            sage: G = sphere() + sphere((1,2,3))
            sage: G.jmol_repr(G.default_render_params())
            [[['isosurface sphere_1  center {0.0 0.0 0.0} sphere 1.0\ncolor isosurface  [102,102,255]']],
             [['isosurface sphere_2  center {1.0 2.0 3.0} sphere 1.0\ncolor isosurface  [102,102,255]']]]
        """
        return [g.jmol_repr(render_params) for g in self.all]

    def threejs_repr(self, render_params):
        r"""
        The three.js representation of a group is the concatenation of the
        representations of its objects.

        EXAMPLES::

            sage: G = point3d((1,2,3)) + point3d((4,5,6)) + line3d([(1,2,3), (4,5,6)])
            sage: G.threejs_repr(G.default_render_params())
            [('point',
              {'color': '#6666ff', 'opacity': 1.0, 'point': (1.0, 2.0, 3.0), 'size': 5.0}),
             ('point',
              {'color': '#6666ff', 'opacity': 1.0, 'point': (4.0, 5.0, 6.0), 'size': 5.0}),
             ('line',
              {'color': '#6666ff',
               'linewidth': 1.0,
               'opacity': 1.0,
               'points': [(1.0, 2.0, 3.0), (4.0, 5.0, 6.0)]})]
        """
        reprs = []
        for g in self.all:
            reprs += g.threejs_repr(render_params)
        return reprs

    def stl_binary_repr(self, render_params):
        r"""
        The stl binary representation of a group is simply the
        concatenation of the representation of its objects.

        The STL binary representation is a list of binary strings,
        one for each triangle.

        EXAMPLES::

            sage: G = sphere() + sphere((1,2,3))
            sage: len(G.stl_binary_repr(G.default_render_params()))
            2736
        """
        data = []
        for sub in self.all:
            data.extend(sub.stl_binary_repr(render_params))
        return data

    def texture_set(self):
        """
        The texture set of a group is simply the union of the textures of
        all its objects.

        EXAMPLES::

            sage: G = sphere(color='red') + sphere(color='yellow')
            sage: [t for t in G.texture_set() if t.color == colors.red] # one red texture
            [Texture(texture..., red, ff0000)]
            sage: [t for t in G.texture_set() if t.color == colors.yellow] # one yellow texture
            [Texture(texture..., yellow, ffff00)]

            sage: T = sage.plot.plot3d.texture.Texture('blue'); T
            Texture(texture..., blue, 0000ff)
            sage: G = sphere(texture=T) + sphere((1, 1, 1), texture=T)
            sage: len(G.texture_set())
            1

        TESTS:

        Check that :issue:`23200` is fixed::

            sage: G = sage.plot.plot3d.base.Graphics3dGroup()
            sage: G.texture_set()
            set()
        """
        if not self.all:
            return set()
        return reduce(set.union, [g.texture_set() for g in self.all])

    def flatten(self):
        """
        Try to reduce the depth of the scene tree by consolidating groups
        and transformations.

        EXAMPLES::

            sage: G = sum([circle((0, 0), t) for t in [1..10]], sphere()); G
            Graphics3d Object
            sage: G.flatten()
            Graphics3d Object
            sage: len(G.all)
            2
            sage: len(G.flatten().all)
            11
        """
        if len(self.all) == 1:
            return self.all[0].flatten()
        all = []
        for g in self.all:
            g = g.flatten()
            if type(g) is Graphics3dGroup:
                all += g.all
            else:
                all.append(g)
        return Graphics3dGroup(all)

    def plot(self):
        return self


class TransformGroup(Graphics3dGroup):
    """
    This class is a container for a group of objects with a common
    transformation.
    """
    def __init__(self, all=[], rot=None, trans=None, scale=None, T=None):
        """
        EXAMPLES::

            sage: sage.plot.plot3d.base.TransformGroup([sphere()], trans=(1,2,3)) + point3d((0,0,0))
            Graphics3d Object

        The are usually constructed implicitly::

            sage: type(sphere((1,2,3)))
            <class 'sage.plot.plot3d.base.TransformGroup'>
            sage: type(dodecahedron().scale(2))
            <class 'sage.plot.plot3d.base.TransformGroup'>
        """
        Graphics3dGroup.__init__(self, all)
        self._rot = rot
        self._trans = trans
        if scale is not None and len(scale) == 1:
            if isinstance(scale, (tuple, list)):
                scale = scale[0]
            scale = (scale, scale, scale)
        self._scale = scale
        if T is not None:
            self.T = T
        self.frame_aspect_ratio(optimal_aspect_ratios([a.frame_aspect_ratio() for a in all]))
        self.aspect_ratio(optimal_aspect_ratios([a.aspect_ratio() for a in all]))
        self._set_extra_kwds(optimal_extra_kwds([a._extra_kwds for a in all if a._extra_kwds is not None]))

    def bounding_box(self):
        """
        Return the bounding box, i.e., the box containing the
        contents of the object after applying the transformation.

        EXAMPLES::

            sage: from math import pi
            sage: G = cube()
            sage: G.bounding_box()
            ((-0.5, -0.5, -0.5), (0.5, 0.5, 0.5))
            sage: G.scale(4).bounding_box()
            ((-2.0, -2.0, -2.0), (2.0, 2.0, 2.0))
            sage: G.rotateZ(pi/4).bounding_box()
            ((-0.7071067811865475, -0.7071067811865475, -0.5),
             (0.7071067811865475, 0.7071067811865475, 0.5))
        """
        try:
            return self._bounding_box
        except AttributeError:
            pass

        cdef Transformation T = self.get_transformation()
        w = sum([T.transform_bounding_box(obj.bounding_box()) for obj in self.all], ())
        self._bounding_box = point_list_bounding_box(w)
        return self._bounding_box

    def x3d_str(self):
        r"""
        To apply a transformation to a set of objects in x3d, simply make them
        all children of an x3d Transform node.

        EXAMPLES::

            sage: sphere((1,2,3)).x3d_str()
            "<Transform translation='1 2 3'>\n<Shape><Sphere radius='1.0'/><Appearance><Material diffuseColor='0.4 0.4 1.0' shininess='1.0' specularColor='0.0 0.0 0.0'/></Appearance></Shape>\n\n</Transform>"
        """
        s = "<Transform"
        if self._rot is not None:
            s += " rotation='%s %s %s %s'" % tuple(self._rot)
        if self._trans is not None:
            s += " translation='%s %s %s'" % tuple(self._trans)
        if self._scale is not None:
            s += " scale='%s %s %s'" % tuple(self._scale)
        s += ">\n"
        s += Graphics3dGroup.x3d_str(self)
        s += "\n</Transform>"
        return s

    def json_repr(self, render_params):
        """
        Transformations are applied at the leaf nodes.

        EXAMPLES::

            sage: G = cube().rotateX(0.2)
            sage: G.json_repr(G.default_render_params())
            [['{"vertices":[{"x":0.5,"y":0.589368,"z":0.390699},...']]
        """
        render_params.push_transform(self.get_transformation())
        rep = [g.json_repr(render_params) for g in self.all]
        render_params.pop_transform()
        return rep

    def stl_binary_repr(self, render_params):
        """
        Transformations are applied at the leaf nodes.

        The STL binary representation is a list of binary strings,
        one for each triangle.

        EXAMPLES::

            sage: G = sphere().translate((1,2,0))
            sage: len(G.stl_binary_repr(G.default_render_params()))
            1368
        """
        render_params.push_transform(self.get_transformation())
        rep = sum([g.stl_binary_repr(render_params) for g in self.all], [])
        render_params.pop_transform()
        return rep

    def tachyon_repr(self, render_params):
        """
        Transformations for Tachyon are applied at the leaf nodes.

        EXAMPLES::

            sage: G = sphere((1,2,3)).scale(2)
            sage: G.tachyon_repr(G.default_render_params())
            [['Sphere center 2.0 4.0 6.0 Rad 2.0 texture...']]
        """
        render_params.push_transform(self.get_transformation())
        rep = [g.tachyon_repr(render_params) for g in self.all]
        render_params.pop_transform()
        return rep

    def obj_repr(self, render_params):
        """
        Transformations for .obj files are applied at the leaf nodes.

        EXAMPLES::

            sage: G = cube().scale(4).translate(1, 2, 3)
            sage: G.obj_repr(G.default_render_params())
            [[['g obj_1',
               'usemtl ...',
               ['v 3 4 5',
                'v -1 4 5',
                'v -1 0 5',
                'v 3 0 5',
                'v 3 4 1',
                'v -1 4 1',
                'v 3 0 1',
                'v -1 0 1'],
               ['f 1 2 3 4',
                'f 1 5 6 2',
                'f 1 4 7 5',
                'f 6 5 7 8',
                'f 7 4 3 8',
                'f 3 2 6 8'],
               []]]]
        """
        render_params.push_transform(self.get_transformation())
        rep = [g.obj_repr(render_params) for g in self.all]
        render_params.pop_transform()
        return rep

    def jmol_repr(self, render_params):
        r"""
        Transformations for jmol are applied at the leaf nodes.

        EXAMPLES::

            sage: G = sphere((1,2,3)).scale(2)
            sage: G.jmol_repr(G.default_render_params())
            [[['isosurface sphere_1  center {2.0 4.0 6.0} sphere 2.0\ncolor isosurface  [102,102,255]']]]
        """
        render_params.push_transform(self.get_transformation())
        rep = [g.jmol_repr(render_params) for g in self.all]
        render_params.pop_transform()
        return rep

    def threejs_repr(self, render_params):
        r"""
        Transformations for three.js are applied at the leaf nodes.

        EXAMPLES::

            sage: G = point3d((1,2,3)) + point3d((4,5,6))
            sage: G = G.translate(-1, -2, -3).scale(10)
            sage: G.threejs_repr(G.default_render_params())
            [('point',
              {'color': '#6666ff', 'opacity': 1.0, 'point': (0.0, 0.0, 0.0), 'size': 5.0}),
             ('point',
              {'color': '#6666ff',
               'opacity': 1.0,
               'point': (30.0, 30.0, 30.0),
               'size': 5.0})]
        """
        render_params.push_transform(self.get_transformation())
        rep = Graphics3dGroup.threejs_repr(self, render_params)
        render_params.pop_transform()
        return rep

    def get_transformation(self):
        """
        Return the current transformation object.

        EXAMPLES::

            sage: G = sphere().scale(100)
            sage: T = G.get_transformation()
            sage: T.get_matrix()
            [100.0   0.0   0.0   0.0]
            [  0.0 100.0   0.0   0.0]
            [  0.0   0.0 100.0   0.0]
            [  0.0   0.0   0.0   1.0]
        """
        try:
            return self.T
        except AttributeError:
            self.T = Transformation(self._scale, self._rot, self._trans)
            return self.T

    def flatten(self):
        """
        Try to reduce the depth of the scene tree by consolidating groups
        and transformations.

        EXAMPLES::

            sage: G = sphere((1,2,3)).scale(100)
            sage: T = G.get_transformation()
            sage: T.get_matrix()
            [100.0   0.0   0.0   0.0]
            [  0.0 100.0   0.0   0.0]
            [  0.0   0.0 100.0   0.0]
            [  0.0   0.0   0.0   1.0]

            sage: G.flatten().get_transformation().get_matrix()
            [100.0   0.0   0.0 100.0]
            [  0.0 100.0   0.0 200.0]
            [  0.0   0.0 100.0 300.0]
            [  0.0   0.0   0.0   1.0]
        """
        G = Graphics3dGroup.flatten(self)
        if isinstance(G, TransformGroup):
            return TransformGroup(G.all, T=self.get_transformation() * G.get_transformation())
        elif isinstance(G, Graphics3dGroup):
            return TransformGroup(G.all, T=self.get_transformation())
        else:
            return TransformGroup([G], T=self.get_transformation())

    def transform(self, **kwds):
        """
        Transforming this entire group can be done by composing transformations.

        EXAMPLES::

            sage: G = dodecahedron(color='red', opacity=.5) + icosahedron(color='blue')
            sage: G
            Graphics3d Object
            sage: G.transform(scale=(2,1/2,1))
            Graphics3d Object
            sage: G.transform(trans=(1,1,3))
            Graphics3d Object
        """
        return Graphics3d.transform(self, **kwds)


class KeyframeAnimationGroup(Graphics3dGroup):
    """A group of objects, each depicting a single frame of animation"""
    def __init__(self, all=(), **kwds):
        r"""
        EXAMPLES::

            sage: frames = [dodecahedron(), icosahedron(), tetrahedron()]
            sage: sage.plot.plot3d.base.KeyframeAnimationGroup(frames)
            Graphics3d Object

        They are usually constructed from an class:`~sage.plot.animate.Animation`::

            sage: type(animate(frames).interactive())
            <class 'sage.plot.plot3d.base.KeyframeAnimationGroup'>
        """
        Graphics3dGroup.__init__(self, all)
        self._extra_kwds.update(kwds)

    def threejs_repr(self, render_params):
        r"""
        Add keyframe information to the representations of the group's contents.

        EXAMPLES::

            sage: a = point3d((0, 0, 1))
            sage: b = point3d((0, 1, 0))
            sage: c = point3d((1, 0, 0))
            sage: g = sage.plot.plot3d.base.KeyframeAnimationGroup([a, b, c])
            sage: g.threejs_repr(g.default_render_params())
            [('point', {..., 'keyframe': 0, ..., 'point': (0.0, 0.0, 1.0), ...}),
             ('point', {..., 'keyframe': 1, ..., 'point': (0.0, 1.0, 0.0), ...}),
             ('point', {..., 'keyframe': 2, ..., 'point': (1.0, 0.0, 0.0), ...})]

        Only top-level objects get a unique keyframe. Nested objects share the
        same keyframe::

            sage: g = sage.plot.plot3d.base.KeyframeAnimationGroup([a + b, c])
            sage: g.threejs_repr(g.default_render_params())
            [('point', {..., 'keyframe': 0, ..., 'point': (0.0, 0.0, 1.0), ...}),
             ('point', {..., 'keyframe': 0, ..., 'point': (0.0, 1.0, 0.0), ...}),
             ('point', {..., 'keyframe': 1, ..., 'point': (1.0, 0.0, 0.0), ...})]
        """
        reprs = []
        for i, g in enumerate(self.all):
            for (kind, desc) in g.threejs_repr(render_params):
                desc['keyframe'] = i
                reprs.append((kind, desc))
        return reprs


class Viewpoint(Graphics3d):
    """
    This class represents a viewpoint, necessary for x3d.

    In the future, there could be multiple viewpoints, and they could have
    more properties. (Currently they only hold a position).
    """
    def __init__(self, *x):
        """
        EXAMPLES::

            sage: sage.plot.plot3d.base.Viewpoint(1, 2, 4).x3d_str()
            "<Viewpoint position='1 2 4'/>"
        """
        if isinstance(x[0], (tuple, list)):
            x = tuple(x[0])
        self.pos = x

    def x3d_str(self):
        """
        EXAMPLES::

            sage: sphere((0,0,0), 100).viewpoint().x3d_str()
            "<Viewpoint position='0 0 6'/>"
        """
        return "<Viewpoint position='%s %s %s'/>" % self.pos


cdef class PrimitiveObject(Graphics3d):
    """
    This is the base class for the non-container 3d objects.
    """
    def __init__(self, **kwds):
        if 'texture' in kwds:
            self.texture = kwds['texture']
            if not isinstance(self.texture, Texture):
                self.texture = Texture(self.texture)
        else:
            self.texture = Texture(kwds)

    def set_texture(self, texture=None, **kwds):
        """
        EXAMPLES::

            sage: G = dodecahedron(color='red'); G
            Graphics3d Object
            sage: G.set_texture(color='yellow'); G
            Graphics3d Object
        """
        if not isinstance(texture, Texture):
            texture = Texture(texture, **kwds)
        self.texture = texture

    def get_texture(self):
        """
        EXAMPLES::

            sage: G = dodecahedron(color='red')
            sage: G.get_texture()
            Texture(texture..., red, ff0000)
        """
        return self.texture

    def texture_set(self):
        """
        EXAMPLES::

            sage: G = dodecahedron(color='red')
            sage: G.texture_set()
            {Texture(texture..., red, ff0000)}
        """
        return set([self.texture])

    def x3d_str(self):
        r"""
        EXAMPLES::

            sage: sphere().flatten().x3d_str()
            "<Transform>\n<Shape><Sphere radius='1.0'/><Appearance><Material diffuseColor='0.4 0.4 1.0' shininess='1.0' specularColor='0.0 0.0 0.0'/></Appearance></Shape>\n\n</Transform>"
        """
        return "<Shape>" + self.x3d_geometry() + self.texture.x3d_str() + "</Shape>\n"

    def tachyon_repr(self, render_params):
        """
        Default behavior is to render the triangulation.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Torus
            sage: G = Torus(1, .5)
            sage: G.tachyon_repr(G.default_render_params())
            ['TRI V0 0 1 0.5
            ...
            'texture...']
        """
        return self.triangulation().tachyon_repr(render_params)

    def obj_repr(self, render_params):
        """
        Default behavior is to render the triangulation.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Torus
            sage: G = Torus(1, .5)
            sage: G.obj_repr(G.default_render_params())
            ['g obj_1',
             'usemtl ...',
             ['v 0 1 0.5',
             ...
              'f ...'],
             []]
        """
        return self.triangulation().obj_repr(render_params)

    def jmol_repr(self, render_params):
        r"""
        Default behavior is to render the triangulation. The actual polygon
        data is stored in a separate file.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Torus
            sage: G = Torus(1, .5)
            sage: G.jmol_repr(G.testing_render_params())
            ['pmesh obj_1 "obj_1.pmesh"\ncolor pmesh  [102,102,255]']
        """
        return self.triangulation().jmol_repr(render_params)

    def threejs_repr(self, render_params):
        r"""
        Default behavior is to render the triangulation.

        EXAMPLES::

            sage: from sage.plot.plot3d.base import PrimitiveObject
            sage: class SimpleTriangle(PrimitiveObject):
            ....:     def triangulation(self):
            ....:         return polygon3d([(0,0,0), (1,0,0), (0,1,0)])
            sage: G = SimpleTriangle()
            sage: G.threejs_repr(G.default_render_params())
            [('surface',
              {'color': '#0000ff',
               'faces': [[0, 1, 2]],
               'opacity': 1.0,
               'vertices': [{'x': 0.0, 'y': 0.0, 'z': 0.0},
                {'x': 1.0, 'y': 0.0, 'z': 0.0},
                {'x': 0.0, 'y': 1.0, 'z': 0.0}]})]
        """
        return self.triangulation().threejs_repr(render_params)


class BoundingSphere(SageObject):
    """
    A bounding sphere is like a bounding box, but is simpler to deal with and
    behaves better under rotations.
    """
    def __init__(self, cen, r):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.base import BoundingSphere
            sage: BoundingSphere((0,0,0), 1)
            Center (0.0, 0.0, 0.0) radius 1
            sage: BoundingSphere((0,-1,5), 2)
            Center (0.0, -1.0, 5.0) radius 2
            """
        self.cen = vector(RDF, cen)
        self.r = r

    def __repr__(self):
        """
        TESTS::

            sage: from sage.plot.plot3d.base import BoundingSphere
            sage: BoundingSphere((0,-1,10), 2)
            Center (0.0, -1.0, 10.0) radius 2
        """
        return "Center %s radius %s" % (self.cen, self.r)

    def __add__(self, other):
        """
        Return the bounding sphere containing both terms.

        EXAMPLES::

            sage: from sage.plot.plot3d.base import BoundingSphere
            sage: BoundingSphere((0,0,0), 1) + BoundingSphere((0,0,0), 2)
            Center (0.0, 0.0, 0.0) radius 2
            sage: BoundingSphere((0,0,0), 1) + BoundingSphere((0,0,100), 1)
            Center (0.0, 0.0, 50.0) radius 51.0
            sage: BoundingSphere((0,0,0), 1) + BoundingSphere((1,1,1), 2)
            Center (0.7886751345948128, 0.7886751345948128, 0.7886751345948128) radius 2.3660254037844384

        Treat None and 0 as the identity::

            sage: BoundingSphere((1,2,3), 10) + None + 0
            Center (1.0, 2.0, 3.0) radius 10
        """
        if other == 0 or other is None:
            return self
        elif self == 0 or self is None:
            return other
        if self.cen == other.cen:
            return self if self.r > other.r else other
        diff = other.cen - self.cen
        dist = (diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]).sqrt()
        diam = dist + self.r + other.r
        off  = diam/2 - self.r
        return BoundingSphere(self.cen + (off/dist)*diff, diam/2)

    def transform(self, T):
        """
        Return the bounding sphere of this sphere acted on by T. This always
        returns a new sphere, even if the resulting object is an ellipsoid.

        EXAMPLES::

            sage: from sage.plot.plot3d.transform import Transformation
            sage: from sage.plot.plot3d.base import BoundingSphere
            sage: BoundingSphere((0,0,0), 10).transform(Transformation(trans=(1,2,3)))
            Center (1.0, 2.0, 3.0) radius 10.0
            sage: BoundingSphere((0,0,0), 10).transform(Transformation(scale=(1/2, 1, 2)))
            Center (0.0, 0.0, 0.0) radius 20.0
            sage: BoundingSphere((0,0,3), 10).transform(Transformation(scale=(2, 2, 2)))
            Center (0.0, 0.0, 6.0) radius 20.0
        """
        return BoundingSphere(T.transform_point(self.cen), self.r * T.max_scale())


class RenderParams(SageObject):
    """
    This class is a container for all parameters that may be needed to
    render triangulate/render an object to a certain format. It can
    contain both cumulative and global parameters.

    Of particular note is the transformation object, which holds the
    cumulative transformation from the root of the scene graph to this
    node in the tree.
    """

    _uniq_counter = 0
    randomize_counter = 0
    force_reload = False
    mesh = False
    dots = False
    antialiasing = 8

    def __init__(self, **kwds):
        """
        EXAMPLES::

            sage: params = sage.plot.plot3d.base.RenderParams(foo='x')
            sage: params.transform_list
            []
            sage: params.foo
            'x'
        """
        self.output_file = tmp_filename()
        self.obj_vertex_offset = 1
        self.transform_list = []
        self.transform = None
        self.ds = 1
        self.crease_threshold = .8
        self.__dict__.update(kwds)
        # for jmol, some things (such as labels) must be attached to atoms
        self.atom_list = []

    def push_transform(self, T):
        """
        Push a transformation onto the stack, updating self.transform.

        EXAMPLES::

            sage: from sage.plot.plot3d.transform import Transformation
            sage: params = sage.plot.plot3d.base.RenderParams()
            sage: params.transform is None
            True
            sage: T = Transformation(scale=(10,20,30))
            sage: params.push_transform(T)
            sage: params.transform.get_matrix()
            [10.0  0.0  0.0  0.0]
            [ 0.0 20.0  0.0  0.0]
            [ 0.0  0.0 30.0  0.0]
            [ 0.0  0.0  0.0  1.0]
            sage: params.push_transform(T)  # scale again
            sage: params.transform.get_matrix()
            [100.0   0.0   0.0   0.0]
            [  0.0 400.0   0.0   0.0]
            [  0.0   0.0 900.0   0.0]
            [  0.0   0.0   0.0   1.0]
        """
        self.transform_list.append(self.transform)
        if self.transform is None:
            self.transform = T
        else:
            self.transform = self.transform * T

    def pop_transform(self):
        """
        Remove the last transformation off the stack, resetting self.transform
        to the previous value.

        EXAMPLES::

            sage: from sage.plot.plot3d.transform import Transformation
            sage: params = sage.plot.plot3d.base.RenderParams()
            sage: T = Transformation(trans=(100, 500, 0))
            sage: params.push_transform(T)
            sage: params.transform.get_matrix()
            [  1.0   0.0   0.0 100.0]
            [  0.0   1.0   0.0 500.0]
            [  0.0   0.0   1.0   0.0]
            [  0.0   0.0   0.0   1.0]
            sage: params.push_transform(Transformation(trans=(-100, 500, 200)))
            sage: params.transform.get_matrix()
            [   1.0    0.0    0.0    0.0]
            [   0.0    1.0    0.0 1000.0]
            [   0.0    0.0    1.0  200.0]
            [   0.0    0.0    0.0    1.0]
            sage: params.pop_transform()
            sage: params.transform.get_matrix()
            [  1.0   0.0   0.0 100.0]
            [  0.0   1.0   0.0 500.0]
            [  0.0   0.0   1.0   0.0]
            [  0.0   0.0   0.0   1.0]
        """
        self.transform = self.transform_list.pop()

    def unique_name(self, desc='name'):
        """
        Return a unique identifier starting with ``desc``.

        INPUT:

        - ``desc`` -- string (default: ``'name'``); the prefix of the names

        EXAMPLES::

            sage: params = sage.plot.plot3d.base.RenderParams()
            sage: params.unique_name()
            'name_1'
            sage: params.unique_name()
            'name_2'
            sage: params.unique_name('texture')
            'texture_3'
        """
        if self.randomize_counter:
            self._uniq_counter = randint(1,1000000)
        else:
            self._uniq_counter += 1
        return "%s_%s" % (desc, self._uniq_counter)


def flatten_list(L):
    """
    This is an optimized routine to turn a list of lists (of lists ...)
    into a single list. We generate data in a non-flat format to avoid
    multiple data copying, and then concatenate it all at the end.

    This is NOT recursive, otherwise there would be a lot of redundant
    copying (which we are trying to avoid in the first place, though at
    least it would be just the pointers).

    EXAMPLES::

        sage: from sage.plot.plot3d.base import flatten_list
        sage: flatten_list([])
        []
        sage: flatten_list([[[[]]]])
        []
        sage: flatten_list([['a', 'b'], 'c'])
        ['a', 'b', 'c']
        sage: flatten_list([['a'], [[['b'], 'c'], ['d'], [[['e', 'f', 'g']]]]])
        ['a', 'b', 'c', 'd', 'e', 'f', 'g']
    """
    if type(L) is not list:
        return [L]
    flat = []
    L_stack = []
    L_pop = L_stack.pop
    i_stack = []
    i_pop = i_stack.pop
    cdef Py_ssize_t i = 0
    while i < PyList_GET_SIZE(L) or PyList_GET_SIZE(L_stack) > 0:
        while i < PyList_GET_SIZE(L):
            tmp = <object>PyList_GET_ITEM(L, i)
            if type(tmp) is list:
                PyList_Append(L_stack, L)
                L = tmp
                PyList_Append(i_stack, i)
                i = 0
            else:
                PyList_Append(flat, tmp)
                i += 1
        if PyList_GET_SIZE(L_stack) > 0:
            L = L_pop()
            i = i_pop()
            i += 1
    return flat


def min3(v):
    """
    Return the componentwise minimum of a list of 3-tuples.

    EXAMPLES::

        sage: from sage.plot.plot3d.base import min3, max3
        sage: min3([(-1,2,5), (-3, 4, 2)])
        (-3, 2, 2)
    """
    return tuple([min([a[i] for a in v]) for i in range(3)])


def max3(v):
    """
    Return the componentwise maximum of a list of 3-tuples.

    EXAMPLES::

        sage: from sage.plot.plot3d.base import min3, max3
        sage: max3([(-1,2,5), (-3, 4, 2)])
        (-1, 4, 5)
    """
    return tuple([max([a[i] for a in v]) for i in range(3)])


def point_list_bounding_box(v):
    """
    Return the bounding box of a list of points.

    EXAMPLES::

        sage: from sage.plot.plot3d.base import point_list_bounding_box
        sage: point_list_bounding_box([(1,2,3),(4,5,6),(-10,0,10)])
        ((-10.0, 0.0, 3.0), (4.0, 5.0, 10.0))
        sage: point_list_bounding_box([(float('nan'), float('inf'), float('-inf')), (10,0,10)])
        ((10.0, 0.0, 10.0), (10.0, 0.0, 10.0))
    """
    cdef point_c low, high, cur
    low.x, low.y, low.z = INFINITY, INFINITY, INFINITY
    high.x, high.y, high.z = -INFINITY, -INFINITY, -INFINITY

    for P in v:
        cur.x, cur.y, cur.z = P
        point_c_update_finite_lower_bound(&low, cur)
        point_c_update_finite_upper_bound(&high, cur)
    return ((low.x, low.y, low.z), (high.x, high.y, high.z))


def optimal_aspect_ratios(ratios):
    """
    Average the aspect ratios.
    compute the elementwise maximum of triples.

    TESTS::

        sage: from sage.plot.plot3d.base import optimal_aspect_ratios
        sage: optimal_aspect_ratios([(2,4,6), (5,4,4), (1,2,7)])
        [5, 4, 7]
    """
    n = len(ratios)
    if n > 0:
        return [max([z[i] for z in ratios]) for i in range(3)]
    return [1.0, 1.0, 1.0]


def optimal_extra_kwds(v):
    """
    Merge a list v of dictionaries such that later
    dictionaries have precedence.

    TESTS::

        sage: from sage.plot.plot3d.base import optimal_extra_kwds
        sage: optimal_extra_kwds([{1:2, 2:3}, {2:4, 3:5}])
        {1: 2, 2: 4, 3: 5}
    """
    a = {}
    for b in v:
        a.update(b)
    return a


def _flip_orientation(v):
    """
    Switch from LH to RH coords to be consistent with Java rendition.

    TESTS::

        sage: from sage.plot.plot3d.base import _flip_orientation
        sage: _flip_orientation((1, 2, 3))
        (1, -2, 3)
    """
    return (v[0], -v[1], v[2])
