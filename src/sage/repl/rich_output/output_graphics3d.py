# sage_setup: distribution = sagemath-repl
r"""
Three-Dimensional Graphics Output Types

This module defines the rich output types for 3-d scenes.
"""
# ****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import os
import importlib.resources

from sage.cpython.string import bytes_to_str, FS_ENCODING
from sage.repl.rich_output.output_basic import OutputBase
from sage.repl.rich_output.buffer import OutputBuffer


class OutputSceneJmol(OutputBase):

    def __init__(self, scene_zip, preview_png):
        """
        JMol Scene.

        By our (Sage) convention, the actual scene is called ``SCENE``
        inside the zip archive.

        INPUT:

        - ``scene_zip`` -- string/bytes; the jmol scene (a zip archive)

        - ``preview_png`` -- string/bytes; preview as png file

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputSceneJmol
            sage: OutputSceneJmol.example()
            OutputSceneJmol container
        """
        self.scene_zip = OutputBuffer(scene_zip)
        self.preview_png = OutputBuffer(preview_png)

    def launch_script_filename(self):
        """
        Return a launch script suitable to display the scene.

        This method saves the scene to disk and creates a launch
        script. The latter contains an absolute path to the scene
        file. The launch script is often necessary to make jmol
        render the 3d scene.

        OUTPUT: string; the file name of a suitable launch script

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputSceneJmol
            sage: rich_output = OutputSceneJmol.example();  rich_output
            OutputSceneJmol container
            sage: filename = rich_output.launch_script_filename();  filename
            '/.../scene.spt'
            sage: with open(filename) as fobj:
            ....:     print(fobj.read())
            set defaultdirectory "/.../scene.spt.zip"
            script SCRIPT
        """
        from sage.misc.temporary_file import tmp_dir
        basedir = tmp_dir()
        scene_filename = os.path.join(basedir, 'scene.spt.zip')
        script_filename = os.path.join(basedir, 'scene.spt')
        self.scene_zip.save_as(scene_filename)
        with open(script_filename, 'w') as f:
            f.write('set defaultdirectory "{0}"\n'.format(scene_filename))
            f.write('script SCRIPT\n')
        return script_filename

    @classmethod
    def example(cls):
        r"""
        Construct a sample Jmol output container.

        This static method is meant for doctests, so they can easily
        construct an example.

        OUTPUT: an instance of :class:`OutputSceneJmol`

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputSceneJmol
            sage: rich_output = OutputSceneJmol.example();  rich_output
            OutputSceneJmol container

            sage: rich_output.scene_zip
            buffer containing 654 bytes
            sage: rich_output.scene_zip.get().startswith(b'PK')
            True

            sage: rich_output.preview_png
            buffer containing 608 bytes
            sage: rich_output.preview_png.get().startswith(b'\x89PNG')
            True
        """
        example_png = importlib.resources.read_binary(__package__, 'example.png')
        scene_zip = importlib.resources.read_binary(__package__, 'example_jmol.spt.zip')
        return cls(scene_zip, example_png)


class OutputSceneCanvas3d(OutputBase):

    def __init__(self, canvas3d):
        """
        Canvas3d Scene.

        INPUT:

        - ``canvas3d`` -- string/bytes; the canvas3d data

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputSceneCanvas3d
            sage: OutputSceneCanvas3d.example()
            OutputSceneCanvas3d container
        """
        self.canvas3d = OutputBuffer(canvas3d)

    @classmethod
    def example(cls):
        r"""
        Construct a sample Canvas3D output container.

        This static method is meant for doctests, so they can easily
        construct an example.

        OUTPUT: an instance of :class:`OutputSceneCanvas3d`

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputSceneCanvas3d
            sage: rich_output = OutputSceneCanvas3d.example();  rich_output
            OutputSceneCanvas3d container

            sage: rich_output.canvas3d
            buffer containing 829 bytes
            sage: rich_output.canvas3d.get_str()
            '[{"vertices":[{"x":1,"y":1,"z":1},...{"x":1,"y":-1,"z":-1}],"faces":[[0,1,2,3]],"color":"008000"}]'
        """
        with importlib.resources.path(__package__, 'example.canvas3d') as filename:
            return cls(OutputBuffer.from_file(filename))


class OutputSceneThreejs(OutputBase):

    def __init__(self, html):
        """
        Three.js Scene.

        INPUT:

        - ``html`` -- string/bytes; the Three.js HTML data

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputSceneThreejs
            sage: OutputSceneThreejs('<html></html>')
            OutputSceneThreejs container
        """
        self.html = OutputBuffer(html)


class OutputSceneWavefront(OutputBase):

    def __init__(self, obj, mtl):
        """
        Wavefront `*.obj` Scene.

        The Wavefront format consists of two files, an ``.obj`` file
        defining the geometry data (mesh points, normal vectors, ...)
        together with a ``.mtl`` file defining texture data.

        INPUT:

        - ``obj`` -- bytes; the Wavefront obj file format describing
          the mesh shape

        - ``mtl`` -- bytes; the Wavefront mtl file format describing
          textures

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputSceneWavefront
            sage: OutputSceneWavefront.example()
            OutputSceneWavefront container
        """
        self.obj = OutputBuffer(obj)
        self.mtl = OutputBuffer(mtl)
        self._check_no_directory(self.mtllib())

    def _check_no_directory(self, filename):
        """
        Verify that ``filename`` does not contain a path.

        We disallow anything but plain filenames since it is a
        potential security issue to point :meth:`mtllib` at random
        paths.

        INPUT:

        - ``filename`` -- string; a filename

        OUTPUT:

        This method returns nothing. A :exc:`ValueError` is raised if
        ``filename`` is not just a plain filename but contains a
        directory (relative or absolute).

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputSceneWavefront
            sage: rich_output = OutputSceneWavefront.example()
            sage: rich_output._check_no_directory('scene.mtl')
            sage: rich_output._check_no_directory('/scene.mtl')
            Traceback (most recent call last):
            ...
            ValueError: must be pure filename, got directory component: /scene.mtl
            sage: rich_output._check_no_directory('relative/scene.mtl')
            Traceback (most recent call last):
            ...
            ValueError: must be pure filename, got directory component: relative/scene.mtl
            sage: rich_output._check_no_directory('/absolute/scene.mtl')
            Traceback (most recent call last):
            ...
            ValueError: must be pure filename, got directory component: /absolute/scene.mtl
        """
        if os.path.split(filename)[0]:
            raise ValueError('must be pure filename, got directory component: {0}'
                             .format(filename))

    def mtllib(self):
        """
        Return the ``mtllib`` filename.

        The ``mtllib`` line in the Wavefront file format (``*.obj``)
        is the name of the separate texture file.

        OUTPUT:

        String. The filename under which ``mtl`` is supposed to be
        saved.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputSceneWavefront
            sage: rich_output = OutputSceneWavefront.example()
            sage: rich_output.mtllib()
            'scene.mtl'
        """
        marker = b'mtllib '
        for line in self.obj.get().splitlines():
            if line.startswith(marker):
                return bytes_to_str(line[len(marker):], FS_ENCODING,
                                    'surrogateescape')
        return 'scene.mtl'

    def obj_filename(self):
        """
        Return the file name of the ``.obj`` file.

        This method saves the object and texture to separate files in
        a temporary directory and returns the object file name. This
        is often used to launch a 3d viewer.

        OUTPUT:

        String. The file name (absolute path) of the saved obj file.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputSceneWavefront
            sage: rich_output = OutputSceneWavefront.example();  rich_output
            OutputSceneWavefront container
            sage: obj = rich_output.obj_filename();  obj
            '/.../scene.obj'
            sage: with open(obj) as fobj:
            ....:     print(fobj.read())
            mtllib scene.mtl
            g obj_1
            ...
            f 3 2 6 8

            sage: path = os.path.dirname(obj)
            sage: mtl = os.path.join(path, 'scene.mtl');  mtl
            '/.../scene.mtl'
            sage: os.path.exists(mtl)
            True
            sage: os.path.dirname(obj) == os.path.dirname(mtl)
            True
            sage: with open(mtl) as fobj:
            ....:     print(fobj.read())
            newmtl texture177
            Ka 0.2 0.2 0.5
            ...
            d 1
        """
        from sage.misc.temporary_file import tmp_dir
        basedir = tmp_dir()
        obj_filename = os.path.join(basedir, 'scene.obj')
        mtl_filename = os.path.join(basedir, self.mtllib())
        self.obj.save_as(obj_filename)
        self.mtl.save_as(mtl_filename)
        return os.path.abspath(obj_filename)

    @classmethod
    def example(cls):
        r"""
        Construct a sample Canvas3D output container.

        This static method is meant for doctests, so they can easily
        construct an example.

        OUTPUT: an instance of :class:`OutputSceneCanvas3d`

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputSceneWavefront
            sage: rich_output = OutputSceneWavefront.example();  rich_output
            OutputSceneWavefront container

            sage: rich_output.obj
            buffer containing 227 bytes
            sage: rich_output.obj.get_str()
            'mtllib scene.mtl\ng obj_1\n...\nf 1 5 6 2\nf 1 4 7 5\nf 6 5 7 8\nf 7 4 3 8\nf 3 2 6 8\n'

            sage: rich_output.mtl
            buffer containing 80 bytes
            sage: rich_output.mtl.get_str()
            'newmtl texture177\nKa 0.2 0.2 0.5\nKd 0.4 0.4 1.0\nKs 0.0 0.0 0.0\nillum 1\nNs 1\nd 1\n'
        """
        with importlib.resources.path(__package__, 'example_wavefront_scene.obj') as filename:
            scene_obj = OutputBuffer.from_file(filename)
        with importlib.resources.path(__package__, 'example_wavefront_scene.mtl') as filename:
            scene_mtl = OutputBuffer.from_file(filename)
        return cls(scene_obj, scene_mtl)
