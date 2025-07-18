# sage_setup: distribution = sagemath-environment
r"""
Feature for testing the presence of ``ffmpeg``
"""
# ****************************************************************************
#       Copyright (C) 2018-2022 Sebastien Labbe <slabqc@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from . import Executable, FeatureTestResult


class FFmpeg(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :ref:`ffmpeg <spkg_ffmpeg>`.

    EXAMPLES::

        sage: from sage.features.ffmpeg import FFmpeg
        sage: FFmpeg().is_present()  # optional - ffmpeg
        FeatureTestResult('ffmpeg', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.ffmpeg import FFmpeg
            sage: isinstance(FFmpeg(), FFmpeg)
            True
        """
        Executable.__init__(self, 'ffmpeg', executable='ffmpeg',
                            spkg='ffmpeg',
                            url='https://www.ffmpeg.org/')

    def is_functional(self):
        r"""
        Return whether command ``ffmpeg`` in the path is functional.

        EXAMPLES::

            sage: from sage.features.ffmpeg import FFmpeg
            sage: FFmpeg().is_functional()   # optional - ffmpeg
            FeatureTestResult('ffmpeg', True)
        """
        # Create the content of 1-pixel png file
        content = b'\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01\x08\x00\x00\x00\x00:~\x9bU\x00\x00\x00\nIDATx\x9cc`\x00\x00\x00\x02\x00\x01H\xaf\xa4q\x00\x00\x00\x00IEND\xaeB`\x82'

        # NOTE:
        #
        # This is how the above content of a 1 pixel PNG was created::
        #
        #    sage: import numpy as np
        #    sage: from PIL import Image
        #    sage: image = Image.fromarray(np.array([[100]], dtype=np.uint8))
        #    sage: image.save('file.png')
        #    sage: with open('file.png', 'rb') as f:
        #    ....:     content = f.read()

        # create a png file with the content
        from sage.misc.temporary_file import tmp_filename
        base_filename_png = tmp_filename(ext='.png')
        with open(base_filename_png, 'wb') as f:
            f.write(content)

        # Set up filenames
        import os
        base, filename_png = os.path.split(base_filename_png)
        filename, _png = os.path.splitext(filename_png)

        # Setting a list of commands (taken from sage/plot/animate.py)
        # The `-nostdin` is needed to avoid the command to hang, see
        # https://stackoverflow.com/questions/16523746/ffmpeg-hangs-when-run-in-background
        commands = []
        for ext in ['.avi', '.flv', '.gif', '.mkv', '.mov',
                    '.mp4', '.ogg', '.ogv', '.webm', '.wmv']:

            cmd = ['ffmpeg', '-nostdin', '-y', '-f', 'image2', '-r', '5',
                   '-i', filename_png, '-pix_fmt', 'rgb24', '-loop', '0',
                   filename + ext]
            commands.append(cmd)

        for ext in ['.avi', '.flv', '.gif', '.mkv', '.mov', '.mpg',
                    '.mp4', '.ogg', '.ogv', '.webm', '.wmv']:

            cmd = ['ffmpeg', '-nostdin', '-y', '-f', 'image2', '-i',
                   filename_png, filename + ext]
            commands.append(cmd)

        # Running the commands and reporting any issue encountered
        from subprocess import run
        for cmd in commands:
            try:
                result = run(cmd, cwd=base, capture_output=True, text=True)
            except OSError as e:
                return FeatureTestResult(self, False, reason='Running command "{}" '
                            'raised an OSError "{}" '.format(' '.join(cmd), e))

            # If an error occurred, return False
            if result.returncode:
                return FeatureTestResult(self, False, reason='Running command "{}" '
                            'returned nonzero exit status "{}" with stderr '
                            '"{}" and stdout "{}".'.format(result.args,
                                                            result.returncode,
                                                            result.stderr.strip(),
                                                            result.stdout.strip()))

        # If necessary, run more tests here
        # ...

        # The command seems functional
        return FeatureTestResult(self, True)


def all_features():
    return [FFmpeg()]
