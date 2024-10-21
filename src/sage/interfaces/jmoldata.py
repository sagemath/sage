# sage_setup: distribution = sagemath-plot
r"""
Interface for extracting data and generating images from Jmol readable files.

JmolData is a no GUI version of Jmol useful for extracting data from files Jmol
reads and for generating image files.

AUTHORS:

- Jonathan Gutow (2012-06-14): complete doctest coverage
- Jonathan Gutow (2012-03-21): initial version
"""

# ******************************************************************************
#       Copyright (C) 2012 Jonathan Gutow (gutow@uwosh.edu)
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ******************************************************************************

from sage.structure.sage_object import SageObject

from sage.features.jmol import JmolDataJar
from sage.misc.temporary_file import tmp_filename
from sage.cpython.string import bytes_to_str

import os
import re
import subprocess
import sys
from pathlib import Path


class JmolData(SageObject):
    r"""
    .. TODO::

        Create an animated image file (GIF) if spin is on and put data
        extracted from a file into a variable/string/structure to return
    """
    def __init__(self):
        """
        EXAMPLES:

        Create a JmolData object::

            sage: from sage.interfaces.jmoldata import JmolData
            sage: JData = JmolData()
        """
        pass

    def is_jvm_available(self):
        """
        Return ``True`` if the Java Virtual Machine is available and ``False`` if not.

        EXAMPLES:

        Check that it returns a boolean::

            sage: from sage.interfaces.jmoldata import JmolData
            sage: JData = JmolData()
            sage: type(JData.is_jvm_available())
            <... 'bool'>
        """
        try:
            version = bytes_to_str(subprocess.check_output(['java', '-version'], stderr=subprocess.STDOUT))
        except (subprocess.CalledProcessError, OSError):
            return False

        java_version_number = int(re.sub(r'.*version "(0\.|1\.)?(\d*)[\s\S]*', r'\2', version, flags=re.S))
        return java_version_number >= 7

    def jmolpath(self):
        """
        Return the path to the jar file.

        EXAMPLES::

            sage: from sage.interfaces.jmoldata import JmolData
            sage: JData = JmolData()
            sage: JData.jmolpath()  # needs jmol
            '.../JmolData.jar'
        """
        jmolpath = JmolDataJar().absolute_filename()

        return jmolpath

    def is_jmol_available(self):
        """
        Return ``True`` if jmol is available and ``False`` if not.

        EXAMPLES:

        Check that it returns a boolean::

            sage: from sage.interfaces.jmoldata import JmolData
            sage: JData = JmolData()
            sage: type(JData.is_jmol_available())
            <... 'bool'>
        """
        if not JmolDataJar().is_present():
            return False

        if not self.is_jvm_available():
            return False

        return True

    def export_image(self,
                     targetfile,
                     datafile,  # name (path) of data file Jmol can read or script file telling it what to read or load
                     datafile_cmd='script',  # "script" or "load"
                     image_type='PNG',  # PNG, JPG, GIF
                     figsize=5,
                     **kwds):
        r"""
        This executes JmolData.jar to make an image file.

        INPUT:

        - ``targetfile`` -- the full path to the file where the image
          should be written

        - ``datafile`` -- full path to the data file Jmol can read or
          text of a script telling Jmol what to read or load

        - ``datafile_cmd`` -- (default: ``'script'``)  ``'load'`` or ``'script'``
          should be ``'load'`` for a data file

        - ``image_type`` -- (default: ``"PNG"``) ``'PNG'`` ``'JPG'`` or ``'GIF'``

        - ``figsize`` -- number (default: 5) equal to (pixels/side)/100

        OUTPUT: image file, .png, .gif or .jpg (default: .png)

        .. NOTE::

            Examples will generate an error message if a functional Java Virtual Machine (JVM)
            is not installed on the machine the Sage instance is running on.

        .. warning::

            Programmers using this module should check that the JVM is
            available before making calls to avoid the user getting
            error messages.  Check for the JVM using the function
            :meth:`is_jvm_available`, which returns ``True`` if a JVM is available.

        EXAMPLES:

        Use Jmol to load a pdb file containing some DNA from a web data
        base and make an image of the DNA. If you execute this in the
        notebook, the image will appear in the output cell::

            sage: from sage.interfaces.jmoldata import JmolData
            sage: JData = JmolData()
            sage: script = "load =1lcd;display DNA;moveto 0.0 { -473 -713 -518 59.94} 100.0 0.0 0.0 {21.17 26.72 27.295} 27.544636 {0.0 0.0 0.0} -25.287832 64.8414 0.0;"
            sage: testfile = tmp_filename(ext="DNA.png")
            sage: JData.export_image(targetfile=testfile,datafile=script,image_type="PNG")  # optional -- java internet
            sage: print(os.path.exists(testfile)) # optional -- java internet
            True

        Use Jmol to save an image of a 3-D object created in Sage.
        This method is used internally by plot3d to generate static images.
        This example doesn't have correct scaling::

            sage: from sage.interfaces.jmoldata import JmolData
            sage: JData = JmolData()
            sage: D = dodecahedron()                                                    # needs sage.plot
            sage: from tempfile import NamedTemporaryFile
            sage: archive = NamedTemporaryFile(suffix='.zip')
            sage: D.export_jmol(archive.name)                                           # needs sage.plot
            sage: archive_native = archive.name
            sage: script  = f'set defaultdirectory "f{archive_native}"\n'
            sage: script += 'script SCRIPT\n'
            sage: with NamedTemporaryFile(suffix='.png') as testfile:   # optional - java, needs sage.plot
            ....:     JData.export_image(targetfile=testfile.name,
            ....:                        datafile=script,
            ....:                        image_type="PNG")
            ....:     os.path.exists(testfile.name)
            True
            sage: archive.close()
        """
        # Set up paths, file names and scripts
        jmolpath = self.jmolpath()
        target_native = targetfile

        launchscript = ""
        if (datafile_cmd != 'script'):
            launchscript = "load "
        launchscript = launchscript + datafile

        imagescript = 'write {} {!r}\n'.format(image_type, target_native)
        size_arg = "%sx%s" % (figsize * 100, figsize * 100)
        # Scratch file for Jmol errors
        scratchout = tmp_filename(ext='.txt')
        with open(scratchout, 'w') as jout:
            # Now call the java application and write the file.
            env = dict(os.environ)
            env['LC_ALL'] = 'C'
            env['LANG'] = 'C'
            subprocess.call(["java", "-Xmx512m", "-Djava.awt.headless=true",
                             "-jar", jmolpath, "-iox", "-g", size_arg,
                             "-J", launchscript, "-j", imagescript],
                            stdout=jout, stderr=jout, env=env)
        if not os.path.isfile(targetfile):
            raise RuntimeError(f"Jmol failed to create file {targetfile}: {Path(scratchout).read_text()}")
        os.unlink(scratchout)
