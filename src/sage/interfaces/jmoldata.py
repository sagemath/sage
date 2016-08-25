r"""
Interface for extracting data and generating images from Jmol readable files.

JmolData is a no GUI version of Jmol useful for extracting data from files Jmol
reads and for generating image files.

AUTHORS:

- Jonathan Gutow (2012-06-14): complete doctest coverage
- Jonathan Gutow (2012-03-21): initial version
"""

#*******************************************************************************
#       Copyright (C) 2012 Jonathan Gutow (gutow@uwosh.edu)
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*******************************************************************************
from __future__ import print_function

from sage.structure.sage_object import SageObject

from sage.env import SAGE_LOCAL
from sage.misc.temporary_file import tmp_filename

import subprocess
import os

class JmolData(SageObject):
    r"""
    .. todo::

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
        Returns True if the Java Virtual Machine is available and False if not.

        EXAMPLES:

        Check that it returns a boolean::

            sage: from sage.interfaces.jmoldata import JmolData
            sage: JData = JmolData()
            sage: type(JData.is_jvm_available())
            <type 'bool'>
        """
        try:
            version = subprocess.check_output(['java', '-version'], stderr=subprocess.STDOUT)
        except (subprocess.CalledProcessError, OSError):
            return False

        import re
        java_version = re.search("version.*[1][.][78]", version)
        return java_version is not None

    def export_image(self,
        targetfile,
        datafile, #name (path) of data file Jmol can read or script file telling it what to read or load
        datafile_cmd='script', #"script" or "load"
        image_type ='PNG', #PNG, JPG, GIF
        figsize=5,
        **kwds):
        r"""
        This executes JmolData.jar to make an image file.

        INPUT:

        - targetfile -- the full path to the file where the image
          should be written.

        - datafile -- full path to the data file Jmol can read or
          text of a script telling Jmol what to read or load.
          If it is a script and the platform is cygwin, the filenames in
          the script should be in native windows format.

        - datafile_cmd -- (default ``'script'``)  ``'load'`` or ``'script'``
          should be ``"load"`` for a data file.

        - image_type -- (default ``"PNG"``) ``'PNG'`` ``'JPG'`` or ``'GIF'``

        - figsize -- number (default 5) equal to (pixels/side)/100

        OUTPUT:

        Image file, .png, .gif or .jpg (default .png)

        .. note::

            Examples will generate an error message if a functional Java Virtual Machine (JVM)
            is not installed on the machine the Sage instance is running on.

        .. warning::

            Programmers using this module should check that the JVM is
            available before making calls to avoid the user getting
            error messages.  Check for the JVM using the function
            :meth:`is_jvm_available`, which returns True if a JVM is available.

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
            sage: D=dodecahedron()
            sage: from sage.misc.misc import SAGE_TMP
            sage: archive_name=os.path.join(SAGE_TMP, "archive.jmol.zip")
            sage: D.export_jmol(archive_name)  #not scaled properly...need some more steps.
            sage: archive_native = archive_name
            sage: import sys
            sage: if sys.platform == 'cygwin':
            ....:     from subprocess import check_output, STDOUT
            ....:     archive_native = check_output(['cygpath', '-w', archive_native],
            ....:                                   stderr=STDOUT).rstrip()
            sage: script = 'set defaultdirectory "{0}"\n script SCRIPT\n'.format(archive_native)
            sage: testfile = os.path.join(SAGE_TMP, "testimage.png")
            sage: JData.export_image(targetfile=testfile, datafile=script, image_type="PNG") # optional -- java
            sage: print(os.path.exists(testfile)) # optional -- java
            True
        """
        # Set up paths, file names and scripts
        jmolpath = os.path.join(SAGE_LOCAL, "share", "jmol", "JmolData.jar")
        target_native = targetfile
        import sys
        if sys.platform == 'cygwin':
            jmolpath      = subprocess.check_output(['cygpath', '-w', jmolpath],
                                                    stderr=subprocess.STDOUT).rstrip()
            target_native = subprocess.check_output(['cygpath', '-w', target_native],
                                                    stderr=subprocess.STDOUT).rstrip()
            if (datafile_cmd != 'script'):
                datafile  = subprocess.check_output(['cygpath', '-w', datafile],
                                                    stderr=subprocess.STDOUT).rstrip()
        launchscript = ""
        if (datafile_cmd!='script'):
            launchscript = "load "
        launchscript = launchscript + datafile
        imagescript = "write " + image_type + " " + target_native + "\n"

        sizeStr = "%sx%s" %(figsize*100,figsize*100)
        # Scratch file for Jmol errors
        scratchout = tmp_filename(ext=".txt")
        with open(scratchout, 'w') as jout:
            # Now call the java application and write the file.
            env = dict(os.environ)
            env['LC_ALL'] = 'C'
            env['LANG'] = 'C'
            subprocess.call(["java", "-Xmx512m", "-Djava.awt.headless=true",
                "-jar", jmolpath, "-iox", "-g", sizeStr,
                "-J", launchscript, "-j", imagescript],
                stdout=jout, stderr=jout, env=env)
        if not os.path.isfile(targetfile):
            raise RuntimeError("Jmol failed to create file %s, see %s for details"%(repr(targetfile), repr(scratchout)))
        os.unlink(scratchout)
