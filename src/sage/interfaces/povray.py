# sage_setup: distribution = sagemath-plot
r"""
POV-Ray, The Persistence of Vision Ray Tracer
"""

from sage.misc.pager import pager
import os


class POVRay:
    """
    POV-Ray The Persistence of Vision Ray Tracer

    INPUT:

    - ``pov_file`` -- complete path to the .pov file you want to be rendered
    - ``outfile`` -- the filename you want to save your result to
    - ``**kwargs`` -- additionally keyword arguments you want to pass to POVRay

    OUTPUT: image is written to the file you specified in outfile

    EXAMPLES:

    AUTHOR:

    Sage interface written by Yi Qiang (yqiang _atNOSPAM_ gmail.com)

    POVRay: http://www.povray.org
    """
    def __repr__(self):
        return 'POV-Ray The Persistence of Vision Ray Tracer'

    def __call__(self, pov_file, outfile='sage.ppm', block=True, **kwargs):
        if not os.path.isfile(pov_file):
            return "%s not found" % (pov_file)

        outfile = os.path.abspath(os.path.expanduser(outfile))

        if not ('W' in kwargs and 'H' in kwargs):
            return "You must specify a width and height."

        cmd = "povray -D +FP +I%s +O%s " % (pov_file, outfile)
        for k, v in kwargs.items():
            cmd += "+%s%s " % (k, v)

        if not block:
            cmd += ' &'
        os.system(cmd)

    def usage(self):
        with os.popen('povray') as f:
            r = f.read()
        pager()(r)


povray = POVRay()
