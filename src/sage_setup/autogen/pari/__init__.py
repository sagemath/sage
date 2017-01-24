import glob
import os
from os.path import join, getmtime, exists

from sage_setup.autogen.pari.generator import PariFunctionGenerator
from sage_setup.autogen.pari.parser import pari_share, sage_src_pari


def rebuild(force=False):
    pari_module_path = sage_src_pari()
    src_files = [join(pari_share(), b'pari.desc')] + \
                 glob.glob(join('sage_setup', 'autogen', 'pari', '*.py'))
    gen_files = [join(pari_module_path, 'auto_gen.pxi')]

    if all(exists(f) for f in gen_files):
        src_mtime = max(getmtime(f) for f in src_files)
        gen_mtime = min(getmtime(f) for f in gen_files)

        if gen_mtime > src_mtime and not force:
            return

    # Clean up old auto-generated files
    old_pari_module_path = join('sage', 'libs', 'pari')
    for filename in glob.glob(join(old_pari_module_path, 'auto_*.pxi')):
        os.remove(filename)

    G = PariFunctionGenerator()
    G()
