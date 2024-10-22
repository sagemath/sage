import _sage_conf
from _sage_conf._conf import *
from _sage_conf.__main__ import _main


# Relocation.  SAGE_ROOT is typically configured to a directory in /var/tmp (sticky).
import os as _os
SAGE_ROOT_ABS = _os.path.join(os.path.dirname(_sage_conf.__file__), 'sage_root')

try:
    _os.symlink(SAGE_ROOT_ABS, SAGE_ROOT, target_is_directory=True)
except FileExistsError:
    pass

if _os.stat(SAGE_ROOT).st_uid != _os.geteuid():
    raise RuntimeError(f'Cannot create symlink {SAGE_ROOT} - it is owned by another user')
