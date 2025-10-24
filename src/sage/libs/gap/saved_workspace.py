"""
LibGAP Workspace Support

The single purpose of this module is to provide the location of the
libgap saved workspace and a time stamp to invalidate saved
workspaces.
"""

import os
import glob
import subprocess
from sage.interfaces.gap_workspace import gap_workspace_file


def timestamp():
    """
    Return a time stamp for (lib)gap.

    OUTPUT:

    Float. Unix timestamp of the most recently changed GAP/LibGAP file(s). In particular, the
    timestamp increases whenever a gap package is added.

    EXAMPLES::

        sage: from sage.libs.gap.saved_workspace import timestamp
        sage: timestamp()   # random output
        1406642467.25684
        sage: type(timestamp())
        <... 'float'>
    """
    libgap_dir = os.path.dirname(__file__)
    libgap_files = glob.glob(os.path.join(libgap_dir, '*'))
    gap_packages = []
    
    # Try to get GAP root paths dynamically
    try:
        import shutil
        gap_exe = shutil.which('gap')
        if gap_exe:
            result = subprocess.run(
                [gap_exe, '-r', '-q', '--bare', '--nointeract', '-c',
                 'Display(JoinStringsWithSeparator(GAPInfo.RootPaths,";"));'],
                capture_output=True, text=True, timeout=5
            )
            if result.returncode == 0:
                gap_root_paths = result.stdout.strip()
                for d in gap_root_paths.split(";"):
                    if d:
                        gap_packages += glob.glob(os.path.join(d, 'pkg', '*'))
    except Exception:
        # If we can't get GAP paths, just use libgap files
        pass

    files = libgap_files + gap_packages
    if len(files) == 0:
        print('Unable to find LibGAP files.')
        return float('inf')
    return max(map(os.path.getmtime, files))


def workspace(name='workspace'):
    """
    Return the filename of the gap workspace and whether it is up to date.

    INPUT:

    - ``name`` -- string; a name that will become part of the
      workspace filename

    OUTPUT:

    Pair consisting of a string and a boolean. The string is the
    filename of the saved libgap workspace (or that it should have if
    it doesn't exist). The boolean is whether the workspace is
    up-to-date. You may use the workspace file only if the boolean is
    ``True``.

    EXAMPLES::

        sage: from sage.libs.gap.saved_workspace import workspace
        sage: ws, up_to_date = workspace()
        sage: ws
        '/.../gap/libgap-workspace-...'
        sage: isinstance(up_to_date, bool)
        True
    """
    workspace = gap_workspace_file("libgap", name)
    try:
        workspace_mtime = os.path.getmtime(workspace)
    except OSError:
        # workspace does not exist
        return (workspace, False)
    return (workspace, workspace_mtime >= timestamp())
