# "sage -python -m sage_docbuild.vendor" updates src/doc/common/_vendor/*.inv
#
import os
import sys

import requests

from .conf import _vendored_inventories_dir, _intersphinx_targets


if __name__ == '__main__':
    if not _vendored_inventories_dir:
        print('Error: sage_docbuild.vendor needs to be able to write to SAGE_SRC', file=sys.stderr)
        sys.exit(1)
    errors = 0
    for key, targets in _intersphinx_targets.items():
        for target in targets:
            if target.startswith('http'):
                inv_url = target + 'objects.inv'
                fname = os.path.join(_vendored_inventories_dir, key + '.inv')
                print(f'Requesting {inv_url}', flush=True)
                try:
                    r = requests.get(inv_url)
                    with open(fname, 'wb') as fd:
                        fd.write(r.content)
                except Exception as e:
                    print(f'Error: {e}')
                    errors += 1
                else:
                    print(f'Updated {fname}')
                    break
    sys.exit(errors)
