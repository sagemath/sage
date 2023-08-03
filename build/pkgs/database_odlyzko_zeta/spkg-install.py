import os
from sage.all import save
from sage.env import SAGE_SHARE

install_root = os.path.join(SAGE_SHARE, 'odlyzko')
target = os.path.join(install_root, 'zeros.sobj')

if __name__ == '__main__':
    os.makedirs(install_root, exist_ok=True)
    print("Creating Odlyzko database.")
    F = [float(x) for x in open("src/zeros6").readlines()]
    save(F, target)
