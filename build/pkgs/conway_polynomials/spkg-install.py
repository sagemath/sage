import os
import pickle

SAGE_SHARE = os.getenv('SAGE_SHARE')
install_root = os.path.join(SAGE_SHARE, 'conway_polynomials')

def create_db():
    db = {}
    from src import conway_polynomials
    for p, n, v in conway_polynomials:
        if not p in db:
            db[p] = {}
        db[p][n] = v
    if not os.path.exists(install_root):
        os.makedirs(install_root)
    with open(os.path.join(install_root, 'conway_polynomials.p'), 'wb') as f:
        pickle.dump(db, f)

if __name__ == '__main__':
    create_db()
