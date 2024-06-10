import os
os.chdir("src")


try:
    SAGE_SHARE = os.environ['SAGE_SHARE']
except KeyError:
    raise RuntimeError("SAGE_SHARE undefined, maybe run `sage -sh`?")


def safe_makedirs(dir):
    try:
        os.makedirs(dir)
    except OSError:
        if not os.path.isdir(dir):
            raise


def install_cremona():
    from sqlite3 import connect

    cremona_root = os.path.join(SAGE_SHARE, 'cremona')
    safe_makedirs(cremona_root)

    target = os.path.join(cremona_root, 'cremona_mini.db')

    print("Creating database {0}".format(target))
    if os.path.exists(target):
        os.remove(target)

    con = connect(target)

    con.execute('CREATE TABLE t_class(rank INTEGER, class TEXT PRIMARY KEY,'
            ' conductor INTEGER)')
    con.execute('CREATE TABLE t_curve(curve TEXT PRIMARY KEY, class TEXT, tors'
            ' INTEGER, eqn TEXT UNIQUE)')
    con.execute('CREATE INDEX i_t_class_conductor ON t_class(conductor)')
    con.execute('CREATE INDEX i_t_curve_class ON t_curve(class)')

    class_data = []
    curve_data = []

    for line in open(os.path.join('common', 'allcurves.00000-09999')):
        N, iso, num, eqn, r, tors = line.split()
        cls = N + iso
        cur = cls + num
        if num == "1":
            class_data.append((N, cls, r))
        curve_data.append((cur, cls, eqn, tors))

    con.executemany('INSERT INTO t_class(conductor,class,rank) VALUES'
            ' (?,?,?)', class_data)
    con.executemany('INSERT INTO t_curve(curve,class,eqn,tors) VALUES'
            ' (?,?,?,?)', curve_data)

    con.commit()


def install_ellcurves():
    ellcurves_root = os.path.join(SAGE_SHARE, 'ellcurves')

    # Remove previous installs (possibly with bad permissions, see
    # https://github.com/sagemath/sage/issues/21641)
    import shutil
    try:
        shutil.rmtree(ellcurves_root)
    except OSError:
        pass

    safe_makedirs(ellcurves_root)

    # Open all source files from the "ellcurves" directory containing
    # a database of curves with given rank. In the install directory,
    # we store all curves from the "allcurves" Cremona database by rank,
    # as well as the curves from "ellcurves".
    rank_src = {'rank0': None, 'rank1': None, 'rank2': None}  # Not in "ellcurves"
    for r in os.listdir('ellcurves'):
        if not r.startswith("."):
            rank_src[r] = open(os.path.join('ellcurves', r))

    # Create files for each rank of elliptic curves in the databases.
    rank_dst = {}
    for r in rank_src:
        filename = os.path.join(ellcurves_root, r)
        print("Creating text file {0}".format(filename))
        rank_dst[r] = open(filename, 'w')

    # Write the curves from the "allcurves" database to the rank file
    for line in open(os.path.join('common', 'allcurves.00000-09999')):
        r = "rank" + line.split()[4]
        rank_dst[r].write(line)

    # Now append the ranks from the "ellcurves" database
    for r, f in rank_src.items():
        if f is not None:
            rank_dst[r].write(f.read())


if __name__ == '__main__':
    install_cremona()
    install_ellcurves()
