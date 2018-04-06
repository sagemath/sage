import os, sys

os.system("echo 'import scipy; scipy.test()' | sage-python23 > test.out 2>test.err")

r = open('test.err').read()

if 'FAILED' in r:
    open('test.out').read()
    print r
    sys.exit(1)
else:
    os.unlink('test.err')
    os.unlink('test.out')
    sys.exit(0)
