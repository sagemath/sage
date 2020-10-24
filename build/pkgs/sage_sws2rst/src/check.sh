#! /bin/sh -x
set -e
sage-sws2rst -h > /dev/null
for a in sage_sws2rst/comments2rst.py sage_sws2rst/results2rst.py; do python3 $a; done
(cd test && for a in *.sws; do sage-sws2rst "$a"; done)
