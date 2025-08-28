#!/bin/sh
expected_num_errors=1

cd builddir/src/doc/latex
echo "All missing character error messages in LaTeX log files:"
num_errors=$(grep -r "Missing character" --include "*.log" | tee /dev/stderr | wc -l)
echo "In total there are $num_errors missing character errors, expecting $expected_num_errors"
[ $num_errors = $expected_num_errors ]
