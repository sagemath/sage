cd local/share/doc/sage/latex
echo "All missing character error messages in LaTeX log files:"
num_errors=$(grep -r "Missing character" --include "*.log" | tee /dev/stderr | wc -l)
expected_num_errors=1
echo "In total there are $num_errors missing character errors, expecting $expected_num_errors"
[ $num_errors = $expected_num_errors ]
