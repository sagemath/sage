#!/bin/bash

# This script traverses files with .py and .pyx suffixes in src/<directory>
# tree of the Sage library and store the module names sorted by their path
# and titles.
#
# EXAMPLE:
#
# $ src/doc/en/reference/make_module_list.sh sage/combinat

if [ $# -ne 1 ]; then
  echo "Usage: $0 <directory>"
  echo "Example: $0 sage/combinat"
  exit 1
fi

module_dir="$1"
input_dir="$SAGE_ROOT/src/$1"

if [ ! -d "$input_dir" ]; then
  echo "Error: Directory '$input_dir' not found."
  exit 1
fi

tmpfile=$(mktemp)

cd "$input_dir"

# Find all .py and .pyx files inside the given directory, excluding __init__.py
find . -type f \( -name "*.py" -o -name "*.pyx" \) ! -name "__init__.py" | while read -r file; do
  # Extract first line after the first triple quotes
  title=$(awk 'BEGIN {found=0}
               /^[rR]?"""/ {if (found == 0) {found=1; next}}
               found && NF {print; exit}' "$file")

  # Format module name: Remove .py or .pyx, replace directory prefix
  mod_name=$(echo "$file" | sed -E "s/\.(pyx|py)$//; s|^\./|    $module_dir/|")
  path_name=$(dirname "$file")

  # Store sorting path/title and module name in temporary file
  printf "%s\t%s\n" "$path_name,$title" "$mod_name" >> "$tmpfile"
done

# Sort and output module names to /tmp/module_list.rst
LC_ALL=C sort "$tmpfile" | cut -f2 > /tmp/module_list

rm -f "$tmpfile"

echo "Sorted module list saved to /tmp/module_list"
