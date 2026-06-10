import pathlib
import re
import sys

check = "-c" in sys.argv

path = pathlib.Path("src/doc/en/reference/references/index.rst")
text = path.read_text()

blocks = re.split(r'(?=^\.\. \[)', text, flags=re.MULTILINE)
blocks = [b for b in blocks if b.strip()]

def key(block):
    m = re.match(r'\.\. \[([^\]]+)\]', block)
    tag = m.group(1) if m else ""
    return (tag.casefold(), tag)

sorted_blocks = sorted(blocks, key=key)

if check:
    if blocks == sorted_blocks:
        sys.exit(0)
    orig_keys = [key(b)[1] for b in blocks]
    sort_keys = [key(b)[1] for b in sorted_blocks]
    for o, s in zip(orig_keys, sort_keys):
        if o != s:
            print(f"sort-references: disorder: [{o}]")
            sys.exit(1)

# Format in-place
path.write_text("".join(sorted_blocks))
