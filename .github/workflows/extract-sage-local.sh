#!/bin/dash
# to be run from $SAGE_ROOT, with arguments sage-local-${{ env.PREVIOUS_STAGES }}.tar

if [ -z "$SAGE_LOCAL" ]; then
    SAGE_LOCAL=$(pwd)/local
fi

# Show all tar files
ls -l $*

for a in $*; do
    echo Extracting $a
    (cd / && tar xf -) < $a
    rm -f $a
done

# We set the installation records to the same mtime so that no rebuilds due to dependencies
# among these packages are triggered.
dummy="$SAGE_LOCAL"/var/lib/sage/installed/.dummy
if [ -f "$dummy" ]; then
    touch "$dummy"
    for tree in "$SAGE_LOCAL" "$SAGE_LOCAL"/var/lib/sage/venv*; do
        inst="$tree"/var/lib/sage/installed
        if [ -d "$inst" ]; then
            # -r is --reference; the macOS version of touch does not accept the long option.
            (cd "$inst" && touch -r "$dummy" .dummy *)
            # Show what has been built already.
            ls -l "$tree" "$inst"
        fi
    done
fi

# Show how we are doing on free space.
df -h
