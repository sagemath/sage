#!/bin/dash
# to be run from $SAGE_ROOT, with arguments sage-local-${{ env.PREVIOUS_STAGES }}.tar

if [ -z "$SAGE_LOCAL" ]; then
    SAGE_LOCAL=$(pwd)/local
fi

# Show all tar files
ls -l $*

# We specifically use the cygwin tar so that symlinks are saved/restored correctly on Windows.
for a in $*; do
    echo Extracting $a
    (cd / && tar xf -) < $a
    rm -f $a
done

# We set the installation records to the same mtime so that no rebuilds due to dependencies
# among these packages are triggered.
(cd "$SAGE_LOCAL"/var/lib/sage/installed/ && touch .dummy && touch --reference=.dummy *)

# Show what has been built already.
ls -l "$SAGE_LOCAL" "$SAGE_LOCAL"/var/lib/sage/installed/
df -h

# Rebase!
exec src/bin/sage-rebase.sh --all "$SAGE_LOCAL"
