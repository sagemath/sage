#!/usr/bin/env bash
# to be run from $SAGE_ROOT, with arguments sage-local-${{ env.PREVIOUS_STAGES }}.tar

# Show all tar files
ls -l $*

# We specifically use the cygwin tar so that symlinks are saved/restored correctly on Windows.
for a in $*; do
    echo Extracting $a
    tar xf $a
    rm -f $a
done

# Also get rid of the stages that were not extracted
rm -f sage-local-*.tar

# We set the installation records to the same mtime so that no rebuilds due to dependencies
# among these packages are triggered.
(cd local/var/lib/sage/installed/ && touch .dummy && touch --reference=.dummy *)

# Show what has been built already.
ls -l local local/var/lib/sage/installed/
df -h

# Rebase!
src/bin/sage-rebase.sh local
