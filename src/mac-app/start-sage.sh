#!/bin/bash

# start-sage.sh
#
# Created by Ivan Andrus on 16/1/10.
# Copyright 2010 Ivan Andrus. All rights reserved.

# Ensure we have enough arguments
if [ $# -lt 3 ]; then
    echo "usage: $0 SAGE_EXECUTABLE LOG [ARGS_FOR_NOTEBOOK]"
    exit 1;
fi

SAGE_EXECUTABLE="$1"
SAGE_LOG="$2"
NB_TYPE="$3"

# Read environment variables
cd $(dirname $SAGE_EXECUTABLE)
source local/bin/sage-env-config
source local/bin/sage-env


# Mac OS X app bundles are *intended* to be moved around, and/or given away
# So always run first the respective script handling this
# This should also catch Intel vs. PPC or 32Bit vs. 64Bit conflicts
echo Checking install location >> "$SAGE_LOG"
# TODO: If relocate-once.py is present run it with some sort of progress
# display, e.g. in a terminal


if [ "$NB_TYPE" = jupyter ]; then

./sage --notebook=jupyter $4 >> "$SAGE_LOG" 2>> "$SAGE_LOG" &


else

echo Checking existence of SageNB directory >> "$SAGE_LOG"
if [ -e $DOT_SAGE/sage_notebook.sagenb/users.pickle ]; then
    echo Starting Notebook >> "$SAGE_LOG"
    # $4 is not quoted because it comes as one argument from the app,
    # so we need the shell to parse it here.
    ./sage --notebook=sagenb $4 >> "$SAGE_LOG" 2>> "$SAGE_LOG"
else
    false
fi

# If it failed to start or it hasn't been run before, hope that we can
# fix it by running in a terminal to allow typing in a password.
if [ $? != 0 ]; then
    # if Terminal.app is not running before it is activated by
    # osascript, then it inherits the environment from osascript.
    # This includes SAGE_ENV_SOURCED which causes problems because
    # PATH can get messed up and we'll end up calling the system
    # python instead of the sage version.  We unset it here so that
    # sage-env will be sourced and PATH set up properly.
    SAGE_ENV_SOURCED=
    echo Starting Notebook in Terminal >> "$SAGE_LOG"
    sage-native-execute osascript \
        -e 'tell app "Terminal"' \
        -e '    activate' \
        -e "    do script \"'$SAGE_EXECUTABLE' --notebook=sagenb\"" \
        -e 'end'
    # We don't include $4 here since this should only happen the first time
    # they run it, and this way we don't have to worry about quoting it.
fi

fi

exit $?
