#!/bin/bash
if [ x"$1" = x"sage-jupyter" ]; then
    # If "sage-jupyter" is given as a first argument, we start a jupyter notebook
    # with reasonable default parameters for running it inside a container.
    shift
    exec sage -n jupyter --no-browser --ip='0.0.0.0' --port=8888 "$@"
else
    exec sage -sh -c "$*"
fi
