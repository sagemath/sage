# If we are loaded with a workspace then \$SAGE will be defined and in
# that case we need to call StartInteract so that the pager will be
# set correctly.  See trac #5043.
if IsBound(\$SAGE) then
    \$SAGE.StartInteract();
fi;
