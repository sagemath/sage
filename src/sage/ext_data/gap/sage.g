#
# Sage support utilities to read into the GAP session.
#

# Disable color in the prompt, it interferes with parsing.
ColorPrompt(false);

# Prevent loading the xgap package; we use the -p flag to GAP in order to
# communicate with it via the pexpect interface; this is normally used by
# for an xgap window to communicate with GAP, so unfortunately setting this
# flag also allows the xgap package to be loaded and for some packages to
# attempt to communicate with a "window handler" that does not exist.
# Therefore we must explicitly disable loading of the xgap package.
#
# Do not use SetUserPreference since that leads to reloading the workspace,
# which is confusing to the pexpect interface
if IsBound(GAPInfo.ExcludeFromAutoload) then
    Append(GAPInfo.ExcludeFromAutoload, "xgap");
else
    GAPInfo.ExcludeFromAutoload := [ "xgap" ];
fi;

\$SAGE := rec();

\$SAGE.OldPager := Pager;


\$SAGE.NewPager :=
          function( data )
    local   str,  lines,  line, fn, start;
    str := OutputTextFile(\$SAGE.tempfile,false);
    start := 1;
    if IsRecord(data) then
        lines := data.lines;
        if IsBound(data.start) then
            start := data.start;
        fi;
    else
        lines := data;
    fi;
    if IsString(lines) then
        lines := SplitString(lines,"\n");
    fi;
    for line in lines do
        WriteLine(str, line);
    od;
    CloseStream(str);
    Print("Page from ",start,"\n");
end;

\$SAGE.StartInteract := function()
    MakeReadWriteGlobal("Pager");
    Pager := \$SAGE.OldPager;
    HELP_VIEWER_INFO.screen.show := \$SAGE.OldPager;
    MakeReadOnlyGlobal("Pager");
end;


\$SAGE.StopInteract := function()
    MakeReadWriteGlobal("Pager");
    Pager := \$SAGE.NewPager;
    HELP_VIEWER_INFO.screen.show := \$SAGE.NewPager;
    MakeReadOnlyGlobal("Pager");
end;


\$SAGE.StopInteract();

#\$SAGE.ErrorHandler := function(m,a,m2,mode)
#    PrintTo("*errout*", m);
#    if a <> fail then
#        PrintTo("*errout*",a);
#    fi;
#    SetErrorHandler(\$SAGE.ErrorHandler);
#    return true;
#end;

#SetErrorHandler(\$SAGE.ErrorHandler);

SetAllInfoLevels(0);

\$SAGE.OperationsAdmittingFirstArgument := function(obj)
    local   myflags, mfi;
    myflags := FlagsObj(obj);
    mfi := function(o)
        local f;
        f := GET_OPER_FLAGS(o);
        return f<>[] and f[1]<>[] and
          IS_SUBSET_FLAGS(myflags, f[1][1]);
    end;
    return Filtered(OPERATIONS, mfi);
end;


\$SAGE.CleanOperationName := function(name)
    local   lt,  ls;
    lt := Length("Tester(");
    if Length(name) > lt and name{[1..lt]} = "Tester(" then
        return Concatenation("Has",name{[lt+1..Length(name)-1]});
    fi;
    ls := Length("Setter(");
    if Length(name) > ls and name{[1..ls]} = "Setter(" then
        return Concatenation("Set",name{[lt+1..Length(name)-1]});
    fi;
    return name;
end;

\$SAGE.HasAtLeastOneMethodAsFirstArgument := function(op,obj)
    local   t,  f,  n,  meths,  i;
    t := TypeObj(obj);
    f := FlagsType(t);
    for n in [1..6] do
        meths := METHODS_OPERATION(op,n);
        for i in [1,n+5..LENGTH(meths)-n-3] do
            if IS_SUBSET_FLAGS(f,meths[i+1]) then
                return true;
            fi;
        od;
    od;
    return false;
end;


\$SAGE.PlausibleTabCompletionsForSage := function(o)
    local   ops,  opnames;
    ops := Filtered(\$SAGE.OperationsAdmittingFirstArgument(o), op ->
                   \$SAGE.HasAtLeastOneMethodAsFirstArgument(op,o));
    opnames := List(ops, op -> \$SAGE.CleanOperationName(NameFunction(op)));
    return Concatenation(opnames, GLOBAL_FUNCTION_NAMES);
end;

# The log below is for debugging only.
# CAREFUL -- do *not* activate this unless you know
# what you are doing.  E.g., if active and the user doesn't
# have write permission to /tmp (e.g., on OS X),
# then gap will completely fail to work for them.   -- WAS
#
# LogTo("/tmp/gapsage.log");
#


# Load the GAP packages that GAP itself tries to autoload in the
# default configuration (see "PackagesToLoad" in lib/package.gi). The
# combination of passing -A to gap and these LoadPackage statements
# allows us to load the usual set of packages, but only if they are
# installed. So most people will get exactly the default behavior,
# but minimal installations won't throw warnings and fail tests.
_autoloads := [ "autpgrp", "alnuth", "crisp", "ctbllib", "factint", "fga",
                "irredsol", "laguna", "polenta", "polycyclic", "resclasses",
                "sophus", "tomlib" ];
for p in _autoloads do
  LoadPackage(p);
od;
