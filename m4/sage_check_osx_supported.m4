AC_DEFUN([SAGE_CHECK_OSX_SUPPORTED], [
    AC_REQUIRE([AC_CANONICAL_HOST])
    AS_CASE([$host],
    [*-apple-darwin*], [
        # Warning: xcodebuild does not seem to be maintained in Xcode 4.3
        # or later, so do not rely on the variable XCODE_VERS with OS X
        # 10.7 or later.
        changequote(<,>)
        XCODE_VERS=`xcodebuild -version 2> /dev/null | grep Xcode | sed -e 's/[A-Za-z ]//g'`
        changequote([,])
        if test -z $XCODE_VERS; then
            XCODE_VERS="2"
        fi
        XCODE_VERS_MAJOR=`echo $XCODE_VERS | cut '-d.' -f1`
        DARWIN_VERSION=`uname -r | cut '-d.' -f1`
        if test $DARWIN_VERSION -gt 10; then
            echo "You are using OS X Lion (or later)."
            echo "You are strongly advised to install Apple's latest Xcode"
            echo "unless you already have it. You can install this using"
            echo "the App Store. Also, make sure you install Xcode's"
            echo "Command Line Tools -- see Sage's README.txt."
        elif test $XCODE_VERS_MAJOR -gt 2; then
            echo "You are using Xcode version $XCODE_VERS."
            echo "You are strongly advised to install Apple's latest Xcode"
            echo "unless you already have it. You can download this from"
            echo "http://developer.apple.com/downloads/."
            echo "If using Xcode 4.3 or later, make sure you install Xcode's"
            echo "Command Line Tools -- see Sage's README.txt."
        else
            echo "You are using Xcode version 1 or 2"
            echo "WARNING: You are strongly advised to install the"
            echo "latest version of Apple's Xcode for your platform,"
            echo "unless you already have it."
            if test $DARWIN_VERSION -eq 10; then
               echo "Probably you need Xcode 3.2.6"
            elif test $DARWIN_VERSION -eq 9; then
               echo "Probably you need Xcode 3.1.4"
            elif test $DARWIN_VERSION -lt 9; then
               echo "Probably you need Xcode 2.5"
            fi
        fi >& AS_MESSAGE_FD

    ])
])
