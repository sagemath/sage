@echo off
rem Authors:
rem * Dmitrii Pasechnik <dimpase@gmail.com>
rem * Jean-Pierre Flori <jean-pierre.flori@ssi.gouv.fr>
rem * Erik M. Bray <erik.bray@lri.fr>
rem
rem Rebase all dlls in the SAGE_LOCAL directory (and its subdirectories),
rem as well as the ones already stored in the system database,
rem and update the database.
rem This system-wide database is located in '/etc/rebase.db.i386' and
rem includes the Cygwin system dlls.
rem
rem Invoke this script from a Windows command prompt
rem and, if Cygwin is installed in a non-standard location,
rem adjusting CYGWIN_ROOT.
rem
rem Ensure that no other Cygwin processes are currently running.
rem Note that you need write access to the system-wide rebase database
rem (which usually means admin rights).

set THIS_BIN=%~dp0

rem SAGE_LOCAL should be one level up from the bin/ this script is in
rem This is about the most elegant way to do this I can find thanks
rem http://stackoverflow.com/a/33404867/982257
call :NORMALIZEPATH "%THIS_BIN%.."
set SAGE_LOCAL=%RETVAL%

rem Cygwin saves its installation root here
rem If there is more than one Cygwin installation on the system this
rem will just pick up the first one
call :READREGISTRY HKEY_LOCAL_MACHINE\SOFTWARE\Cygwin\setup rootdir
set CYGWIN_ROOT=%RETVAL%

rem Make sure dash can be called from MSDOS prompt:
path %CYGWIN_ROOT%\bin;%path%
rem Suppress warning about MSDOS-style path:
set CYGWIN=%CYGWIN% nodosfilewarning
rem Call the dash script to do the real work:
cd %SAGE_LOCAL%
dash bin\sage-rebaseall.sh


:: ========== FUNCTIONS ==========
exit /B

:READREGISTRY
  for /F "usebackq skip=2 tokens=3" %%V in (`reg query %1 /v %2 2^>nul`) do (
    set RETVAL=%%V
    break
  )
  exit /B

:NORMALIZEPATH
  set RETVAL=%~dpfn1
  exit /B
