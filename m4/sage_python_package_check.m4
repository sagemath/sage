#
# SYNOPSIS
#
#   SAGE_PYTHON_PACKAGE_CHECK(package)
#
# DESCRIPTION
#
#   Determine if the system copy of a python package can be used by sage.
#
#   This macro uses importlib.metadata and packaging to check that the
#   "version_requirements.txt" file (or entry in "pyproject.toml") for
#   the named package and its dependencies are satisfied. It can fail when:
#
#     1. If --enable-system-site-packages was not passed to ./configure,
#
#     2. If we are not using the system python (no $PYTHON_FOR_VENV),
#
#     3. If we are unable to create a venv with the system python,
#
#     4. If packaging is not available to the system python,
#
#     5. If the contents of version_requirements.txt (or entry in
#        "pyproject.toml") are not met (wrong version, no version,
#        etc.) by the system python.
#
#   In any of those cases, we set sage_spkg_install_$package to "yes"
#   so that the corresponding SPKG is installed. Otherwise, we do
#   nothing, since the default value of sage_spkg_install_$package
#   is "no" (to use the system copy).
#
#   The $PYTHON_FOR_VENV variable is set in configure.ac.
#   The check is run inside a new venv, and with the PYTHONUSERBASE
#   variable poisoned in the same manner as sage-env poisons it, to
#   ensure that the ./configure- and run-time views of the system
#   are as similar as possible.
#
#   To avoid suggesting these system packages to users who have not
#   set --enable-system-site-packages, this macro also changes the
#   default for --with-system-foo from "yes" to "no" in that case.
#

AC_DEFUN([SAGE_PYTHON_PACKAGE_CHECK], [
  AC_MSG_CHECKING([if --enable-system-site-packages was used])
  AS_IF([test "${enable_system_site_packages}" = "yes"], [
    AC_MSG_RESULT(yes)

    dnl We run this check inside a python venv, because that's ultimately
    dnl how the system $PYTHON_FOR_VENV will be used.
    AC_MSG_CHECKING([if we can create a python venv in config.venv])

    dnl Use --clear because ./configure typically clobbers its output files.
    AS_IF(["${PYTHON_FOR_VENV}" build/bin/sage-venv      dnl
                                  --system-site-packages dnl
                                  --clear                dnl
                                  config.venv            dnl
                                  2>&AS_MESSAGE_LOG_FD], [
      AC_MSG_RESULT(yes)
      dnl SAGE_PKG_VERSPEC is in the format of a toml list, but
      dnl without surrounding brackets, of single-quoted strings,
      dnl with any double-quotes escaped by backslash.
      AS_VAR_SET([SAGE_PKG_VERSPEC], ["SPKG_INSTALL_REQUIRES_]$1["])
      AC_MSG_CHECKING([for python package $1 (${SAGE_PKG_VERSPEC%,})])

      WITH_SAGE_PYTHONUSERBASE([dnl
        AS_IF(
          [config.venv/bin/python3 -c dnl
             "from importlib.metadata import distribution
from packaging.requirements import Requirement
from packaging.utils import canonicalize_name

pending = list((Requirement(text), ()) for text in (${SAGE_PKG_VERSPEC}))
visited = set()
while pending:
    requirement, parent_extras = pending.pop()
    if requirement.marker is not None and not any(
        requirement.marker.evaluate(dict(extra=extra))
        for extra in ('', *parent_extras)
    ):
        continue
    installed = distribution(requirement.name)
    # Like pkg_resources.require, accept installed prereleases.
    if not requirement.specifier.contains(installed.version, prereleases=True):
        raise RuntimeError(f'{requirement}: installed version is {installed.version}')
    extras = frozenset(map(canonicalize_name, requirement.extras))
    provided_extras = set(map(canonicalize_name,
                             installed.metadata.get_all('Provides-Extra') or ()))
    if not extras.issubset(provided_extras):
        raise RuntimeError(f'{requirement}: unknown extras {extras - provided_extras}')
    key = (canonicalize_name(requirement.name), extras)
    if key in visited:
        continue
    visited.add(key)
    pending.extend((Requirement(text), extras) for text in installed.requires or ())" dnl
           2>&AS_MESSAGE_LOG_FD],
          [AC_MSG_RESULT(yes)],
          [AC_MSG_RESULT(no); sage_spkg_install_$1=yes]
        )
      ])
    ], [
      dnl failed to create a venv for some reason
      AC_MSG_RESULT(no)
      sage_spkg_install_$1=yes
    ])

    dnl Clean up config.venv, but only if we could have created it.
    dnl (The --clear flag to pyvenv will not clobber a plain file.)
    AS_IF([test -d config.venv], [rm -rf config.venv])
  ], [
    dnl System site packages are disabled.
    AC_MSG_RESULT(no; skipping check)
    sage_spkg_install_$1=yes

    dnl We have to retroactively hack the --with-system-foo={no,yes,force}
    dnl mechanism here because it wasn't designed with the ability to
    dnl disable arbitrary chunks of system packages in mind. The easy cases
    dnl are "no" and "force" which require no action; "no" means we won't
    dnl suggest the package anyway, and "force" will raise an error when
    dnl the system-package check fails.
    dnl
    dnl The default of "yes" is more troubling because it is the default. To
    dnl avoid prompting users to install packages that won't be used, we want
    dnl to ignore "yes" when reporting the "hint: install these packages..."
    dnl at the end of ./configure. To accomplish that, we change "yes" to
    dnl "no" here, essentially changing the default for packages using this
    dnl macro when --enable-system-site-packages is disabled. Packages with
    dnl "no" are not suggested to the user.
    AS_IF([test "${sage_use_system_$1}" = "yes"],[sage_use_system_$1=no])
  ])
])


AC_DEFUN([WITH_SAGE_PYTHONUSERBASE], [dnl
  dnl To prevent user-site (pip install --user) packages from being
  dnl detected as "system" packages, we poison PYTHONUSERBASE. The
  dnl sage-env script also does this at runtime; we mimic that
  dnl implementation to ensure that the behaviors at ./configure and
  dnl runtime are identical. Beware that (as in sage-env) the poisoning
  dnl is skipped if PYTHONUSERBASE is non-empty. In particular, if the
  dnl user points PYTHONUSERBASE to any path (even the default), then
  dnl his local pip packages will be detected.
  PYTHONUSERBASE_SAVED="${PYTHONUSERBASE}"
  AS_IF([test -z "${PYTHONUSERBASE}"], [dnl
    PYTHONUSERBASE="${HOME}/.sage/local"
    export PYTHONUSERBASE
  ])
  $1
  PYTHONUSERBASE="${PYTHONUSERBASE_SAVED}"
])
