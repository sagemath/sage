#!/usr/bin/env bash
set -euo pipefail

: "${CI_SOURCE:?}" "${CI_RESULTS:?}" "${CI_RUN_DIR:?}" "${CI_HELPERS:?}"
export PATH="/opt/homebrew/bin:/opt/homebrew/sbin:$PATH"
export DOT_SAGE="$CI_RUN_DIR/dot-sage"
export SAGE_DATA_PATH="$CI_RUN_DIR/data"
export IPYTHONDIR="$CI_RUN_DIR/ipython"
export MPLCONFIGDIR="$CI_RUN_DIR/matplotlib"
export XDG_CACHE_HOME="$CI_RUN_DIR/cache"
export CCACHE_DIR="${RUNNER_TEMP:-$CI_RUN_DIR}/sage-homebrew-ccache"
export CCACHE_MAXSIZE=1G
export PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME SAGE_LOCAL SAGE_ROOT SAGE_VENV
cd "$CI_SOURCE"

activate() {
    # shellcheck disable=SC1091
    source "$CI_RUN_DIR/venv/bin/activate"
    # Includes the native compiler/linker paths from sagemath-deps.
    # shellcheck disable=SC1090,SC1091
    source "$(brew --prefix dimpase/tap/sagemath-docbuild-deps)/share/sagemath-docbuild-deps/python3.14-env.sh"
}

case "${1:?Expected dependencies, venv, build, verify, pytest, or doctest}" in
    dependencies)
        brew tap dimpase/tap
        brew tap dimpase/gap
        brew tap macaulay2/tap
        brew trust --tap dimpase/tap dimpase/gap macaulay2/tap
        # Fetch the tap metadata for the published dependency bottles.
        brew update
        # These formulae declare conflicts with their core counterparts.
        # Keep the core kegs installed for other jobs on the runner.
        for formula in numpy scipy python-matplotlib; do
            if [[ -L "$(brew --prefix)/var/homebrew/linked/$formula" ]]; then
                brew unlink "$formula"
            fi
        done
        brew install --force-bottle dimpase/tap/sagemath-deps dimpase/tap/sagemath-docbuild-deps
        brew list --versions > "$CI_RESULTS/brew-versions.txt"
        brew info --json=v2 dimpase/tap/sagemath-deps dimpase/tap/sagemath-docbuild-deps \
            > "$CI_RESULTS/dependency-bottles.json"
        for tap in dimpase/tap dimpase/gap macaulay2/tap; do
            printf '%s ' "$tap"
            git -C "$(brew --repo "$tap")" rev-parse HEAD
        done | tee "$CI_RESULTS/tap-commits.txt"
        sw_vers
        df -h .
        ;;
    venv)
        "$(brew --prefix python@3.14)/bin/python3.14" -m venv --system-site-packages "$CI_RUN_DIR/venv"
        activate
        # Freeze the bottled Python dependencies so adding pytest cannot upgrade them.
        python -m pip freeze --all --exclude sagemath > "$CI_RESULTS/bottled-python.txt"
        python "$CI_HELPERS/check.py" requirements > "$CI_RESULTS/test-requirements.txt"
        python -m pip install --only-binary=:all: -c "$CI_RESULTS/bottled-python.txt" \
            -r "$CI_RESULTS/test-requirements.txt"
        python "$CI_HELPERS/check.py" data | tee "$CI_RESULTS/data-paths.json"
        python -m pip list --format=json > "$CI_RESULTS/python-packages.json"
        ;;
    build)
        activate
        setup_args=(-Dbuild-docs=false -Ddefer_feature_checks=true --wrap-mode=nofallback)
        for feature in bliss brial coxeter3 eclib libbraiding libhomfly mcqd meataxe rankwidth sirocco tdlib; do
            setup_args+=("-D$feature=enabled")
        done
        config_args=()
        for arg in "${setup_args[@]}"; do
            config_args+=("--config-settings=setup-args=$arg")
        done
        python -m pip install --no-build-isolation --no-deps --editable . -v \
            "--config-settings=build-dir=$CI_RUN_DIR/build" \
            "--config-settings=compile-args=-j${CI_JOBS:-4}" "${config_args[@]}"
        ;;
    verify)
        activate
        python "$CI_HELPERS/check.py" provenance | tee "$CI_RESULTS/provenance.json"
        python -m pip check
        ;;
    pytest)
        activate
        python -m pytest -n "${CI_JOBS:-4}" -rfEs -s \
            --junitxml="$CI_RESULTS/pytest.xml" src/sage
        ;;
    doctest)
        activate
        python "$CI_HELPERS/check.py" doctest
        ;;
    *)
        echo "Unknown phase: $1" >&2
        exit 2
        ;;
esac
