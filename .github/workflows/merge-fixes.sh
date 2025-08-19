#!/bin/sh
# Apply open PRs labeled "p: CI Fix" from sagemath/sage as patches.
# (policy set by vote in 2024-03,
#  https://groups.google.com/g/sage-devel/c/OKwwUGyKveo/m/vpyCXYBqAAAJ)
#
# This script is invoked by various workflows in .github/workflows
#
# The repository variable SAGE_CI_FIXES_FROM_REPOS can be set
# to customize this for CI runs in forks:
#
# - If set to a whitespace-separated list of repositories, use them instead of sagemath/sage.
# - If set to "none", do not apply any PRs.
#
# https://docs.github.com/en/actions/learn-github-actions/variables#creating-configuration-variables-for-a-repository
export GIT_AUTHOR_NAME="ci-sage workflow"
export GIT_AUTHOR_EMAIL="ci-sage@example.com"
export GIT_COMMITTER_NAME="$GIT_AUTHOR_NAME"
export GIT_COMMITTER_EMAIL="$GIT_AUTHOR_EMAIL"
mkdir -p upstream
for REPO in ${SAGE_CI_FIXES_FROM_REPOSITORIES:-sagemath/sage}; do
    case $REPO in
        none)
            echo "Nothing to do for 'none' in SAGE_CI_FIXES_FROM_REPOSITORIES"
            ;;
        */*)
            echo "Getting open PRs with 'p: CI Fix' label from https://github.com/$REPO/pulls?q=is%3Aopen+label%3A%22p%3A+CI+Fix%22"
            GH="gh -R $REPO"
            REPO_FILE="upstream/ci-fixes-${REPO%%/*}-${REPO##*/}"
            PRs="$($GH pr list --label "p: CI Fix" --json number --jq '.[].number' | tee $REPO_FILE)"
            date -u +"%Y-%m-%dT%H:%M:%SZ" > $REPO_FILE.date  # Record the date, for future reference
            if [ -z "$PRs" ]; then
                echo "Nothing to do: Found no open PRs with 'p: CI Fix' label in $REPO."
            else
                echo "Found open PRs with 'p: CI Fix' label in $REPO: $(echo $PRs)"
                git tag -f test_base
                git commit -q -m "Uncommitted changes" --no-allow-empty -a
                for a in $PRs; do
                    # We used to pull the branch and merge it (after unshallowing), but when run on PRs that are
                    # based on older releases, it has the side effect of updating to this release,
                    # which may be confusing.
                    #
                    # Here, we obtain the "git am"-able patch of the PR branch.
                    # This also makes it unnecessary to unshallow the repository.
                    #
                    # Considered alternative: Use https://github.com/$REPO/pull/$a.diff,
                    # which squashes everything into one diff without commit metadata.
                    PULL_URL="https://github.com/$REPO/pull/$a"
                    PULL_SHORT="$REPO#$a"
                    PULL_FILE="$REPO_FILE-$a"
                    PATH=build/bin:$PATH build/bin/sage-download-file --quiet "$PULL_URL.patch" $PULL_FILE.patch
                    date -u +"%Y-%m-%dT%H:%M:%SZ" > $PULL_FILE.date  # Record the date, for future reference
                    LAST_SHA=$(sed -n -E '/^From [0-9a-f]{40}/s/^From ([0-9a-f]{40}).*/\1/p' $PULL_FILE.patch | tail -n 1)
                    COMMITS_URL="https://github.com/$REPO/commits/$LAST_SHA"
                    echo "::group::Applying PR $PULL_URL @ $COMMITS_URL as a patch"
                    export GIT_COMMITTER_NAME="$GIT_AUTHOR_NAME applying $PULL_URL @ $COMMITS_URL"
                    if git am --signoff --empty=keep < $PULL_FILE.patch; then
                        echo "---- Applied patch --------------------------------------------------------------------------------"
                        cat $PULL_FILE.patch
                        echo "--------------------------------------------------------------------8<-----------------------------"
                        echo "::endgroup::"
                    elif git am --abort \
                            && if git fetch --unshallow --all > /dev/null 2>&1; then echo "Unshallowed"; fi \
                            && echo "Retrying with 3-way merge" \
                            && git am --empty=keep --3way < $PULL_FILE.patch; then
                        echo "---- Applied patch --------------------------------------------------------------------------------"
                        cat $PULL_FILE.patch
                        echo "--------------------------------------------------------------------8<-----------------------------"
                        echo "::endgroup::"
                    else
                        echo "---- Failure applying patch -----------------------------------------------------------------------"
                        git am --signoff --show-current-patch=diff
                        echo "--------------------------------------------------------------------8<-----------------------------"
                        echo "::endgroup::"
                        echo "Failure applying $PULL_SHORT as a patch, resetting"
                        git am --signoff --abort
                    fi
                done
            fi
            ;;
        *)
            echo "Repository variable SAGE_CI_FIXES_FROM_REPOSITORIES, if set, should be a whitespace-separated list of repositories or 'none'"
            echo "Got: $SAGE_CI_FIXES_FROM_REPOSITORIES"
            exit 1
            ;;
    esac
done
