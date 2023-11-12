#!/bin/sh
# Merge open PRs from sagemath/sage labeled "blocker".
REPO="sagemath/sage"
GH="gh -R $REPO"
mkdir -p upstream
PRs="$($GH pr list --label "p: blocker / 1" --json number --jq '.[].number' | tee upstream/ci-fixes.txt)"
date -u +"%Y-%m-%dT%H:%M:%SZ" > upstream/ci-fixes.date  # Record the date, for future reference
if [ -z "$PRs" ]; then
    echo 'Nothing to do: Found no open PRs with "blocker" status.'
else
    echo "Found PRs: " $PRs
    export GIT_AUTHOR_NAME="ci-sage workflow"
    export GIT_AUTHOR_EMAIL="ci-sage@example.com"
    export GIT_COMMITTER_NAME="$GIT_AUTHOR_NAME"
    export GIT_COMMITTER_EMAIL="$GIT_AUTHOR_EMAIL"
    git tag -f test_base
    git commit -q -m "Uncommitted changes" --no-allow-empty -a
    for a in $PRs; do
        git tag -f test_head
        # We used to pull the branch and merge it (after unshallowing), but when run on PRs that are
        # based on older releases, it has the side effect of updating to this release,
        # which may be confusing.
        #
        # Here, we obtain the "git am"-able patch of the PR branch.
        # This also makes it unnecessary to unshallow the repository.
        #
        # Considered alternative: Use https://github.com/$REPO/pull/$a.diff,
        # which squashes everything into one diff without commit metadata.
        PATH=build/bin:$PATH build/bin/sage-download-file "https://github.com/$REPO/pull/$a.patch" upstream/$a.patch
        date -u +"%Y-%m-%dT%H:%M:%SZ" > upstream/$a.date  # Record the date, for future reference
        LAST_SHA=$(sed -n -E '/^From [0-9a-f]{40}/s/^From ([0-9a-f]{40}).*/\1/p' upstream/$a.patch | tail -n 1)
        echo "::group::Applying PR https://github.com/$REPO/pull/$a @ https://github.com/$REPO/commit/$LAST_SHA as a patch"
        export GIT_COMMITTER_NAME="$GIT_AUTHOR_NAME applying https://github.com/$REPO/pull/$a @ https://github.com/$REPO/commit/$LAST_SHA"
        if git am --signoff --empty=keep < upstream/$a.patch; then
            echo "---- Applied patch ------------------------------------------------------------"
            cat upstream/$a.patch
            echo "--------------------------------------------------------------------8<---------"
            echo "::endgroup::"
        elif git am --abort \
                && if git fetch --unshallow --all > /dev/null 2>&1; then echo "Unshallowed"; fi \
                && echo "Retrying with 3-way merge" \
                && git am --empty=keep --3way < upstream/$a.patch; then
            echo "---- Applied patch ------------------------------------------------------------"
            cat upstream/$a.patch
            echo "--------------------------------------------------------------------8<---------"
            echo "::endgroup::"
        else
            echo "---- Failure applying patch ---------------------------------------------------"
            git am --signoff --show-current-patch=diff
            echo "--------------------------------------------------------------------8<---------"
            echo "::endgroup::"
            echo "Failure applying #$a as a patch, resetting"
            git am --signoff --abort
        fi
    done
    #git log test_base..HEAD
fi
