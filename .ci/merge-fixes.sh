#!/bin/sh
# Merge open PRs from sagemath/sage labeled "blocker".
REPO="sagemath/sage"
GH="gh -R $REPO"
PRs="$($GH pr list --label "p: blocker / 1" --json number --jq '.[].number')"
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
        echo "::group::Merging PR https://github.com/$REPO/pull/$a"
        git tag -f test_head
        $GH pr checkout -b pr-$a $a
        git fetch --unshallow --all
        git checkout -q test_head
        if git merge --no-edit --squash -q pr-$a; then
            echo "::endgroup::"
            if git commit -q -m "Merge https://github.com/$REPO/pull/$a" -a --no-allow-empty; then
                echo "Merged #$a"
            else
                echo "Empty, skipped"
            fi
        else
            echo "::endgroup::"
            echo "Failure merging #$a, resetting"
            git reset --hard
        fi
    done
    git log test_base..HEAD
fi
