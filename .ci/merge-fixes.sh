#!/bin/sh
# Merge open PRs from sagemath/sage labeled "blocker".
GH="gh -R sagemath/sage"
PRs="$($GH pr list --label "p: blocker / 1" --json number --jq '.[].number')"
if [ -z "$PRs" ]; then
    echo 'Nothing to do: Found no open PRs with "blocker" status.'
else
    echo "Found PRs: " $PRs
    export GIT_AUTHOR_NAME="ci-sage workflow"
    export GIT_AUTHOR_EMAIL="ci-sage@example.com"
    export GIT_COMMITTER_NAME="$GIT_AUTHOR_NAME"
    export GIT_COMMITTER_EMAIL="$GIT_AUTHOR_EMAIL"
    git commit -q -m "Uncommitted changes" --no-allow-empty -a
    git tag -f test_head
    for a in $PRs; do
        echo "::group::Merging PR #$a"
        $GH pr checkout $a
        git checkout -q test_head
        if git merge --no-edit -q FETCH_HEAD; then
            echo "::endgroup::"
            echo "Merged #$a"
        else
            echo "::endgroup::"
            echo "Failure merging #$a, resetting"
            git reset --hard
        fi
    done
fi
