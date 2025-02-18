#!/bin/sh
if [ $# != 2 ]; then
    echo >&2 "usage: $0 WORKTREE_NAME WORKTREE_DIRECTORY"
    echo >&2 "Ensures that the current working directory is a git repository,"
    echo >&2 "then makes WORKTREE_DIRECTORY a git worktree named WORKTREE_NAME."
fi
WORKTREE_NAME="$1"
WORKTREE_DIRECTORY="$2"

export GIT_AUTHOR_NAME="ci-sage workflow"
export GIT_AUTHOR_EMAIL="ci-sage@example.com"
export GIT_COMMITTER_NAME="$GIT_AUTHOR_NAME"
export GIT_COMMITTER_EMAIL="$GIT_AUTHOR_EMAIL"

set -e

# Set globally for other parts of the workflow
git config --global user.name "$GIT_AUTHOR_NAME"
git config --global user.email "$GIT_AUTHOR_EMAIL"

set -x

# If actions/checkout downloaded our source tree (in current directory) using the
# GitHub REST API instead of with git, we first make the source tree a repo.
if [ ! -d .git ]; then git init; fi

# Commit and tag this state of the source tree "new". The current directory
# contains "new" source files with which we want to build and test.
git add -A && git commit --quiet --allow-empty -m "new"
git tag -f new

# Our container image contains a source tree in $WORKTREE_DIRECTORY with a full build of Sage.
# But $WORKTREE_DIRECTORY is not a git repository.
# We make $WORKTREE_DIRECTORY a worktree whose index is at tag "new".
if [ -L $WORKTREE_NAME ]; then
    rm -f $WORKTREE_NAME
fi
git worktree prune --verbose
git worktree add --detach $WORKTREE_NAME
rm -rf $WORKTREE_DIRECTORY/.git && mv $WORKTREE_NAME/.git $WORKTREE_DIRECTORY/
rm -rf $WORKTREE_NAME && ln -s $WORKTREE_DIRECTORY $WORKTREE_NAME

# Copy "new" .gitignore (in current directory) to the worktree if there is not one already.
if [ ! -f $WORKTREE_NAME/.gitignore ]; then cp .gitignore $WORKTREE_NAME/; fi
# We then commit the worktree and set the tag "old". (This keeps all modification times unchanged.)
(cd $WORKTREE_NAME && git add -A && git commit --quiet --allow-empty -m "old" && git tag -f old)
# Then we update the worktree and index with "git checkout new". (This keeps modification times of "old" files
# and makes modification times of "new" files newer, which triggers incremental build.)
(cd $WORKTREE_NAME && git checkout -f new && \
 # For unknown reason, "build/" directory is untracked by git. Hence we force
 # it tracked, to avoid removal of "build/make/Makefile" by "git clean -fd"
 git add build/ && git clean -fd)
