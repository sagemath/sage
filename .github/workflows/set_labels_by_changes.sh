#!/usr/bin/env bash

echo 'set_labels_by_changes.sh called with environment:'
echo "BASE SHA: $PR_BASE_SHA" 
echo "HEAD SHA: $PR_HEAD_SHA" 
echo "SMALL THRESHOLD $SMALL_THRESHOLD"
echo "MODERATE THERESHOLD: $MODERATE_THRESHOLD"
echo "LARGE THRESHOLD: $LARGE_THRESHOLD"

# get all the changes made and changed files
CHANGES=$(git diff --ignore-all-space $PR_BASE_SHA $PR_HEAD_SHA)

# ignore blank lines
CHANGES=$(echo "$CHANGES" | grep -vE '^[\+\-]\s*$')

# ignore non necessary lines from git diff
CHANGES=$(echo "$CHANGES" | grep -E '^[+\-]' | grep -vE '^\+\+\+|^\-\-\-')

# count total no of lines
CHANGES=$(echo "$CHANGES" | wc -l)

echo "CHANGES MADE: $CHANGES"

AUTH_HEADER="Authorization: Bearer $GITHUB_TOKEN"

MINIMAL="v: minimal"
SMALL="v: small"
MODERATE="v: moderate"
LARGE="v: large"

DELETE_LABELS=("$MINIMAL" "$SMALL" "$MODERATE" "$LARGE")

if [ "$CHANGES" -gt "$LARGE_THRESHOLD" ]; then
    SIZE_LABEL="$LARGE"
elif [ "$CHANGES" -gt "$MODERATE_THRESHOLD" ]; then
    SIZE_LABEL="$MODERATE"
elif [ "$CHANGES" -gt "$SMALL_THRESHOLD" ]; then
    SIZE_LABEL="$SMALL"
else
    SIZE_LABEL="$MINIMAL"
fi

DELETE_LABELS=("${DELETE_LABELS[@]//${SIZE_LABEL}/}")

# API for adding labels on the Pull Request 
API_URL="https://api.github.com/repos/$REPOSITORY/issues/$PR_NUMBER/labels"

echo "Adding label: ${SIZE_LABEL[@]}"
for LABEL in "${SIZE_LABEL[@]}"; do
    curl -X POST \
        -H "$AUTH_HEADER" \
        -H "Accept: application/vnd.github+json" \
        -H "X-GitHub-Api-Version: 2022-11-28" \
        -d "{\"labels\":[\"$LABEL\"]}" \
        "$API_URL" >/dev/null
done

echo "Deleting Labels:"

for DELETE_LABEL in "${DELETE_LABELS[@]}"; do
    ENCODED_LABEL=$(echo "$DELETE_LABEL" | sed 's/ /%20/g')
    curl -X DELETE \
        -H "Accept: application/vnd.github+json" \
        -H "$AUTH_HEADER" \
        -H "X-GitHub-Api-Version: 2022-11-28" \
        "$API_URL/$ENCODED_LABEL" >/dev/null
done
