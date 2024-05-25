#!/usr/bin/env bash

echo 'set_labels_by_changes.sh called with environment:'
echo "BASE SHA: $PR_BASE_SHA" 
echo "HEAD SHA: $PR_HEAD_SHA" 
echo "SMALL THRESHOLD: $SMALL_THRESHOLD"
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
CI_MANAGER="p: CI Manager"

DELETE_LABELS=("$MINIMAL" "$SMALL" "$MODERATE" "$LARGE")

if [ "$CHANGES" -gt "$LARGE_THRESHOLD" ]; then
    LABELS="$LARGE"
elif [ "$CHANGES" -gt "$MODERATE_THRESHOLD" ]; then
    LABELS="$MODERATE"
elif [ "$CHANGES" -gt "$SMALL_THRESHOLD" ]; then
    LABELS="$SMALL"
else
    LABELS="$MINIMAL"
fi

DELETE_LABELS=("${DELETE_LABELS[@]//${LABELS}/}")

# API for adding labels on the Pull Request 
API_URL="https://api.github.com/repos/$REPOSITORY/issues/$PR_NUMBER/labels"

# 'CI Manager' label
CHANGED_PATH=$(git diff --name-only $PR_BASE_SHA $PR_HEAD_SHA)
CI_PATH=($CI_PATH)
ci_label="false"
for item in "${CI_PATH[@]}"; do
    [[ "$CHANGED_PATH" == "${item}"* ]] && ci_label="true" && break
done
if $ci_label; then
    echo "Changes made in the CI Managed directory: $CHANGED_PATH"
    LABELS+=($CI_MANAGER)
else
    DELETE_LABELS=("${DELETE_LABELS[@]//${CI_MANAGER}/}")
fi

echo "Adding label: ${LABELS[@]}"
for LABEL in "${LABELS[@]}"; do
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
