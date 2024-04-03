#!/usr/bin/env bash

# get all the changes made and changed files
CHANGES=$(git diff --ignore-all-space  $PR_BASE_SHA $PR_HEAD_SHA)

# ignore blank lines
CHANGES=$(echo "$CHANGES" | grep -vE '^[\+\-]\s*$')

# ignore non necessary lines from git diff
CHANGES=$(echo "$CHANGES" | grep -E '^[+\-]' | grep -vE '^\+\+\+|^\-\-\-')

# count total no of lines
CHANGES=$(echo "$CHANGES" | wc -l)

echo "CHANGES MADE: $CHANGES"

# get changed files
CHANGED_PATH=$(git diff --name-only $PR_BASE_SHA $PR_HEAD_SHA)

# extract changed dir
CHANGED_DIR=$(echo "$CHANGED_PATH" | awk -F'src/sage/' '{print $2}' | cut -d'/' -f1 | sed 's/\([^ ]\+\)/c: \1/g')

AUTH_HEADER="Authorization: Bearer $GITHUB_TOKEN"

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

DELETE_LABELS=("${DELETE_LABELS[@]//${SIZE_LABEL}}")

readarray -t LABELS <<< "$CHANGED_DIR"
LABELS+=("$SIZE_LABEL")

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