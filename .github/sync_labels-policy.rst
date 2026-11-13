Label synchronization and Actions policies
==========================================

The ``sync_labels.yml`` workflow needs ``pull_request_target`` to update
labels and reviews on pull requests from forks. It executes only
``.github/sync_labels.py`` from the repository's default branch. Pull request
code is never checked out or executed; event fields and GitHub API responses
are treated as data. The checkout must retain its explicit default-branch
``ref``, including for the separate ``pull_request_review`` trigger.

Administrator setup
-------------------

GitHub will enforce its default restriction on ``pull_request_target`` in
affected public repositories on December 2, 2026. Merging a workflow change
does not configure the repository's Actions policy.

Before enforcement, a repository administrator or organization owner should:

1. Review **Settings > Actions > Policy insights** to identify affected runs
   and any applicable organization policies.
2. Under **Settings > Actions > Policies**, create or update an event policy
   scoped to ``.github/workflows/sync_labels.yml`` in this repository. Allow
   its four event types: ``issues``, ``pull_request_target``,
   ``pull_request_review``, and ``schedule``. Preserve the ability to respond
   to contributions from forks when configuring any actor restrictions.
3. Verify in policy insights that the workflow's intended runs are permitted.
   A repository policy cannot override a restriction imposed by an applicable
   organization or enterprise policy; involve its owner if necessary.

Do not enable ``allow-unsafe-pr-checkout`` or broaden this exception to
workflows that build or test pull request code.

The workflow limits ``GITHUB_TOKEN`` to reading contents and writing issues
and pull requests. If ``SYNC_LABELS_BOT_TOKEN`` is configured, its permissions
must be reviewed separately: workflow ``permissions`` does not restrict a
personal access token. Preserve the bot's required organization membership
access for triage checks. This change does not alter the existing limitations
on secrets and write tokens for fork ``pull_request_review`` runs.

References:

* `Securely using pull_request_target
  <https://docs.github.com/en/actions/reference/security/securely-using-pull_request_target#default-policy-for-pull_request_target>`_
* `Governing execution of GitHub Actions workflows
  <https://docs.github.com/en/actions/how-tos/administer/control-workflow-execution>`_

Offline validation
------------------

Run ``python3 -m unittest build.test.test_sync_labels``. These tests use a
fake ``gh`` executable and mocked handlers; they do not contact GitHub or
change labels, comments, reviews, or repository policies.
The ordinary lint workflow also runs these tests without a bot token.
