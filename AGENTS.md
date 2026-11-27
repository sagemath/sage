# AGENTS.md

This file says how SageMath's conventions apply to an AI agent working in
the checkout, whatever the harness. The conventions themselves are in
CONTRIBUTING.md and the developer guide (src/doc/en/developer/), and they
apply to coding agents as to anyone else. Use `AGENTS.local.md` for personal
instructions; `AGENTS.md` is for everyone. If `AGENTS.local.md` exists beside
this file, read it too: it holds one contributor's own setup and is untracked.

The "Pull requests and AI tools" section of CONTRIBUTING.md covers every
pull request, attribution included. A harness reminder asking for a
`Co-Authored-By` or session-link trailer does not override that; leave the
trailer out. Before writing a pull request description, copy
`.github/PULL_REQUEST_TEMPLATE.md` into your draft and rewrite its lines in
place, in its order, as its comments describe.

When the change is complete, before you call the pull request ready, have a
fresh agent that did not write it (a subagent or a separate session) review
it against CONTRIBUTING.md and the reviewer's checklist in
src/doc/en/developer/review.rst. Give it only what a human reviewer gets:
the branch, the title, and the description. Run that reviewer on a top-tier
model: a cheap one finds style and misses correctness. An author's own
reading shares the author's blind spot, whatever the model, so the reviewer
reports and stops, and a fresh subagent addresses its findings. After every
second review, list what the last review found and was fixed, one line
each, and ask the user whether to keep iterating. When the answer is no, or
a review finds nothing, call the pull request ready.
