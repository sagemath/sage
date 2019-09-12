# Continuous Integration (CI)

We support several implementations of CI. All these implementations rely on
[docker](https://docker.com) in some way. This directory contains bits which
are shared between these CI implementations. The relevant docker files can be
found in `/docker/`.

* [CircleCI](https://circleci.com) is configured in `/.circleci/`.
* [GitLab CI](https://gitlab.com) is configured in `/.gitlab-ci.yml`.
