#!/usr/bin/env bash

# Exit on error
set -e

# Setup trac as remote
## In order to push to trac, generate a new key with `ssh-keygen -f tempkey` and save the private key to gitpod `gp env PRIVATE_SSH_KEY="$(<tempkey)"` (or by following https://www.gitpod.io/docs/environment-variables#using-the-account-settings)
## then follow https://doc.sagemath.org/html/en/developer/trac.html#linking-your-public-key-to-your-trac-account to register the public key with trac.
## Afterwards, create a new gitpod workspace.
git remote remove trac 2> /dev/null || true # might still exists from a previous run/prebuild
if [[ -n "${PRIVATE_SSH_KEY}" ]]; then
  # Setup ssh key for authentication with trac
  mkdir -p ~/.ssh
  echo $PRIVATE_SSH_KEY | sed 's/\(-----\(BEGIN\|END\) OPENSSH PRIVATE KEY-----\)/\n\1\n/g' > ~/.ssh/id_rsa
  sed -i '/^$/d' ~/.ssh/id_rsa
  chmod 600 ~/.ssh/id_rsa
  echo "PubkeyAcceptedKeyTypes +ssh-rsa" > ~/.ssh/config
  ssh-keyscan -H trac.sagemath.org >> ~/.ssh/known_hosts

  # Setup trac repo
  git remote add trac git@trac.sagemath.org:sage.git -t master -t develop -t $(git branch --show-current)
  git remote set-url --push trac git@trac.sagemath.org:sage.git
  git fetch trac
  git branch -u trac/$(git branch --show-current)
else
  # Fallback to sagemath mirror
  git remote add trac https://github.com/sagemath/sagetrac-mirror.git -t master -t develop
  git remote set-url --push trac pushing-needs-ssh-key
fi
