#!/bin/bash
# Necessary to push files while working in UKB RAP
dx download -r /.ssh
eval "$(ssh-agent -s)"
chmod 600 .ssh/id_ed25519
ssh-add .ssh/id_ed25519
ssh-keyscan github.com >> .ssh/known_hosts
git push
