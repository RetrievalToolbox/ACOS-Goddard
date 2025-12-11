#!/bin/bash

## This is for re-setting the gh-pages branch and
## pushing a fresh history.

git branch --delete gh-pages # Delete local branch
git push github --delete gh-pages # Delete remote/github branch

git checkout --orphan gh-pages # Create new orphaned branch

git rm -rf . # Unstage all staged files

git checkout main -- docs # Grab the docs from `main`

git add docs # Stage docs files

# Make the commit
git commit -m "Documentation update / fresh gh-pages"

# Push docs up, but use docs/build as root
git subtree push --prefix docs/build github gh-pages
