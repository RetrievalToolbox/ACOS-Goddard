#!/bin/bash

## This is for updating the gh-pages branch based on
## updates on the main branch.

git switch gh-pages

git checkout main -- docs # Grab the docs from `main`

julia --project=./docs ./docs/make.jl

git add docs # Stage docs files

# Make the commit
commit_hash=`git log -n 1 main --pretty=format:"%H"`
git commit -m "Documentation update from commit ${commit_hash}"

# Push docs up, but use docs/build as root
git subtree push --prefix docs/build github gh-pages
