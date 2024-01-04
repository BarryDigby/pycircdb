#!/usr/bin/env bash

find . -type f -size +20M | sed 's/^\.\.\?\///g' > .gitignore
echo "results/old_res/" >> .gitignore
