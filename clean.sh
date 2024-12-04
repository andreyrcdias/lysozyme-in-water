#!/bin/bash
set -x

echo "Cleaning target files..."
find . -maxdepth 1 ! -name 'clean.sh' ! -name 'run.sh' ! -name 'README.md' -exec rm -f {} +
