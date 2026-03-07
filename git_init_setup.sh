#!/bin/bash

git stash && git pull

# Init Setup - Convert to Unix line endings & set executable permissions

cd "$(dirname "${BASH_SOURCE[0]}")"

# Convert all text files to Unix line endings
find . -type f \( -name "*.sh" -o -name "*.py" -o -name "*.txt" -o -name "*.md" \
    -o -name "*.R" -o -name "*.pl" -o -name "*.yaml" -o -name "*.yml" \
    -o -name "*.json" -o -name "*.csv" -o -name "*.pt" \) -exec dos2unix {} + 2>/dev/null || true

# Set executable permissions on scripts
find . -type f \( -name "*.sh" -o -name "*.py" -o -name "*.pl" -o -name "*.R" \) \
    -exec chmod +x {} + 2>/dev/null || true

echo "Setup complete."
