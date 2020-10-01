#!/bin/bash

# Add FastQC files
if [[ -d data/raw/FastQC ]]; then
  echo "Adding raw FastQC files"
  git add data/raw/FastQC
fi
if [[ -d data/trimmed/FastQC ]]; then
  echo "Adding trimmed FastQC files"
  git add data/trimmed/FastQC
fi

# Add trimming logs
if [[ -d data/trimmed/logs ]]; then
  echo "Adding trimming logs"
  git add -f data/trimmed/logs
fi

# Add STAR logs & counts
if [[ -d data/aligned ]]; then
  echo "Adding STAR logs and updating counts"
  git add data/aligned
fi
