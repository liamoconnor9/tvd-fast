#!/usr/bin/env bash
#
# cleanup_slicepoints.sh
#
# Deletes:
#   1. Any directory named "slicepoints" older than 4 months (based on mtime)
#   2. Any directory named "cleanx", "cleany", or "cleanz" (regardless of age)
#
# Does NOT touch files named cleanx.mp4, cleany.mp4, cleanz.mp4 since -type d
# only matches directories, never files.
#
# Usage:
#   ./cleanup_slicepoints.sh /path/to/parent          # dry run (default, safe)
#   ./cleanup_slicepoints.sh /path/to/parent --delete # actually deletes

set -euo pipefail

TARGET_DIR="${1:-}"
MODE="${2:-}"

if [[ -z "$TARGET_DIR" ]]; then
  echo "Usage: $0 /path/to/parent [--delete]"
  exit 1
fi

if [[ ! -d "$TARGET_DIR" ]]; then
  echo "Error: '$TARGET_DIR' is not a directory."
  exit 1
fi

DRY_RUN=true
if [[ "$MODE" == "--delete" ]]; then
  DRY_RUN=false
fi

# 4 months ~ 120 days (adjust if you want exact calendar months)
DAYS_OLD=120

echo "=== Scanning $TARGET_DIR ==="
echo

echo "--- 'slicepoints' directories older than $DAYS_OLD days ---"
while IFS= read -r -d '' dir; do
  if $DRY_RUN; then
    echo "[DRY RUN] Would delete: $dir"
  else
    echo "Deleting: $dir"
    rm -rf -- "$dir"
  fi
done < <(find "$TARGET_DIR" -type d -name "slicepoints" -mtime +"$DAYS_OLD" -print0)

echo
echo "--- 'cleanx', 'cleany', 'cleanz' directories (any age) ---"
while IFS= read -r -d '' dir; do
  if $DRY_RUN; then
    echo "[DRY RUN] Would delete: $dir"
  else
    echo "Deleting: $dir"
    rm -rf -- "$dir"
  fi
done < <(find "$TARGET_DIR" -type d \( -name "cleanx" -o -name "cleany" -o -name "cleanz" \) -print0)

echo
if $DRY_RUN; then
  echo "=== Dry run complete. No files were deleted. ==="
  echo "=== Re-run with --delete to actually remove these directories. ==="
else
  echo "=== Deletion complete. ==="
fi