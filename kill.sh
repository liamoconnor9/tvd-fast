#!/bin/bash
set -uo pipefail

CFG="$1/options.cfg"
# CFG="${1:-options.cfg}"

if [[ ! -f "$CFG" ]]; then
    echo "Error: $CFG not found" >&2
    exit 1
fi

# shellcheck disable=SC1090
source "$CFG"
if [[ -z "${jobid:-}" ]]; then
    echo "Error: jobid not set in $CFG" >&2
    exit 1
fi

echo "Cancelling job $jobid..."
scancel "$jobid"