#!/bin/bash
# For every job directory in tvd-fast whose name matches PATTERN (running or
# completed), reads the latest avg(Ke1)/avg(Ke2) from its most recent
# goldfish.o* file and prints the min/max/median/mean of the ratio.
# Usage: bash ke_ratio.sh [pattern ...] [-w|--watch] [-l|--list]
#   pattern: bash glob matched against directory names, e.g.
#   'skd15_Ro3p5_Rm1e3_Ny16_Lz2pi' or 'skd*' (default: '*'). Multiple patterns
#   OR together -- this also makes an unquoted glob (which your shell
#   pre-expands into many literal names) work correctly, since each expanded
#   name is just matched as its own pattern.
#   -l|--list: print every job's ratio (sorted) once, instead of the
#   min/max/median/mean summary. Ignores -w/--watch.

PATTERNS=()
WATCH=false
LIST=false
for arg in "$@"; do
    if [ "$arg" = "-w" ] || [ "$arg" = "--watch" ]; then
        WATCH=true
    elif [ "$arg" = "-l" ] || [ "$arg" = "--list" ]; then
        LIST=true
    else
        PATTERNS+=("$arg")
    fi
done
[ ${#PATTERNS[@]} -eq 0 ] && PATTERNS=("*")
LIST_PY=False
$LIST && LIST_PY=True

matching_dirs () {
    for d in */; do
        d="${d%/}"
        for p in "${PATTERNS[@]}"; do
            if [[ "$d" == $p ]]; then
                echo "$d"
                break
            fi
        done
    done
}

dirs=$(matching_dirs)
echo "matched $(echo -n "$dirs" | grep -c .) job dir(s)"

run_once () {
    echo "$dirs" | while read -r name; do
        [ -z "$name" ] && continue
        f=$(ls -t "$name"/goldfish.o* 2>/dev/null | head -1)
        [ -z "$f" ] && continue
        line=$(grep "avg(Ke1)" "$f" | tail -1)
        ke1=$(echo "$line" | grep -oP 'avg\(Ke1\)=\K[0-9.eE+-]+')
        ke2=$(echo "$line" | grep -oP 'avg\(Ke2\)=\K[0-9.eE+-]+')
        [ -z "$ke1" ] || [ -z "$ke2" ] && continue
        python3 -c "print('$name', $ke1/$ke2)"
    done > /tmp/ke_ratios.txt

    python3 -c "
vals = []
for l in open('/tmp/ke_ratios.txt'):
    parts = l.split()
    if len(parts) != 2:
        continue
    try:
        vals.append((float(parts[1]), parts[0]))
    except ValueError:
        continue
if not vals:
    print('no jobs with avg(Ke1)/avg(Ke2) data yet')
else:
    vals.sort()
    n = len(vals)
    if $LIST_PY:
        for r, name in vals:
            print(f'{r:.4g}  {name}')
    else:
        median = vals[n//2][0] if n % 2 else (vals[n//2-1][0] + vals[n//2][0]) / 2
        mean = sum(r for r, _ in vals) / n
        print(f'min={vals[0][0]:.4g} ({vals[0][1]})  max={vals[-1][0]:.4g} ({vals[-1][1]})  median={median:.4g}  mean={mean:.4g}')
"
}

if $LIST; then
    run_once
elif $WATCH; then
    while true; do
        printf '%s  ' "$(date +%H:%M:%S)"
        run_once
        sleep 5
    done
else
    run_once
fi
