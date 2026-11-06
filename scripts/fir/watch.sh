#!/bin/bash

# Exits on the first new job problem, or when no job is left, so that a background run wakes the session.
# Rerun after handling a problem: tasks already reported are remembered in seen.txt, next to this log.
# Usage: bash watch.sh <start date, YYYY-MM-DD> <first job id> [folder for seen.txt]
START=${1:?Usage: watch.sh <start date> <first job id> [folder]}
FIRST=${2:?Usage: watch.sh <start date> <first job id> [folder]}
SEEN=${3:-$HOME}/watch_seen.txt; touch "$SEEN"
NAMES='preprocess|host_index|metaphlan|motus3|motus4|sylph|sourmash|kraken'

while true; do
	new=$(sacct -u "$USER" -S "$START" -X -n -P -o JobID,JobName,State,ExitCode,Elapsed,NodeList 2>/dev/null \
		| awk -F'|' -v n="^($NAMES)" -v first="$FIRST" '{split($1, a, "_"); if (a[1] + 0 >= first && $2 ~ n) print}' \
		| grep -E '\|(FAILED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|BOOT_FAIL|DEADLINE|PREEMPTED)' \
		| grep -vFf <(cut -d'|' -f1 "$SEEN"; echo @none@))
	stuck=$(squeue -u "$USER" -h -t PD -o '%A %j %r' 2>/dev/null \
		| grep -E 'DependencyNeverSatisfied|BadConstraints|InvalidAccount|InvalidQOS|JobHeldAdmin')
	if [[ -n $new || -n $stuck ]]; then
		[[ -n $new ]] && echo "$new" >> "$SEEN"
		echo "PROBLEM $(date '+%F %T')"
		[[ -n $new ]] && { echo "$new" | awk -F'|' '{split($1, a, "_"); print a[1], $2, $3, $4}' | sort | uniq -c; echo "$new" | head -8; }
		[[ -n $stuck ]] && echo "STUCK: $stuck"
		exit 0
	fi
	left=$(squeue -u "$USER" -h -o '%j' 2>/dev/null | grep -cE "^($NAMES)")
	if [[ $left == 0 ]]; then
		sleep 120
		[[ $(squeue -u "$USER" -h -o '%j' 2>/dev/null | grep -cE "^($NAMES)") == 0 ]] && { echo "ALL DONE $(date '+%F %T')"; exit 0; }
	fi
	sleep 300
done
