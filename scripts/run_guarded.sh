#!/bin/bash
# usage: scripts/run_guarded.sh <kill_gb> <logfile> <command...>
# Runs the command in the background, watches ITS pid (not a wrapper shell) every 0.5 s, kills it above <kill_gb>
# GB RSS, and appends the peak RSS and exit code to the log.
lim=$1; log=$2; shift 2
"$@" > "$log" 2>&1 &
pid=$!
peak=0
while kill -0 $pid 2>/dev/null; do
  rss=$(ps -o rss= -p $pid 2>/dev/null | tr -d ' '); [ -z "$rss" ] && break
  gb=$(echo "scale=2; $rss/1000000" | bc); (( $(echo "$gb > $peak" | bc) )) && peak=$gb
  if (( $(echo "$gb > $lim" | bc) )); then echo "WATCHDOG: killing $pid at $gb GB" >> "$log"; kill -9 $pid; break; fi
  sleep 0.5
done
wait $pid 2>/dev/null; code=$?
echo "guarded: pid $pid exit code $code, peak RSS $peak GB" >> "$log"
