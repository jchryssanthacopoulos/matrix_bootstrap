#!/bin/bash
# usage: scripts/watchdog.sh <pid> <kill_gb> [poll_s]   -- kills <pid> if its RSS exceeds <kill_gb>; logs the peak
pid=$1; lim=$2; poll=${3:-3}; peak=0
while kill -0 $pid 2>/dev/null; do
  rss=$(ps -o rss= -p $pid 2>/dev/null | tr -d ' '); [ -z "$rss" ] && break
  gb=$(echo "scale=2; $rss/1000000" | bc); (( $(echo "$gb > $peak" | bc) )) && peak=$gb
  if (( $(echo "$gb > $lim" | bc) )); then echo "WATCHDOG: killing $pid at $gb GB"; kill -9 $pid; break; fi
  sleep $poll
done
echo "watchdog: pid $pid finished, peak RSS $peak GB"
