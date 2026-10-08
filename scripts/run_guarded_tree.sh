#!/bin/bash
# usage: scripts/run_guarded_tree.sh <kill_gb> <logfile> <command...>
# Like run_guarded.sh, but watches the TOTAL RSS of the command and all its descendants (forked worker pools
# included) every 0.5 s, kills the whole tree above <kill_gb> GB, and appends the peak total RSS and exit code.
lim=$1; log=$2; shift 2
"$@" > "$log" 2>&1 &
pid=$!
peak=0
tree_rss() {   # total RSS (KB) of $1 and all descendants
  ps -A -o pid=,ppid=,rss= | awk -v root="$1" '
    { par[$1]=$2; rss[$1]=$3 }
    END { tot=0; for (p in par) { q=p; while (q!="" && q!=root && q!=1 && q!=0) q=par[q]; if (q==root) tot+=rss[p] } print tot }'
}
tree_pids() {
  ps -A -o pid=,ppid= | awk -v root="$1" '
    { par[$1]=$2 }
    END { for (p in par) { q=p; while (q!="" && q!=root && q!=1 && q!=0) q=par[q]; if (q==root) print p } }'
}
while kill -0 $pid 2>/dev/null; do
  rss=$(tree_rss $pid); [ -z "$rss" ] && break
  gb=$(echo "scale=2; $rss/1000000" | bc); (( $(echo "$gb > $peak" | bc) )) && peak=$gb
  if (( $(echo "$gb > $lim" | bc) )); then
    echo "WATCHDOG: killing tree of $pid at $gb GB" >> "$log"
    for q in $(tree_pids $pid); do kill -9 $q 2>/dev/null; done
    break
  fi
  sleep 0.5
done
wait $pid 2>/dev/null; code=$?
echo "guarded-tree: pid $pid exit code $code, peak total RSS $peak GB" >> "$log"
