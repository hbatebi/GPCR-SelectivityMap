#!/usr/bin/env bash
set -euo pipefail
ROOT="$(pwd)"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-root) ROOT="$(cd "$2" && pwd)"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done
OUT="$ROOT/revision_bjp_2026"
[[ -d "$OUT" ]] || { echo "Missing $OUT" >&2; exit 1; }
(
  cd "$ROOT"
  find revision_bjp_2026 -type f -print0 | sort -z | xargs -0 sha256sum > revision_bjp_2026/SHA256SUMS.txt
  tar -czf BJP_revision_results_$(date +%Y%m%d_%H%M%S).tar.gz revision_bjp_2026
  ls -lh BJP_revision_results_*.tar.gz | tail -1
)
