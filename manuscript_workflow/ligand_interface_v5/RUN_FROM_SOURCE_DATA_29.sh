#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
ARCHIVE="$ROOT/source_data/natcomms_v5/Source_Data_29.zip"
OUT="${GPCR_LIGAND_WORKDIR:-$HERE/reproduction_workdir}"
[ -r "$ARCHIVE" ] || { echo "Missing $ARCHIVE" >&2; exit 2; }
mkdir -p "$OUT"
python3 - "$ARCHIVE" "$OUT" <<'PYINNER'
import zipfile,sys,os
with zipfile.ZipFile(sys.argv[1]) as z:
    target=os.path.realpath(sys.argv[2]);
    for info in z.infolist():
        path=os.path.realpath(os.path.join(target,info.filename))
        if not (path==target or path.startswith(target+os.sep)):
            raise SystemExit('Unsafe path in Source Data 29: '+info.filename)
    z.extractall(target)
PYINNER
(cd "$OUT" && bash RUN_ALL.sh)
echo "Complete results at $OUT/GPCR_Ligand_Function_v5_COMPLETE_RESULTS.zip"
