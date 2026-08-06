from __future__ import annotations
from pathlib import Path
import json
import pandas as pd

def write_json(path: str | Path, payload: object) -> None:
    p = Path(path); p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")

def read_table(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    sep = "\t" if p.suffix.lower() in {".tsv", ".tab"} else ","
    return pd.read_csv(p, sep=sep, low_memory=False)

def write_table(frame: pd.DataFrame, path: str | Path) -> None:
    p = Path(path); p.parent.mkdir(parents=True, exist_ok=True)
    sep = "\t" if p.suffix.lower() in {".tsv", ".tab"} else ","
    frame.to_csv(p, sep=sep, index=False)
