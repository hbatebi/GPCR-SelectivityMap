from __future__ import annotations
from html import escape
from pathlib import Path
from typing import Any
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLOCKS = {
    "Contact network": "contact_",
    "Surface chemistry": "surface_",
    "Approximate electrostatics": "electro_",
    "Mechanical susceptibility": "mechanical_",
}

def selected_summary(row: dict[str, Any]) -> pd.DataFrame:
    requested = [
        "contact_network_density", "contact_salt_bridge_count", "contact_hbond_candidate_count",
        "surface_intracellular_total", "surface_positive_fraction", "surface_negative_fraction",
        "electro_cavity_potential_mean", "electro_cavity_potential_std",
        "mechanical_total_softness", "mechanical_tm6_opening_softness",
        "mechanical_orthosteric_interface_coupling_mean",
    ]
    records = []
    for feature in requested:
        if feature in row:
            value = row[feature]
            if isinstance(value, (int, float, np.integer, np.floating)) and np.isfinite(value):
                records.append({"feature": feature, "value": float(value)})
    return pd.DataFrame(records)

def profile_svg(row: dict[str, Any], output: str | Path) -> None:
    table = selected_summary(row)
    if table.empty:
        return
    p = Path(output); p.parent.mkdir(parents=True, exist_ok=True)
    values = table["value"].to_numpy(float)
    scaled = (values - np.nanmedian(values)) / (np.nanstd(values) + 1e-9)
    fig, ax = plt.subplots(figsize=(7.2, max(2.7, 0.28 * len(table))))
    y = np.arange(len(table))[::-1]
    ax.barh(y, scaled[::-1])
    ax.set_yticks(y); ax.set_yticklabels(table["feature"].tolist()[::-1], fontsize=8)
    ax.axvline(0, color="0.4", lw=0.8)
    ax.set_xlabel("Within profile standardized display value")
    ax.set_title("Selected GPCR receptor chain descriptors", loc="left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(p, format="svg", bbox_inches="tight")
    plt.close(fig)

def html_report(row: dict[str, Any], output: str | Path, figure_name: str | None = None) -> None:
    p = Path(output); p.parent.mkdir(parents=True, exist_ok=True)
    availability = []
    for key in ("contact", "surface", "electro", "mechanical"):
        availability.append((key, row.get(f"quality_{key}", "unknown"), row.get(f"error_{key}", "")))
    selected = selected_summary(row)
    rows = "".join(
        f"<tr><td>{escape(str(r.feature))}</td><td>{r.value:.6g}</td></tr>" for r in selected.itertuples()
    )
    avail_rows = "".join(
        f"<tr><td>{escape(name)}</td><td>{escape(str(status))}</td><td>{escape(str(error))}</td></tr>"
        for name, status, error in availability
    )
    image = f'<img src="{escape(figure_name)}" alt="Selected receptor descriptor profile">' if figure_name else ""
    html = f"""<!doctype html><html><head><meta charset='utf-8'><title>GPCR SelectivityMap profile</title>
    <style>body{{font-family:Arial,sans-serif;max-width:1000px;margin:32px auto;color:#111}}table{{border-collapse:collapse;width:100%}}th,td{{border:1px solid #bbb;padding:6px;text-align:left}}th{{background:#eee}}code{{background:#f4f4f4;padding:2px 4px}}img{{max-width:100%}}.warn{{background:#fff3cd;padding:12px}}</style></head><body>
    <h1>GPCR SelectivityMap structure profile</h1>
    <p><b>PDB identifier:</b> {escape(str(row.get('pdb_id','')))} &nbsp; <b>Chain:</b> {escape(str(row.get('chain_id','')))} &nbsp; <b>Receptor:</b> {escape(str(row.get('receptor_name','')))}</p>
    <div class='warn'><b>Scope:</b> These are static receptor chain descriptors. Mechanical values are elastic network susceptibility proxies, not free energies, trajectories, transition rates, or validated functional coupling probabilities.</div>
    <h2>Feature block audit</h2><table><tr><th>Block</th><th>Status</th><th>Error</th></tr>{avail_rows}</table>
    <h2>Selected descriptors</h2><table><tr><th>Feature</th><th>Value</th></tr>{rows}</table>{image}
    <h2>Reproducibility</h2><p>The complete machine readable profile is stored beside this report.</p></body></html>"""
    p.write_text(html, encoding="utf-8")
