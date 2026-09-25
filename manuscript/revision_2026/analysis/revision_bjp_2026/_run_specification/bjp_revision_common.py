#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import itertools
import json
import math
import os
import re
import sys
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

PRIMARY_CLASSES = ("Gs", "Gi/o", "Gq/11")


def ensure_project_imports(root: Path):
    src = root / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))
    from gpcr_icl2.microswitches.config import ProjectConfig
    from gpcr_icl2.microswitches.prospective_validation import _analysis_settings
    from gpcr_icl2.microswitches.publication_controls import (
        ExperimentSpec,
        FeatureSpec,
        _dynamic_splits,
        _metric_record,
        cross_validate_experiment,
    )
    return ProjectConfig, _analysis_settings, ExperimentSpec, FeatureSpec, _dynamic_splits, _metric_record, cross_validate_experiment


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(obj: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, sort_keys=True, default=str)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def unique_join(values: Iterable[Any], sep: str = ";") -> str:
    out = []
    for v in values:
        if pd.isna(v):
            continue
        s = str(v).strip()
        if not s or s.lower() == "nan":
            continue
        out.append(s)
    return sep.join(sorted(set(out)))


def norm_bool_series(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y"})


def discover_inputs(root: Path) -> dict[str, Path]:
    preferred = {
        "bound_frame": root / "results/microswitches/publication_controls_v9.2_final/publication_control_analysis_frame.csv",
        "active_none": root / "results/microswitches/publication_controls_v9.2_final/active_transducer_none_structures.csv",
        "strict_frame": root / "data/manual/strict_polymer_only_structures.csv",
        "entity_audit": root / "data/manual/receptor_only_entity_audit_final.csv",
        "active_transducers": root / "results/microswitches/coupling_specificity/data/active_transducer_structures_all.csv",
        "same_receptor": root / "results/microswitches/coupling_specificity/statistics/same_receptor_paired_comparisons.csv",
        "functional_provenance": root / "results/microswitches/sequence_interface_audit_v1/label_provenance/functional_label_provenance.tsv",
        "sequence_predictions": root / "results/microswitches/sequence_interface_audit_v1/predictions/out_of_fold_predictions.tsv",
        "categories": root / "results/microswitches/classification/microswitch_categories.csv",
    }
    imprint_candidates = list(root.glob("results/microswitches/manuscript_claims_nhr/**/transducer_imprint/transducer_imprint_matched_units.csv"))
    if not imprint_candidates:
        imprint_candidates = list(root.glob("results/**/transducer_imprint_matched_units.csv"))
    if imprint_candidates:
        matched = max(imprint_candidates, key=lambda p: p.stat().st_mtime)
        preferred["p0_matched"] = matched
        preferred["p0_imprint_dir"] = matched.parent
        preferred["p0_predictions"] = matched.parent / "transducer_imprint_model_predictions.csv"
        preferred["p0_audit"] = matched.parent / "transducer_imprint_cross_validation_audit.csv"
        preferred["p0_performance"] = matched.parent / "transducer_imprint_model_performance.csv"
    return preferred


def require_paths(paths: dict[str, Path], keys: Sequence[str]) -> None:
    missing = [f"{k}: {paths.get(k)}" for k in keys if k not in paths or not paths[k].exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs:\n" + "\n".join(missing))


def p0_base_features(paths: dict[str, Path]) -> list[str]:
    require_paths(paths, ["p0_matched"])
    cols = pd.read_csv(paths["p0_matched"], nrows=1).columns
    feats = [c[len("bound__"):] for c in cols if c.startswith("bound__")]
    if len(feats) != 92:
        raise RuntimeError(f"Expected frozen P0 feature set of 92 features; found {len(feats)}")
    return feats


def make_reference(free: pd.DataFrame, features: Sequence[str]) -> pd.DataFrame:
    work = free.copy()
    for f in features:
        if f not in work:
            work[f] = np.nan
        work[f] = pd.to_numeric(work[f], errors="coerce")
    agg = work.groupby("receptor_name", as_index=False)[list(features)].median(numeric_only=True)
    counts = work.groupby("receptor_name", as_index=False).agg(
        receptor_only_structure_count=("pdb_id", "nunique"),
        receptor_only_pdb_ids=("pdb_id", lambda x: ";".join(sorted(set(x.astype(str))))),
        receptor_only_chains=("chain_id", lambda x: unique_join(x)) if "chain_id" in work else ("receptor_name", lambda x: ""),
    )
    return agg.merge(counts, on="receptor_name", how="left")


def make_paired(bound: pd.DataFrame, free: pd.DataFrame, features: Sequence[str]) -> pd.DataFrame:
    b = bound.loc[bound["transducer_family"].isin(PRIMARY_CLASSES)].copy()
    common = sorted(set(b["receptor_name"].astype(str)) & set(free["receptor_name"].astype(str)))
    b = b.loc[b["receptor_name"].astype(str).isin(common)].copy()
    free = free.loc[free["receptor_name"].astype(str).isin(common)].copy()
    ref = make_reference(free, features)
    ids = [c for c in ("pdb_id", "chain_id", "receptor_name", "transducer_family") if c in b]
    for f in features:
        if f not in b:
            b[f] = np.nan
    raw = b[ids + list(features)].merge(ref, on="receptor_name", how="inner", suffixes=("__bound_raw", "__free_raw"))
    out = raw[ids + ["receptor_only_structure_count", "receptor_only_pdb_ids", "receptor_only_chains"]].copy()
    feature_cols = {}
    for f in features:
        bv = pd.to_numeric(raw[f"{f}__bound_raw"], errors="coerce")
        fv = pd.to_numeric(raw[f"{f}__free_raw"], errors="coerce")
        feature_cols[f"bound__{f}"] = bv
        feature_cols[f"receptor_only__{f}"] = fv
        feature_cols[f"delta__{f}"] = bv - fv
    return pd.concat([out, pd.DataFrame(feature_cols, index=out.index)], axis=1)


def classification_margin(row: pd.Series) -> float:
    true = str(row.get("transducer_family", ""))
    if true not in PRIMARY_CLASSES:
        return np.nan
    probs = {c: float(row[f"probability_{c}"]) for c in PRIMARY_CLASSES}
    if not all(np.isfinite(list(probs.values()))):
        return np.nan
    return probs[true] - max(v for c, v in probs.items() if c != true)


def exact_or_mc_signflip(values: Sequence[float], seed: int = 20291729, mc: int = 100000) -> float:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return np.nan
    obs = abs(float(np.mean(vals)))
    if len(vals) <= 20:
        null = np.fromiter(
            (abs(float(np.mean(vals * np.asarray(s, dtype=float)))) for s in itertools.product((-1.0, 1.0), repeat=len(vals))),
            dtype=float,
        )
        return float(np.mean(null >= obs - 1e-15))
    rng = np.random.default_rng(seed)
    hits = 0
    batch = 5000
    done = 0
    while done < mc:
        n = min(batch, mc - done)
        signs = rng.choice((-1.0, 1.0), size=(n, len(vals)))
        stats = np.abs((signs * vals[None, :]).mean(axis=1))
        hits += int(np.sum(stats >= obs - 1e-15))
        done += n
    return float((hits + 1) / (mc + 1))


def bootstrap_mean_ci(values: Sequence[float], iterations: int = 2000, seed: int = 20291729) -> tuple[float, float]:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    means = np.array([rng.choice(vals, len(vals), replace=True).mean() for _ in range(iterations)])
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def receptor_margin_table(predictions: pd.DataFrame, bound_name: str = "bound_geometry", free_name: str = "receptor_only_geometry") -> tuple[pd.DataFrame, dict[str, Any]]:
    p = predictions.loc[predictions["experiment"].isin([bound_name, free_name])].copy()
    p["classification_margin"] = p.apply(classification_margin, axis=1)
    keys = [c for c in ("pdb_id", "chain_id", "receptor_name", "transducer_family") if c in p]
    unit = p.pivot_table(index=keys, columns="experiment", values="classification_margin").reset_index()
    if bound_name not in unit or free_name not in unit:
        return pd.DataFrame(), {"status": "missing_prediction_branch"}
    unit["margin_gain"] = unit[bound_name] - unit[free_name]
    rec = unit.groupby("receptor_name", as_index=False).agg(
        margin_gain=("margin_gain", "mean"),
        bound_margin=(bound_name, "mean"),
        partner_free_margin=(free_name, "mean"),
        n_units=("transducer_family", "size"),
        captured_classes=("transducer_family", unique_join),
    )
    vals = rec["margin_gain"].to_numpy(float)
    lo, hi = bootstrap_mean_ci(vals)
    summary = {
        "status": "ok",
        "n_receptors": int(len(rec)),
        "n_positive": int(np.sum(vals > 0)),
        "n_negative": int(np.sum(vals < 0)),
        "n_zero": int(np.sum(vals == 0)),
        "mean_margin_gain": float(np.mean(vals)),
        "median_margin_gain": float(np.median(vals)),
        "mean_ci_low": lo,
        "mean_ci_high": hi,
        "signflip_p_two_sided": exact_or_mc_signflip(vals),
    }
    return rec, summary


def run_representation_cv(
    paired: pd.DataFrame,
    base_features: Sequence[str],
    root: Path,
    seed: int,
    quick: bool,
    nested: bool = True,
    solver_override: str | None = None,
    names: tuple[str, str] = ("bound_geometry", "receptor_only_geometry"),
) -> dict[str, Any]:
    ProjectConfig, _analysis_settings, ExperimentSpec, FeatureSpec, _, _, cross_validate_experiment = ensure_project_imports(root)
    config = ProjectConfig.from_dir(root / "config")
    settings = _analysis_settings(config, quick)
    if solver_override:
        settings["solver"] = solver_override
    experiments = [
        ExperimentSpec(names[0], FeatureSpec(names[0], tuple(f"bound__{f}" for f in base_features)), nested_tuning=nested),
        ExperimentSpec(names[1], FeatureSpec(names[1], tuple(f"receptor_only__{f}" for f in base_features)), nested_tuning=nested),
    ]
    preds=[]; perf=[]; audits=[]; coefs=[]
    for ex in experiments:
        res=cross_validate_experiment(paired, ex, settings, seed, minimum_fold_coverage=0.0, bootstrap_iterations=100 if quick else 2000)
        if not res["predictions"].empty: preds.append(res["predictions"])
        perf.append(dict(res["performance"]))
        if not res["audit"].empty: audits.append(res["audit"])
        if not res.get("coefficients", pd.DataFrame()).empty: coefs.append(res["coefficients"])
    predictions=pd.concat(preds,ignore_index=True) if preds else pd.DataFrame()
    return {
        "predictions": predictions,
        "performance": pd.DataFrame(perf),
        "audit": pd.concat(audits,ignore_index=True) if audits else pd.DataFrame(),
        "coefficients": pd.concat(coefs,ignore_index=True) if coefs else pd.DataFrame(),
    }


def run_pca_cv(paired: pd.DataFrame, base_features: Sequence[str], k: int, root: Path, seed: int, quick: bool) -> dict[str, Any]:
    _, _, _, _, _dynamic_splits, _metric_record, _ = ensure_project_imports(root)
    work=paired.reset_index(drop=True).copy()
    y=work["transducer_family"].astype(str).to_numpy()
    groups=work["receptor_name"].astype(str).to_numpy()
    splits=_dynamic_splits(y,groups,5,seed)
    classes=tuple(sorted(np.unique(y)))
    rows=[]; audit=[]
    for repr_name,prefix in [(f"bound_pca{k}","bound__"),(f"receptor_only_pca{k}","receptor_only__")]:
        probs=np.full((len(work),len(classes)),np.nan); fold_ids=np.full(len(work),-1)
        cols=[prefix+f for f in base_features]
        for fold,(tr,te) in enumerate(splits):
            if set(np.unique(y[tr]))!=set(classes):
                continue
            pipe=Pipeline([
                ("imputer",SimpleImputer(strategy="median",keep_empty_features=True)),
                ("scale",StandardScaler()),
                ("pca",PCA(n_components=min(k, len(tr)-1, len(cols)), random_state=seed+fold)),
                ("model",LogisticRegression(solver="lbfgs",penalty="l2",C=1.0,class_weight="balanced",max_iter=5000,random_state=seed+fold)),
            ])
            Xtr=work.loc[tr,cols].apply(pd.to_numeric,errors="coerce"); Xte=work.loc[te,cols].apply(pd.to_numeric,errors="coerce")
            pipe.fit(Xtr,y[tr]); local=pipe.predict_proba(Xte)
            for j,c in enumerate(pipe.named_steps["model"].classes_): probs[te,classes.index(str(c))]=local[:,j]
            fold_ids[te]=fold
            audit.append({"experiment":repr_name,"fold":fold,"train_n":len(tr),"test_n":len(te),"train_receptors":len(set(groups[tr])),"test_receptors":len(set(groups[te]))})
        pred=work[[c for c in ("pdb_id","chain_id","receptor_name","transducer_family") if c in work]].copy()
        pred["experiment"]=repr_name; pred["outer_fold"]=fold_ids
        for j,c in enumerate(classes): pred[f"probability_{c}"]=probs[:,j]
        rows.append(pred)
    preds=pd.concat(rows,ignore_index=True)
    perf=[]
    for exp in preds.experiment.unique():
        part=preds[preds.experiment==exp].copy(); mask=part[[f"probability_{c}" for c in classes]].notna().all(axis=1)
        y0=part.loc[mask,"transducer_family"].astype(str).to_numpy(); pr=part.loc[mask,[f"probability_{c}" for c in classes]].to_numpy(float)
        try: auc=roc_auc_score(y0,pr,labels=list(classes),multi_class="ovr",average="macro")
        except ValueError: auc=np.nan
        perf.append({"experiment":exp,"macro_roc_auc":auc,"n":int(mask.sum())})
    return {"predictions":preds,"performance":pd.DataFrame(perf),"audit":pd.DataFrame(audit)}


def paired_auc_from_predictions(pred: pd.DataFrame, bound_name: str, free_name: str) -> dict[str, float]:
    classes=sorted(PRIMARY_CLASSES)
    out={}
    for name in [bound_name,free_name]:
        p=pred[pred.experiment==name].copy(); p=p[p[[f"probability_{c}" for c in classes]].notna().all(axis=1)]
        try: auc=roc_auc_score(p.transducer_family.astype(str),p[[f"probability_{c}" for c in classes]].to_numpy(float),labels=classes,multi_class="ovr",average="macro")
        except ValueError: auc=np.nan
        out[name]=float(auc) if np.isfinite(auc) else np.nan
    out["delta_auc"] = out.get(bound_name,np.nan)-out.get(free_name,np.nan)
    return out


def input_manifest(root: Path, paths: dict[str, Path], used: Sequence[str], analysis: str, extra: dict[str, Any] | None = None) -> dict[str, Any]:
    files={}
    for k in used:
        p=paths.get(k)
        if p and p.exists():
            files[k]={"path":str(p.relative_to(root) if p.is_relative_to(root) else p),"sha256":sha256(p),"size":p.stat().st_size}
    prespec=root/"bjp_revision_2026_suite"/"PRESPEC.md"
    if not prespec.exists():
        alt=root/"PRESPEC.md"
        prespec=alt if alt.exists() else prespec
    obj={"analysis":analysis,"project_root":str(root),"inputs":files}
    if prespec.exists(): obj["prespec_sha256"]=sha256(prespec)
    if extra: obj.update(extra)
    return obj


def bh_fdr(pvalues: Sequence[float]) -> np.ndarray:
    p=np.asarray(pvalues,dtype=float); out=np.full(len(p),np.nan)
    mask=np.isfinite(p)
    vals=p[mask]
    if not len(vals): return out
    order=np.argsort(vals); ranked=vals[order]; n=len(ranked)
    q=ranked*n/np.arange(1,n+1); q=np.minimum.accumulate(q[::-1])[::-1]; q=np.clip(q,0,1)
    back=np.empty(n); back[order]=q; out[np.where(mask)[0]]=back
    return out


def sign_test_two_sided(pos: int, neg: int) -> float:
    n=pos+neg
    if n==0: return np.nan
    k=min(pos,neg)
    # exact binomial p for p=0.5
    prob=sum(math.comb(n,i) for i in range(k+1))/(2**n)
    return min(1.0,2*prob)


def spearman_safe(x: pd.Series, y: pd.Series) -> tuple[float,float,int]:
    z=pd.DataFrame({"x":pd.to_numeric(x,errors="coerce"),"y":pd.to_numeric(y,errors="coerce")}).dropna()
    if len(z)<3 or z.x.nunique()<2 or z.y.nunique()<2: return np.nan,np.nan,len(z)
    r,p=spearmanr(z.x,z.y)
    return float(r),float(p),len(z)
