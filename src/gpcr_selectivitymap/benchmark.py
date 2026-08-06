from __future__ import annotations
from pathlib import Path
from typing import Any
import json
import numpy as np
import pandas as pd
import yaml
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, log_loss, balanced_accuracy_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from .io import read_table, write_json, write_table

DEFAULT_BLOCKS = {
    "contact_network": ["contact_", "availability_contact"],
    "surface_electrostatics": ["surface_", "electro_", "availability_surface", "availability_electro"],
    "mechanical_susceptibility": ["mechanical_", "availability_mechanical"],
    "all_new_static": ["contact_", "surface_", "electro_", "mechanical_", "availability_"],
}

def _macro_auc(y: np.ndarray, prob: np.ndarray, classes: np.ndarray) -> float:
    return float(roc_auc_score(y, prob, labels=classes, multi_class="ovr", average="macro"))

def _model_columns(frame: pd.DataFrame, spec: dict[str, Any]) -> list[str]:
    columns=[]
    for c in spec.get("columns", []):
        if c in frame.columns: columns.append(c)
    for prefix in spec.get("prefixes", []):
        columns.extend([c for c in frame.columns if c.startswith(prefix)])
    return sorted(set(columns))

def _build_pipeline(frame: pd.DataFrame, columns: list[str]) -> Pipeline:
    numeric=[c for c in columns if pd.api.types.is_numeric_dtype(frame[c])]
    categorical=[c for c in columns if c not in numeric]
    transformers=[]
    if numeric:
        transformers.append(("numeric",Pipeline([("impute",SimpleImputer(strategy="median",add_indicator=True)),("scale",StandardScaler())]),numeric))
    if categorical:
        transformers.append(("categorical",Pipeline([("impute",SimpleImputer(strategy="most_frequent")),("onehot",OneHotEncoder(handle_unknown="ignore"))]),categorical))
    prep=ColumnTransformer(transformers, remainder="drop")
    model=LogisticRegression(max_iter=20000,solver="lbfgs",C=1.0)
    return Pipeline([("preprocess",prep),("model",model)])

def _predictions(frame: pd.DataFrame, columns: list[str], target: str, group: str, folds: int, seed: int) -> pd.DataFrame:
    X=frame[columns]; y=frame[target].astype(str).to_numpy(); groups=frame[group].astype(str).to_numpy()
    classes=np.array(sorted(np.unique(y)))
    splitter=StratifiedGroupKFold(n_splits=folds,shuffle=True,random_state=seed)
    records=[]
    for fold,(train,test) in enumerate(splitter.split(X,y,groups)):
        pipe=_build_pipeline(frame,columns)
        pipe.fit(X.iloc[train],y[train])
        probability=pipe.predict_proba(X.iloc[test])
        model_classes=pipe.named_steps["model"].classes_
        aligned=np.zeros((len(test),len(classes)))
        for j,label in enumerate(model_classes): aligned[:,np.where(classes==label)[0][0]]=probability[:,j]
        predicted=classes[np.argmax(aligned,axis=1)]
        for row_index, actual, pred, probs in zip(test,y[test],predicted,aligned):
            record={"row_index":int(row_index),"fold":fold,"receptor_name":str(groups[row_index]),"actual":actual,"predicted":pred}
            record.update({f"probability_{label}":float(probs[j]) for j,label in enumerate(classes)})
            records.append(record)
    return pd.DataFrame(records), classes

def _cluster_bootstrap(pred: pd.DataFrame, classes: np.ndarray, iterations: int, seed: int) -> tuple[float,float]:
    rng=np.random.default_rng(seed); groups=pred["receptor_name"].unique(); vals=[]
    pcols=[f"probability_{c}" for c in classes]
    for _ in range(iterations):
        chosen=rng.choice(groups,size=len(groups),replace=True)
        parts=[]
        for k,g in enumerate(chosen):
            part=pred[pred.receptor_name.eq(g)].copy(); part["_boot_group"]=k; parts.append(part)
        sample=pd.concat(parts,ignore_index=True)
        try: vals.append(_macro_auc(sample.actual.to_numpy(),sample[pcols].to_numpy(),classes))
        except ValueError: pass
    if not vals: return float("nan"),float("nan")
    return float(np.quantile(vals,.025)),float(np.quantile(vals,.975))

def _group_permutation(frame: pd.DataFrame, columns: list[str], target: str, group: str, folds: int, seed: int, observed: float, iterations: int) -> float:
    if iterations <= 0: return float("nan")
    rng=np.random.default_rng(seed+331); map_table=frame[[group,target]].drop_duplicates(group)
    null=[]
    for i in range(iterations):
        shuffled=map_table[target].to_numpy().copy(); rng.shuffle(shuffled)
        mapping=dict(zip(map_table[group].astype(str),shuffled))
        perm=frame.copy(); perm[target]=perm[group].astype(str).map(mapping)
        try:
            pred,classes=_predictions(perm,columns,target,group,folds,seed+i+1)
            pcols=[f"probability_{c}" for c in classes]
            null.append(_macro_auc(pred.actual.to_numpy(),pred[pcols].to_numpy(),classes))
        except Exception: continue
    return (1+sum(v>=observed for v in null))/(1+len(null)) if null else float("nan")

def run_benchmark(*, features: str | Path, target: str, group: str, output: str | Path,
                  config: str | Path | None = None, folds: int = 5, seed: int = 20272729,
                  bootstraps: int = 500, permutations: int = 100) -> pd.DataFrame:
    frame=read_table(features)
    for c in (target,group):
        if c not in frame: raise ValueError(f"Missing required column: {c}")
    frame=frame.dropna(subset=[target,group]).reset_index(drop=True)
    if config:
        cfg=yaml.safe_load(Path(config).read_text()) or {}
        models=cfg.get("models",cfg)
    else:
        models={name:{"prefixes":prefixes} for name,prefixes in DEFAULT_BLOCKS.items()}
    out=Path(output); out.mkdir(parents=True,exist_ok=True)
    metrics=[]; all_predictions=[]
    for name,spec in models.items():
        columns=_model_columns(frame,spec or {})
        if not columns: continue
        pred,classes=_predictions(frame,columns,target,group,folds,seed)
        pcols=[f"probability_{c}" for c in classes]
        auc=_macro_auc(pred.actual.to_numpy(),pred[pcols].to_numpy(),classes)
        lo,hi=_cluster_bootstrap(pred,classes,bootstraps,seed)
        ll=float(log_loss(pred.actual,pred[pcols],labels=classes))
        bal=float(balanced_accuracy_score(pred.actual,pred.predicted))
        pvalue=_group_permutation(frame,columns,target,group,folds,seed,auc,permutations)
        metrics.append({"model":name,"n_rows":len(frame),"n_groups":frame[group].nunique(),"n_features":len(columns),"macro_roc_auc":auc,"ci_low":lo,"ci_high":hi,"log_loss":ll,"balanced_accuracy":bal,"grouped_permutation_p":pvalue})
        pred.insert(0,"model",name); all_predictions.append(pred)
    metric_frame=pd.DataFrame(metrics).sort_values("macro_roc_auc",ascending=False)
    write_table(metric_frame,out/"model_metrics.tsv")
    if all_predictions: write_table(pd.concat(all_predictions,ignore_index=True),out/"out_of_fold_predictions.tsv")
    write_json(out/"analysis_manifest.json",{"target":target,"group":group,"folds":folds,"seed":seed,"bootstraps":bootstraps,"permutations":permutations,"models":models,"scope":"Static receptor chain association benchmark. Outputs are not validated functional coupling probabilities."})
    return metric_frame
