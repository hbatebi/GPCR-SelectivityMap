from pathlib import Path
import numpy as np
import pandas as pd
from gpcr_selectivitymap.benchmark import run_benchmark

def test_grouped_benchmark_smoke(tmp_path: Path) -> None:
    rng=np.random.default_rng(3); rows=[]; labels=["Gs","Gi/o","Gq/11"]
    for i in range(18):
        label=labels[i%3]
        for rep in range(2):
            rows.append({"receptor_name":f"r{i}","transducer_family":label,
                         "contact_signal":float(i%3)+rng.normal(0,.2),
                         "surface_signal":float(i%3)+rng.normal(0,.2),
                         "mechanical_signal":rng.normal()})
    table=tmp_path/"features.tsv"; pd.DataFrame(rows).to_csv(table,sep="\t",index=False)
    cfg=tmp_path/"models.yaml"; cfg.write_text("models:\n  combined:\n    prefixes: [contact_, surface_, mechanical_]\n")
    result=run_benchmark(features=table,target="transducer_family",group="receptor_name",config=cfg,output=tmp_path/"out",folds=3,bootstraps=20,permutations=2,seed=7)
    assert len(result)==1
    assert 0 <= result.iloc[0].macro_roc_auc <= 1
    assert (tmp_path/"out/model_metrics.tsv").exists()
