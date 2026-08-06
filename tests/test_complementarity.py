import numpy as np
import pandas as pd

from gpcr_selectivitymap.complementarity import (
    FAMILIES,
    build_family_templates,
    normalize_generic_number,
    pair_components,
    score_receptor_against_template,
    score_receptor_rows,
)


def test_generic_number_normalisation_preserves_trailing_zero() -> None:
    assert normalize_generic_number(1.50) == "1x50"
    assert normalize_generic_number(34.50) == "34x50"
    assert normalize_generic_number("5x68") == "5x68"


def test_pair_components_are_finite() -> None:
    values = pair_components("L", "I")
    assert set(("hydropathy", "charge", "combined_physics")).issubset(values)
    assert all(np.isfinite(values[key]) for key in ("hydropathy", "charge", "volume", "combined_physics"))
    assert values["hydropathy"] > pair_components("L", "D")["hydropathy"]


def test_receptor_template_scoring_and_ranking() -> None:
    template_gs = pd.DataFrame([
        {
            "receptor_generic_position": "34x51",
            "galpha_reference_position": 376,
            "galpha_kd": 4.5,
            "galpha_charge": 0.0,
            "galpha_volume": 166.7,
            "galpha_aromatic": 0.0,
            "template_weight": 1.0,
        }
    ])
    template_gio = template_gs.copy()
    template_gio["galpha_kd"] = -3.5
    template_gio["galpha_charge"] = -1.0
    template_gq = template_gio.copy()
    templates = {"Gs": template_gs, "Gi/o": template_gio, "Gq/11": template_gq}
    frame = pd.DataFrame([
        {
            "pdb_id": "TEST",
            "receptor_name": "receptor_x",
            "transducer_family": "Gs",
            "gpcrdb_34x51_expected_aa": "I",
            "residue_34x51_identity": "I",
        }
    ])
    prediction = score_receptor_rows(frame, templates, positions=("34x51",))
    combined = prediction[prediction["model"].eq("combined_physics")]
    assert combined.loc[combined["candidate_family"].eq("Gs"), "compatibility_score"].iloc[0] > combined.loc[combined["candidate_family"].eq("Gi/o"), "compatibility_score"].iloc[0]
    assert combined["predicted_family"].iloc[0] == "Gs"


def test_family_template_builds_training_only_edges() -> None:
    contacts = pd.DataFrame([
        {"pdb_id": "A", "transducer_family": "Gs", "receptor_generic_position": "34x51", "galpha_reference_position": 376, "min_distance": 3.0, "atom_contact_count": 4, "contact_weight": 1.0},
        {"pdb_id": "B", "transducer_family": "Gs", "receptor_generic_position": "34x51", "galpha_reference_position": 376, "min_distance": 3.2, "atom_contact_count": 3, "contact_weight": 1.0},
        {"pdb_id": "C", "transducer_family": "Gi/o", "receptor_generic_position": "34x51", "galpha_reference_position": 376, "min_distance": 3.1, "atom_contact_count": 2, "contact_weight": 1.0},
    ])
    profiles = pd.DataFrame([
        {"pdb_id": "A", "transducer_family": "Gs", "galpha_reference_position": 376, "galpha_kd": 4.5, "galpha_charge": 0.0, "galpha_volume": 166.7, "galpha_aromatic": 0.0},
        {"pdb_id": "B", "transducer_family": "Gs", "galpha_reference_position": 376, "galpha_kd": 3.8, "galpha_charge": 0.0, "galpha_volume": 166.7, "galpha_aromatic": 0.0},
        {"pdb_id": "C", "transducer_family": "Gi/o", "galpha_reference_position": 376, "galpha_kd": -3.5, "galpha_charge": -1.0, "galpha_volume": 111.1, "galpha_aromatic": 0.0},
    ])
    receptor_profiles = {"A": {"34x51": "I"}, "B": {"34x51": "L"}, "C": {"34x51": "D"}}
    templates = build_family_templates(contacts, profiles, receptor_profiles, ["A", "B", "C"], minimum_edge_count=1)
    assert len(templates["Gs"]) == 1
    assert templates["Gs"]["n_training_complexes"].iloc[0] == 2
    score = score_receptor_against_template({"34x51": "I"}, templates["Gs"])
    assert np.isfinite(score["combined_physics"])
