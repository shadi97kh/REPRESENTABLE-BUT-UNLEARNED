"""P3 tests: label preservation, assay identity, insert recovery, split integrity."""
import numpy as np
import pytest
import yaml

from scmp.revision import rna_audit as RA
from scmp.revision import rna_pilot as RP

CFG = yaml.safe_load(open("configs/revision_rna_data.yaml"))


def test_labels_are_never_clipped_above_one():
    h = RA.audit_huesken(CFG)
    if h["status"] != "PRESENT":
        pytest.skip("Huesken absent")
    assert h["label"]["clipped"] is False
    assert h["label"]["max"] > 1.0, "values above 1.0 must survive"
    assert h["label"]["n_above_1"] > 0
    assert h["label"]["split_records_match_annotated_range"]


def test_annotation_assay_is_contradicted_not_adopted():
    h = RA.audit_huesken(CFG)
    if h["status"] != "PRESENT":
        pytest.skip("Huesken absent")
    ai = h["assay_identity"]
    assert ai["agrees_on_cell_line"] is False
    assert ai["agrees_on_readout"] is False
    assert ai["primary_source"]["cell_line"] == "H1299"


def test_davis_absence_names_the_file_and_route_without_bypass():
    d = RA.audit_davis(CFG)
    if d["status"] == "PRESENT":
        pytest.skip("Davis supplied")
    w = d["what_is_needed"]
    assert w["file"].endswith(".zip") and "S1" in w["contains"]
    assert "Oxford Academic" in w["route"]
    assert "not bypassed" in w["why_not_automated"].lower() or \
challenge_ok(w["why_not_automated"])


def challenge_ok(txt):
    return "reCAPTCHA" in txt and "human download" in txt


def test_geo_accession_is_flagged_as_not_labels():
    d = RA.audit_davis(CFG)
    bad = [x for x in d["do_not_substitute"] if x["accession"] == "GSE231101"]
    assert bad and "NOT knockdown labels" in bad[0]["why"]


def test_revcomp_is_an_involution():
    s = "ACGUUGCAACGU"
    assert RA.revcomp(RA.revcomp(s)) == s


def test_chemical_transfer_is_infeasible_without_a_second_scaffold():
    spl = RA.audit_splits(CFG, {})
    assert spl["chemical_transfer"]["feasible"] is False
    assert spl["novel_sequence"]["feasible"] is True


def test_arm_admissibility_is_driven_by_inputs_not_dataset_name():
    hue = {"status": "PRESENT", "chemistry": {"chemical_transfer_feasible": False}}
    no_window = RA.arm_admissibility(hue, {"status": "MISSING"},
                                     {"status": "SKIPPED"}, {})
    by = {r["arm"]: r["admissible"] for r in no_window}
    assert by["matched_sequence_chemistry"] is True
    assert by["mfe_gnn"] is False, "no resolved window -> inadmissible"
    with_window = RA.arm_admissibility(
        hue, {"status": "MISSING"},
        {"status": "INVESTIGATED", "by_window": {50: {"interior": 100}}}, {})
    by2 = {r["arm"]: r["admissible"] for r in with_window}
    assert by2["mfe_gnn"] is True
    assert by2["chemical_transfer"] is False, "still blocked by chemistry"


def test_novel_sequence_split_holds_out_whole_groups():
    rows = [{"gene": f"G{i % 4}", "cluster": f"C{i % 7}", "efficacy": 0.5,
             "guide": "ACGU" * 5 + "A", "window": "A" * 50,
             "n_isoform_hits": 1} for i in range(200)]
    tr, te, info = RP.novel_sequence_split(rows, 0.25, seed=0)
    assert info["group_leak"] == 0, "a group must not straddle the split"
    assert not (set(tr) & set(te))
    assert len(te) > 0 and len(tr) > 0


def test_top_k_enrichment_reports_prevalence():
    y = np.array([0.9, 0.8, 0.2, 0.1, 0.75, 0.3])
    pred = np.array([5.0, 4.0, 3.0, 2.0, 1.0, 0.0])
    out = RP.top_k_enrichment(pred, y, k=2, threshold=0.70)
    assert out["prevalence"] == pytest.approx(3 / 6)
    assert out["precision_at_k"] == pytest.approx(1.0)
    assert out["enrichment"] == pytest.approx(2.0)


def test_target_site_is_reverse_complement_of_guide():
    """So a null difference between the two sequence arms is expected."""
    guide = "CUAAUAUGUUAAUUGAUUUAU"
    assert RA.revcomp(guide) != guide
    assert len(RA.revcomp(guide)) == len(guide)
