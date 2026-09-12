"""Revision cycle (protocol v2) tooling.

Everything here follows the post-run action plan. Three rules hold throughout:

  * historical records are READ-ONLY. Nothing in this package rewrites a file
    under `runs/` that a previous run produced. Corrections are written as NEW
    records that cite the record they correct.
  * every claim carries an explicit evidence kind (see `EVIDENCE_KINDS`), so an
    exact integer identity, a real-arithmetic proof, a float diagnostic, an
    enclosed numerical result, a learned synthetic result and a measured RNA
    result are never averaged into one sentence.
  * a number that came from the run report is a CLAIM until it is traced to the
    immutable record that produced it. `audit_report` does the tracing.
"""

EXACT_INTEGER = "exact_integer"          # integer/rational arithmetic, no rounding
EXACT_RATIONAL = "exact_rational"        # Fraction arithmetic, no rounding
REAL_PROOF = "real_arithmetic_proof"     # proved in exact real arithmetic, by hand
FLOAT_DIAGNOSTIC = "float_diagnostic"    # float64, tolerance-based, NOT certified
FLOAT_ENCLOSED = "float_enclosed"        # rigorous enclosure via directed rounding
LEARNED_SYNTHETIC = "learned_synthetic"  # fitted parameters, synthetic task
MEASURED_RNA = "measured_rna"            # fitted parameters, measured RNA labels

EVIDENCE_KINDS = (EXACT_INTEGER, EXACT_RATIONAL, REAL_PROOF, FLOAT_DIAGNOSTIC,
                  FLOAT_ENCLOSED, LEARNED_SYNTHETIC, MEASURED_RNA)

__all__ = ["EVIDENCE_KINDS", "EXACT_INTEGER", "EXACT_RATIONAL", "REAL_PROOF",
           "FLOAT_DIAGNOSTIC", "FLOAT_ENCLOSED", "LEARNED_SYNTHETIC", "MEASURED_RNA"]
