"""F4: RETIRED 2026-09-04.

Exact top-k rank certification compares the k-th best candidate's worst case against the
(k+1)-th best candidate's BEST case. That needs a sound lower bound on every candidate.
f3b showed the lattice minimum is not one: it forces mandatory pairs present while real
structures omit them, and actual sampled structures scored below it. The lower bound was
removed from certmp.certify, so this certificate has no foundation and is withdrawn.

It certified nothing when it did run: margins were negative at every k, from -170 at k=1
to -5024 at k=10, because the reachable intervals straddled every cut.

Reinstating this needs a genuinely sound lower bound over the ensemble, which certmp does
not currently have. This file is kept so the retirement is visible in the repository
rather than being a silent deletion.
"""
from certmp.provenance import new_run, record

REASON = ("withdrawn: exact top-k certification requires a sound ensemble lower bound; "
          "certmp.certify no longer claims one (see f3b, f7)")

def main():
    run = new_run("f4_topk_stability", dict(status="retired"))
    print("F4 is RETIRED.")
    print(f"  {REASON}")
    print("  prior result, for the record: no ranking certified at k = 1, 3, 5 or 10; "
          "margins -170.4 to -5024.2")
    record(run, "f4_topk_stability", {"status": "retired", "reason": REASON,
                                      "prior_result": {"certified_any_k": False,
                                                       "margins": {1: -170.4270, 3: -1786.7254,
                                                                   5: -2949.8848, 10: -5024.2232}}})

if __name__ == "__main__": main()
