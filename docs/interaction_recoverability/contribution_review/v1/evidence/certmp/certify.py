"""Certificates a designer would act on.

ONE-SIDED BY CONSTRUCTION. This module makes upper-bound claims only.

The lower bound was removed on 2026-09-04. experiments/f3b_relaxation_slack.py found that
the lattice minimum forced every mandatory pair present while real structures omit them,
so an actual Boltzmann-sampled structure scored BELOW the certified minimum on 2 of 20
sequences. A bound that real data violates is not a bound. It is gone from the API and
from the claims rather than being reported with a caveat.

reach.exact_interval still returns both lattice endpoints, because the extremality theorem
is a statement about both and f1/f5 test both. That is lattice arithmetic. Turning an
endpoint into a claim about the RNA ensemble is this module's job, and only the upper
endpoint earns one.
"""
from .reach import exact_interval


def certify_threshold(model, n, mandatory, optional, X, tau):
    """Upper-bound certificate: no structure in the lattice scores above tau.

    For the claim to cover the whole Boltzmann ensemble, build the lattice with
    ensemble.sound_lattice(seq, floor=0.0), which leaves mandatory empty and discards no
    pair. Then every secondary structure over the sequence is inside the lattice and
    monotonicity makes 'worst_case' a genuine upper bound on all of them; f7 confirms this
    empirically at 0 violations in 20000 sampled structures.

    Soundness is not free. At floor 0 the bound ran a median 355x the true sampled
    ensemble maximum in f7, so a certificate that passes is strong evidence while one that
    fails may only mean the relaxation is loose. Raising the floor tightens the bound and
    forfeits coverage; f7 has the curve.
    """
    _, hi, passes = exact_interval(model, n, mandatory, optional, X)
    return {"certified": hi <= tau, "worst_case": hi, "tau": tau,
            "forward_passes": passes, "k": len(optional),
            "lattice_size_log10": len(optional) * 0.30103,
            "margin": tau - hi,
            "sound_over_ensemble": len(mandatory) == 0}


def certify_topk(*args, **kwargs):
    """REMOVED 2026-09-04. Exact top-k rank certification needs a sound LOWER bound on
    each candidate, to show the k-th best worst case clears the (k+1)-th best best case.
    certmp no longer claims a lower bound over the ensemble, so this certificate has no
    foundation and is withdrawn rather than left standing on an invalid bound.
    experiments/f4_topk_stability.py is retired for the same reason; it certified nothing
    even when it ran.
    """
    raise NotImplementedError(certify_topk.__doc__)
