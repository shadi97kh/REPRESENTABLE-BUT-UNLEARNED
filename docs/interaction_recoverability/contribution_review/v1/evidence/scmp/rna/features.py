"""Fixed covariates X, and target-window construction with explicit isoforms.

SEPARATION OF ROLES. Two different objects, and conflating them is the error this module
exists to prevent:

  X  -- FIXED covariates: guide sequence, guide chemistry, assay identity. These never
        depend on the latent structure.
  z  -- UNCERTAIN adjacency: optional base pairs of the TARGET-RNA local window.

The energy model is applied to the target window only. The target window is native,
unmodified RNA, so a nearest-neighbour model is in scope for it. The guide's chemistry
(2'-OMe, 2'-F, phosphorothioate, ...) enters **only** through X and is never handed to a
folding routine, because no validated nearest-neighbour parameters exist for those
chemistries. `assert_chemistry_not_folded` enforces that at run time.

ISOFORMS ARE EXPLICIT. A guide that matches several transcript isoforms produces one
window per isoform, each carrying its own accession, strand and coordinates. Exon unions
are refused: a union transcript is a sequence no molecule has, and folding it produces an
ensemble that does not exist.

CONTEXT SUBSTITUTION IS REFUSED BY DEFAULT. If the assay used a reporter construct, the
native transcript is not the assayed context. Producing native windows for such a dataset
requires an explicit, recorded override; there is no silent path.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

NUC = "ACGU"
UNMODIFIED = "unmodified"


class ContextSubstitutionRefused(Exception):
    pass


class ExonUnionRefused(Exception):
    pass


class ChemistryFoldingRefused(Exception):
    pass


def assert_chemistry_not_folded(chemistry: str, what: str):
    """Guard: the folding routine may only ever see native, unmodified target RNA."""
    if chemistry and chemistry != UNMODIFIED:
        raise ChemistryFoldingRefused(
            f"refusing to fold {what} with chemistry {chemistry!r}: no validated "
            f"nearest-neighbour parameters exist for modified chemistry. Guide "
            f"chemistry belongs in X, not in the energy model.")


# --------------------------------------------------------------------- covariates
def guide_features(seq: str) -> np.ndarray:
    """Non-negative, position-resolved guide covariates (H4)."""
    seq = seq.upper().replace("T", "U")
    L = len(seq)
    X = np.zeros(4 * L)
    for t, c in enumerate(seq):
        if c in NUC:
            X[4 * t + NUC.index(c)] = 1.0
    return X


def chemistry_features(chemistry: str, vocabulary) -> np.ndarray:
    """One-hot over a DECLARED chemistry vocabulary. An unseen chemistry is an error,
    not a silent zero vector, because a silent zero would make a modified guide look
    unmodified."""
    vocab = list(vocabulary)
    if chemistry not in vocab:
        raise KeyError(f"chemistry {chemistry!r} not in declared vocabulary {vocab}")
    v = np.zeros(len(vocab))
    v[vocab.index(chemistry)] = 1.0
    return v


def assay_features(assay: str, vocabulary) -> np.ndarray:
    vocab = list(vocabulary)
    if assay not in vocab:
        raise KeyError(f"assay {assay!r} not in declared vocabulary {vocab}")
    v = np.zeros(len(vocab))
    v[vocab.index(assay)] = 1.0
    return v


def build_covariates(guide: str, chemistry: str, assay: str,
                     chem_vocab, assay_vocab) -> np.ndarray:
    """X = [guide one-hot | chemistry one-hot | assay one-hot], all non-negative."""
    return np.concatenate([guide_features(guide),
                           chemistry_features(chemistry, chem_vocab),
                           assay_features(assay, assay_vocab)])


def per_node_covariates(X_flat: np.ndarray, n_nodes: int) -> np.ndarray:
    """Broadcast the fixed covariate vector onto every window node.

    The covariates describe the guide and the assay, not any individual nucleotide of
    the target, so they are constant across nodes. Positional identity of the target
    window is carried by the adjacency, not by X.
    """
    return np.tile(X_flat[None, :], (n_nodes, 1))


# ------------------------------------------------------------------------ windows
@dataclass
class TargetWindow:
    sequence: str
    length: int
    accession: str | None
    isoform: str | None
    species: str | None
    strand: str | None
    start: int | None
    context_source: str          # native_transcript | reporter_construct | declared
    provenance: dict = field(default_factory=dict)

    def __post_init__(self):
        if len(self.sequence) != self.length:
            raise ValueError("window length does not match the sequence")


def make_declared_window(sequence: str, length: int, name="declared") -> TargetWindow:
    """A window whose sequence is given outright -- used for tiny audits and tests,
    where no transcript resolution is involved."""
    return TargetWindow(sequence[:length], length, None, None, None, None, None,
                        "declared", {"name": name})


def windows_for_isoforms(guide_site, isoforms, length: int, *,
                         context_source: str,
                         allow_context_substitution: bool = False,
                         substitution_justification: str = ""):
    """One window per isoform. Never a union.

    `isoforms` is a list of dicts with accession, isoform, species, strand and the
    transcript sequence. Each is resolved separately and carries its own coordinates.
    """
    if context_source == "reporter_construct" and not allow_context_substitution:
        raise ContextSubstitutionRefused(
            "the assayed context is a reporter construct; native transcript windows are "
            "not the assayed context. Set allow_context_substitution with a written "
            "justification to proceed, and the substitution will be recorded in every "
            "output.")
    out = []
    for iso in isoforms:
        seq = iso["sequence"]
        pos = seq.find(guide_site)
        if pos < 0:
            continue
        centre = pos + len(guide_site) // 2
        lo = max(0, centre - length // 2)
        win = seq[lo:lo + length]
        if len(win) < length or set(win) - set(NUC):
            continue
        out.append(TargetWindow(
            win, length, iso.get("accession"), iso.get("isoform"),
            iso.get("species"), iso.get("strand"), lo, context_source,
            {"substituted": bool(allow_context_substitution),
             "justification": substitution_justification,
             "site_offset_in_window": pos - lo}))
    return out


def refuse_exon_union(*_args, **_kwargs):
    raise ExonUnionRefused(
        "folding an exon union is refused: the union is a sequence no transcript has, "
        "so its ensemble does not correspond to any molecule. Model isoforms explicitly.")
