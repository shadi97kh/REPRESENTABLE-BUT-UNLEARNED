# certmp -- certified reachability for monotone message passing

Two forward passes give the EXACT worst case over an exponentially large uncertainty set,
provided the aggregation is sum or max and the network is sign-constrained.

Application: RNA secondary structure is a Boltzmann ensemble, not a single structure.
Current siRNA/ASO pipelines take the minimum-free-energy structure and treat it as certain.
This repo instead certifies over the whole ensemble: "no structure in the ensemble drives
predicted off-target risk above tau", at the cost of two forward passes.

    certmp/models.py     monotone MPNN (numpy reference + torch trainable)
    certmp/reach.py      exact interval reachability + brute-force oracle
    certmp/ensemble.py   ViennaRNA base-pair probabilities -> edge lattice
    certmp/certify.py    threshold certificate + exact top-k rank stability
    certmp/data.py       dataset integrity verification against published statistics
    certmp/void.py       void detection
    certmp/provenance.py per-run provenance

    make setup && make test && make all
