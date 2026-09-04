# certmp -- certified reachability for monotone message passing

Two forward passes give the EXACT worst case over an exponentially large uncertainty set,
provided the aggregation is sum or max and the network is sign-constrained.

Application: RNA secondary structure is a Boltzmann ensemble, not a single structure.
Current siRNA/ASO pipelines take the minimum-free-energy structure and treat it as certain.
This repo instead certifies over the whole ensemble: "no structure in the ensemble drives
predicted off-target risk above tau", at the cost of two forward passes.

The upper bound is sound over the lattice and loose by a measured median factor of 1.82
against Boltzmann-sampled structures. The lower bound is NOT sound over the ensemble.
See RESULTS.md, section F3b, before quoting any number from here.

    certmp/models.py     monotone MPNN (numpy reference + torch trainable)
    certmp/reach.py      exact interval reachability + brute-force oracle
    certmp/ensemble.py   ViennaRNA base-pair probabilities -> edge lattice
    certmp/certify.py    threshold certificate + exact top-k rank stability
    certmp/void.py       void detection
    certmp/provenance.py per-run provenance

    certmp/data.py       dataset integrity verification against published statistics

    make setup && make test && make all      # theorem + RNA experiments, no download
    make data && make f6                     # fetch Huesken et al. 2005, train, evaluate

`make data` pulls the Huesken table in its published 2182/249 split from the DSIR
redistribution and refuses it unless it matches the statistics in the paper. The data is
not vendored here.
