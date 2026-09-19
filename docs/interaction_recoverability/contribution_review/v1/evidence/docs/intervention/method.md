# Query-conditioned graph intervention responses — preparation specification

This implementation is a graph conditional neural process (CNP). The response-only candidate and its capacity-matched generic CNP control are the same Python class and have 17,409 parameters. A graph representation and auxiliary losses currently establish an application and an empirical question, rather than a distinct ML operator. The [novelty matrix](novelty_matrix.md) makes the source comparisons explicit; the closest operator is the original [CNP](https://proceedings.mlr.press/v80/garnelo18a.html).

For a fixed predictor, construction procedure and assay context, write the complete input as (X(G,c)=(G,d(G),c)). An action is valid only if its registered transformation reconstructs every input required by that predictor. The estimand is

\[
v_f(a)=f(T(G,a),d(T(G,a)),c)-f(G,d(G),c),\qquad v_f(0)=0.
\]

For compatible commuting actions with all four states registered,

\[
I_f(a,b)=v_f(a\cup b)-v_f(a)-v_f(b).
\]

This is a reference finite difference. It is neither automatically a Shapley interaction nor a measured biological effect. Noncommuting transformations would require an ordered estimand such as (v_f(b\circ a)-v_f(a)-v_f(b)); this preparation does not implement an ordered-action experiment. Overlapping binary switches and component axes are rejected. The low-level `interaction` helper expects states obtained through the compatible constructor; it cannot infer action commutativity from four arbitrary tensors.

`GraphInput` carries atom/node features, adjacency, derived descriptors, context, construction ID, molecular/graph identity and optional E/Z bond channels. `FrozenQueryOracle` rejects a construction mismatch, freezes a supplied PyTorch predictor, and hashes its declared model identity plus every input tensor, construction and context. Effects reject context changes. Its cache records real uncached query counts. A model hash is accounting metadata and never reaches the explainer. Callers must hash the actual predictor configuration/weights before using the wrapper; the wrapper does not prove that a supplied string matches those weights.

The supported transformations are deliberately small and explicit:

| Interface | Reconstructed inputs and validity contract | Limit |
|---|---|---|
| `BinaryGraphFamily` | Three distinct node-state or edge switches; all eight states; four derived graph statistics rebuilt on every call | Analytic graph controls, with no chemical or RNA interpretation |
| `molecular_input` | RDKit sanitization; eight atom features; bond-order adjacency; absolute atom R/S and bond E/Z features; molecular weight, logP, TPSA and rotatable bonds | Four-descriptor candidate schema; incompatible with AGILE Mordred and TransMA 3D inputs |
| `replace_mapped_atom` | Complete unique atom maps, expected reference element, fixed connectivity, sanitization, descriptors and stereo regenerated | Rejects changes at specified tetrahedral centers; broader edits require an explicitly specified stereochemical product |
| `ProductLibrary` | Component tuple maps to an explicitly registered product SMILES; all four products must exist | No reaction prediction or atom correspondence for fragment substitutions; parsing alone cannot establish mixture composition or assay comparability |
| `rebuild_rna` | Valid full edited A/C/G/U sequence passed to the registered feature builder | Caller must implement and version all sequence/structure reconstruction; no RNA folding engine or RNA predictor is supplied |
| `fingerprint_descriptor_features` | Chirality-aware Morgan radius-2, 2,048 bits plus the same four freshly rebuilt descriptors | Candidate tabular representation, without fitted scaling or published-checkpoint compatibility |

Canonical isomeric SMILES identifies the molecular state while atom maps preserve correspondence separately. No salt stripping, tautomer correction, protonation normalization or stereochemical imputation is performed. Atom features use absolute CIP labels, avoiding a dependence on SMILES traversal direction. Unspecified stereochemistry stays unspecified; neither the representation nor canonicalization resolves mixtures. The audited lipid data make this limitation material.

The encoder applies two layers of sum-aggregated message passing with tanh activations and 32 hidden units. Node embeddings are equivariant to relabeling; their sum is invariant. An edit embedding concatenates pooled reference and edited graphs and their difference, reference and edited descriptors and their difference, and fixed context. It does not require a guessed fragment atom alignment. Distinct edits can still collide under this finite representation; invariance does not prove expressiveness.

For observed actions and signed responses (Q=\{(a_t,v_t)\}), the network computes

\[
z_Q=\operatorname{mean}_t h(e(G,a_t,c),v_t),\quad
\hat v(a\mid Q)=D(e(G,a,c),z_Q,\log(1+|Q|))
-D(e(G,0,c),z_Q,\log(1+|Q|)).
\]

An empty transcript has a zero set representation and count. Both the paired observation set and final response are invariant to query order. Identity centering is exact in this deterministic implementation. No model ID, predictor weights, training labels, assay labels or test-source embedding is an input. Query responses are indispensable observations, not biological supervision. Source labels remain outside this transcript interface.

The implemented loss terms are normalized response MSE, MSE of compatible second differences, and optional MSE of differences between two models' predicted response vectors. Model contrast requires their reference outputs to differ by at most 0.01 and uses model-specific transcripts. Loss weights are 0.25 for pair and contrast terms. Response, response-plus-pair and response-plus-pair-plus-contrast are separate arms; a favorable combination would still require a separate novelty argument. The gradient checks take backward passes without an optimizer or parameter update.

Two independent controls specify the target functions directly. `edge_switch` changes graph edges; `node_state` changes declared node attributes. Each has three switches and an exhaustive eight-state cube. Fixed analytic coefficient tables cover additive, pairwise and third-order responses. Another pair of functions agrees at the reference and every proper subset, but differs at the triple edit. Since their observed transcripts are identical, any deterministic estimate (m) at that edit obeys

\[
\max(|m-y_1|,|m-y_2|)\ge |y_1-y_2|/2.
\]

This triangle-inequality consequence is an information limit, not a new theorem. These predictors are named `AnalyticPredictor`; they are not trained GNNs. `GraphDescriptorRegressor` supplies an actual message-passing predictor interface but remains untrained, with no capability or R² claim.

The prospective source generator samples seven independent standard-normal monomial coefficients. Each independent seed creates two matched variants sharing first-order coefficients; four higher-order coefficients differ. Split seeds and the exact `SeedSequence` construction are frozen in [configuration v2](../../configs/intervention_screen_v2.json). Only a lazy generator exists: the main source-function collection and teacher cache have not been materialized. The 80-row exported table comes from fixed correctness fixtures and is not held-out screening evidence.

This generator has a known conditional mean: with action design matrix (X_Q), independent standard-normal coefficients and noiseless observations, (E[\theta\mid Q]=X_Q^+y_Q). `PolynomialGaussianResponseSurrogate` implements that mean with the full seven-monomial basis. At identical observations, it minimizes expected squared error under this declared prior. A CNP cannot have a systematic expected-error advantage over this correctly specified reference on these controls. Graph-feature seed variation also leaves the coefficient law unchanged; it does not create architecture generalization. A substantive learned-method study would need a justified additional mechanism and problem setting before a new preregistration.

| Baseline or ablation | Executable interface/status | Comparison contract |
|---|---|---|
| Additive/pairwise ridge and pairwise LASSO | `SparseResponseSurrogate`; exact least-squares correctness checks | Same signed responses and fixed/random action schedule; fit time charged |
| RBF GP posterior mean | `KernelResponseSurrogate`; checked on a tiny transcript | Identity-centered kernel; development-only hyperparameters |
| Known-prior polynomial GP | `PolynomialGaussianResponseSurrogate`; exact conditional calculation checked | Mandatory privileged control for the declared analytic family |
| Direct boosted surrogate | `BoostedResponseSurrogate`; unfitted interface | Refit on queried responses; identity-centered output; not a claimed ProxySPEX reproduction |
| Exact finite differences/cache | `FrozenQueryOracle` plus full enumeration | Charge first-time evaluation/cache construction; report warm reuse separately |
| Generic CNP and candidate arms | Same `GraphConditionalProcess` class and capacity; loss interfaces complete | Same queries, source splits and optimization schedule; no favorable seed selection |
| Responses zeroed/shuffled | Transcript interfaces checked | Same action states; distinguish input dependence from learned held-out utility |
| PGExplainer/FastSHAP/selective refinement/GraphSHAP-IQ | Applicable author source audited; no executed adapters | Native subgraph or attribution evaluations only; their native outputs do not equal signed molecular effects |
| SPEX/ProxySPEX | Source and valid binary-game adapter contract reviewed; external integration unexecuted | Query a binary vector by mapping its selected indices to registered edits, recover full values, then derive reference differences; include all query, fitting and extraction costs |
| Learned query policy | Excluded | No optimization, generated policy data or policy benefit claimed |

The 2026 [Extension Sufficiency Test](https://arxiv.org/abs/2601.20815) belongs to the native subgraph track. For the response track, matched transcripts, response removal/shuffling and held-out function evaluation test different information-flow failures. None of these checks alone certifies explanation quality.

The primary endpoint is mean squared signed-response error, divided by a train-only per-task squared RMS scale, integrated exactly as a step function over log total runtime from 0.001 to 1 second. Total time includes each method's allocated upfront cost divided by 1,000 deployments. Zero effects apply before a complete output; later outputs replace the entire prediction vector, including estimates for unqueried actions. Deadlines retain the latest completed output. Sensitivity to other deployment counts is secondary.

Continuation requires at least 10% relative improvement on each task over its development-selected strongest applicable baseline, with the lower endpoint of a paired 95% improvement interval above zero. Thirty-two independent function seeds per task are the units: matched variants, graph seeds and initialization runs are averaged inside each unit before a 10,000-resample bootstrap. No choice of endpoint, interval, split, baseline exclusion or seed may depend on held-out results. The generic CNP and exact/known-prior controls remain mandatory. This configuration is reviewable preparation; it is not an authorized training run or a claim that these controls can establish a new method.

The final [18-check result](../../runs/intervention/tests-v4.json) covers valid and invalid transformations, full descriptor rebuilding, atom/bond stereochemistry, equivariance and set invariance, identity centering, transcript response dependence, indistinguishability, unchanged weights after backward passes, exact polynomial responses, cache input keys, zero queries, runtime deadlines and the known-prior control. Tests do not establish generalization, predictive quality, convergence or biological utility.
