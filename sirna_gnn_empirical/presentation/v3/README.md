# Current workflow and result figures

These are the existing campaign v3 figures used by manuscript v8, **Do Graphs Learn the Chemistry? Controlled Evidence from siRNA Prediction**. No fitting or scientific reevaluation was performed to publish them.

![Architecture and separate supervision protocols](figures/workflow.png)

The supplied 3300 × 1950 PNG is preserved byte for byte. The standalone PDF embeds that image without a Figure caption; it is a raster PDF, not a vector drawing. The superseded SVG is not presented as its source. In the paper, Figure 1 is on page two with a separate caption below it.

The self transform is drawn with column vectors; the code uses the equivalent transposed row-vector convention. R0 uses eight typed transforms. R1 uses a shared directional backbone base and a residual gated by positive training occurrence. The no-message control substitutes receiver states while retaining degrees. Activity loss is equal-component standardized squared error. Pair supervision adds lambda times equal-component squared error of `(z1-z0) - (y1-y0)/s_train`, with shared weights; grouped development selects lambda. Outer held-out labels enter evaluation only; inner validation selects models. Exact equations are implemented in [architecture.py](../../v3/architecture.py) and [engine.py](../../v3/engine.py).

| Figure | Content |
|---|---|
| [Main support](figures/main_support.pdf) | Directed relation counts and training gradients |
| [Main activity](figures/main_activity.pdf) | Ten-seed errors, ensembles and calibration |
| [Main B2/ranking](figures/main_pair_ranking.pdf) | Tie-aware APP ranking and measured chemistry differences |
| `appendix_*.pdf` | Eight existing supporting figures: learning, calibration, seeds, ranking, B2, B3 and source reconciliation |

Each result figure has PDF, SVG and PNG exports. [Plotting code](../../v3/figures.py) is preserved as executed; it reads original run paths and copies figures into the historical local manuscript directory, which is not included in GitHub. It requires omitted row-level inputs for some panels and is not a standalone reproduction command for this aggregate checkout. Aggregate plotted marks are in `figure_data/`; full marks and predictions remain in the local scientific archive. Source paths/hashes are recorded in [figure provenance](../../results/v3/figure_provenance.json). Figure assets visualize recorded measurements/predictions; they are not new measurements.
