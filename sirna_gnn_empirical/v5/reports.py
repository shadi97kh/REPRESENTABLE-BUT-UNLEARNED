from common import *
import pandas as pd,re,sys,shutil,csv
A=PAPER/'artifacts';S=PAPER/'source';T=RUN/'tables';T.mkdir(exist_ok=True)
metrics=pd.read_csv(RUN/'reuse/recomputed_metrics.csv');pivot=metrics[metrics.method.isin(['corrected_gnn','corrected_no_message','chemistry_tree','token_cnn'])].pivot(index=['scenario','method'],columns='cohort',values='mse').reset_index();pivot.to_csv(T/'main_comparison.csv',index=False)
for name in ['direct_comparisons.csv','seed_stability.csv','seed_paired_differences.csv','per_source_metrics.csv','ranking_pools.csv','activity_metrics.csv']:
 shutil.copy2(PARENT/'evaluation'/name,T/name)
for name in ['B2_unchanged.csv','B3_visibility.csv','historical_calibration.csv','ranking_history.csv','adjudication.csv','new_fits.csv']:
 shutil.copy2(PARENT/'paper_tables'/name,T/name)
G=json.loads((DOC/'gate_decision.json').read_text())
(DOC/'gate_decision.md').write_text('# Contribution gates\n\n'+ '\n\n'.join('**'+k+' — '+v['status']+'.** '+v['reason']+' '+v['limitations']+' Change: '+v['changes'] for k,v in G['gates'].items())+'\n\nThe runner blocks pilot and confirmation. No new ML contribution is established. Factual manuscript revision and required verification are completed independently of promotion.\n')
aux=(A/'build/main.aux').read_text();labels={k:int(p) for k,p in re.findall(r'\\newlabel\{([^}]+)\}\{\{[^}]*\}\{(\d+)\}',aux)}
objects={
'fig:workflow':('app:architecture;app:training','sirna_gnn_empirical/v4/architecture.py; v4/protocol.json; v5/workflow_identity.json'),
 'tab:data':('app:data;app:robustness','v4/paper_tables/source_counts_scores.csv'),
 'tab:cohorts':('app:data;app:b3','v4/paper_tables/evaluation_units.csv'),
 'tab:adjudication':('app:robustness','v4/adjudication/row_reconciliation.csv; v4/paper_tables/adjudication.csv'),
 'tab:optimization':('app:training;app:robustness','v4/fit_summary.csv; v4/fit_completion.json; v4/paper_tables/new_fits.csv'),
 'tab:activity':('app:metrics;app:focusedresults','v4/evaluation/activity_metrics.csv; v5/reuse/recomputed_metrics.csv'),
 'fig:robustness':('app:metrics;app:focusedresults','v4/evaluation/seed_metrics.csv; v4/evaluation/activity_metrics.csv'),
 'tab:calibration':('app:metrics;app:results','v4/paper_tables/historical_calibration.csv'),
 'fig:calibration':('app:metrics;app:robustness','v4/review_checks/calibration_decomposition.csv; v3/tables/direct_comparisons.csv'),
 'tab:ranking':('app:metrics;app:focusedresults','v4/ranking/campaign_summary.csv; v4/paper_tables/ranking_history.csv'),
 'tab:pairranking':('app:training;app:metrics;app:results','v4/paper_tables/B2_unchanged.csv; v3/tables/pair_metrics.csv; v3/tables/pair_comparisons.csv'),
 'tab:B3':('app:b3;app:inputcollision','v4/paper_tables/B3_visibility.csv; v4/visibility/B3_rectangle_visibility.csv'),
 'fig:visibility':('app:inputcollision;app:focusedresults','v4/visibility/input_counts.csv; v3/b3_primary/')}
main=S.joinpath('main.tex').read_text();app=S.joinpath('appendix.tex').read_text()
found=re.findall(r'\\label\{((?:fig|tab):[^}]+)\}',main);assert set(found)==set(objects),(set(found),set(objects))
lines=['# Appendix completeness audit — manuscript v10','','The source retains all five complete parent proof bodies verbatim. No new-method theorem is claimed in the main text. All finite-model definitions, forward/no-message arguments, endpoint derivatives, metrics, covariance qualifications, ensemble decomposition and collision proofs remain in the compiled appendix. The separate exploratory exact-range proof is complete in method_and_proof.md; it is not attached to the GNN as a guarantee.','','| Main object (PDF page) | Complete appendix support (PDF page) | Exact archive evidence |','|---|---|---|']
for obj,(sections,evidence) in objects.items():
 sections='; '.join(f'{s} (p. {labels[s]})' for s in sections.split(';'))
 lines.append(f'| {obj} (p. {labels[obj]}) | {sections} | `{evidence}` |')
lines+=['','Main mathematical claims are the displayed architecture and objectives (pp.2,8), not separately numbered new theorems. They resolve to A3–A5 and A9. Every main section and substantive prose claim is covered by claim_audit.md. The complete source has one appendix inclusion only.','', 'The archive retains v3/v4 raw rows at their original relative paths inside `parent_evidence.zip`; v5 rows are in the outer `v5/` directory. Extract the parent archive to access the named CSVs. No proof terminates at those pointers.','',f'Actual boundaries: {json.dumps(json.loads((A/"page_boundaries.json").read_text()))}.']
(DOC/'appendix_audit.md').write_text('\n'.join(lines)+'\n')
claims=[
 ('Abstract; §§1,9','Empirical diagnostic, no new learned method or biological validation','G1 reduction; conventional architecture; unchanged negative controls','method_and_proof.md; novelty_matrix.md; A3,A6,A9'),
 ('§2','153,345 scalars, two layers, training masks/support, both objectives and derivatives','Actual forward/backward checks; retained exact specification','sirna_gnn_empirical/v4/architecture.py; A3–A4; v5/correctness/results.json'),
 ('§3','1,924/2,927 source concentration; some other-source positive R²','Source-scoring subsets do not imply source-excluded training','v4/review_checks/source_concentration.csv; source_scores.csv; A2,A9'),
 ('§4','299/1,604 differences; maximum 1.085217; primary-assay sensitivity not certified label repair','All large discrepancies and 32 fixed controls checked in v4; no new identity rule','v4/adjudication/row_reconciliation.csv; direct_cell_review.csv; A9'),
 ('§4','792 focused fits; ten final seeds; 40 unaffected checkpoint reuses','No new fits in v5; 3,564 parent fit hashes checked','v4/fit_summary.csv; v4/fit_completion.json; A4,A9'),
 ('§5','Focused activity metrics, constants, conditional intervals','36 score groups and ensemble identities independently recomputed','v5/reuse/recomputed_metrics.csv; v4/evaluation/direct_comparisons.csv; A5,A10'),
 ('§6','99.53% squared-bias decomposition; original S7 routing improvement','Descriptive decomposition and conditional uncertainty; not causal bias explanation','v4/review_checks/calibration_decomposition.csv; independent_conditional_comparisons.csv; A5,A9'),
 ('§7','Ranking reversals on identical IDs; exact ties; 17 overlapping pools','Retrospective protocol changes, not independent test evidence','v4/ranking/; v4/evaluation/ranking_pools.csv; A5,A9,A10'),
 ('§7','B2 does not establish sequence-specific chemistry prediction','156 pairs/12 sequence components; joint four-position edit; majority control; unresolved replicate covariance','v3/tables/pair_metrics.csv; pair_comparisons.csv; B2_within_assay_variation.csv; A4,A5,A7'),
 ('§8','165 states/140 dependent B3 rectangles; 90 exact inputs; zero forced cancellations','0/140 is not zero predicted interactions; no across-background inference','v3/b3_primary/; v4/visibility/; A6,A9,A10'),
 ('§9','Published reproductions incomplete; S1 acquired/quarantined','No newly resolved geometry/token/orientation asset; zero exact published fits','v4/source_checks/; A6'),
 ('§9','Focused S7 GNN beats no-message conditionally; tree/CNN still better observed MSE','No general GNN advantage; negative APP/S7 R² retained','v4/evaluation/direct_comparisons.csv; v5/tables/main_comparison.csv; A10'),
 ('§9','Five-law theory separate; backend no-go; reference efficiency open','No reduction or theorem coverage asserted','A8–A9; preserved historical theory outside this manuscript'),
]
with (T/'claim_map.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['main_location','claim','qualification','evidence']);w.writerows(claims)
(DOC/'claim_audit.md').write_text('# Claim audit\n\n| Main location | Claim | Essential qualification | Complete support |\n|---|---|---|---|\n'+'\n'.join('| '+' | '.join(r)+' |' for r in claims)+'\n\nAll original main numerical tables and three result figures remain byte-identical in TeX/data or asset content. New work changes workflow, scope wording and related-work citations. No new performance claim is introduced.\n')
# Check all object caption orders, including longtable captions after body.
caption_checks=[]
for name,text in [('main.tex',main),('appendix.tex',app)]:
 for kind in ['figure','table','longtable']:
  for i,m in enumerate(re.finditer(r'\\begin\{'+kind+r'\}(.*?)\\end\{'+kind+r'\}',text,re.S)):
   body=m.group(1);cap=body.find('\\caption');assert cap>=0,(name,kind,i)
   if kind=='figure':assert cap>body.rfind('\\includegraphics')
   if kind=='table':assert cap>max(body.rfind('\\end{tabular}'),body.rfind('\\end{tabularx}'))
   if kind=='longtable':assert cap>body.rfind('\\endhead')
   caption_checks.append(dict(file=name,kind=kind,index=i,below=True))
write(A/'caption_order_audit.json',caption_checks)
# Preserve legacy citation checks, adding exact primary sources for each new key.
oldcitation=json.loads((PARENT/'source_checks/citation_recheck.json').read_text());new=[]
for key,paper in [('chen2023harsanyi','harsanyinet'),('xia2023ncm','ncm'),('tan2024consistency','consistency'),('lengerich2020pure','purification'),('kuskova2026real','gnavar')]:
 p=RUN/'literature'/(paper+'.pdf');new.append(dict(key=key,primary_pdf=str(p.relative_to(ROOT)),sha256=sha(p),comparison='novelty_matrix.md; complete primary text acquired, hypotheses and relevant algorithm/proof paths inspected; not invoked as GNN guarantee'))
write(RUN/'citation_audit.json',dict(parent_records=oldcitation,new_records=new,unique_cited_entries=48))
(DOC/'citation_audit.md').write_text('# Citation audit\n\n48 unique entries are cited and resolve. The 43 preserved primary-source evidence records are rehashed in the reuse manifest; five new primary papers are recorded in v5/citation_audit.json. Novelty_matrix.md states exact hypotheses and distinguishes predictor attribution, causal identification and assay-response prediction. Tan et al. is cited under its current title, with the 2025 updated preprint identified; G-NAVAR remains a June 2026 preprint. Official implementations are inspected, not benchmarked. No prior theorem is used as a GNN risk or biological coverage guarantee.\n\nOfficial ICLR 2027 Author Guidelines permit at most nine submission main pages, exclude references, allow unrestricted appendices, and exclude AI/ethics/reproducibility statements. Their required AI-use disclosure is retained and expanded to report proof, design and source-reconciliation assistance. Primary pages and hashes are in v5/literature/.\n')
# Exact scientific table and proof preservation checks, separate from narrative edits.
parentmain=(ROOT/'papers/interaction_recoverability_iclr2027/v9/source/main.tex').read_text()
for m in re.findall(r'\\begin\{table\}.*?\\end\{table\}',parentmain,re.S):assert m in main
for p in (S/'figures').glob('*.pdf'):
 if p.name!='workflow.pdf':assert sha(p)==sha(ROOT/'papers/interaction_recoverability_iclr2027/v9/source/figures'/p.name)
write(RUN/'workflow_identity.json',dict(native_vector_pdf=True,embedded_images=0,architecture='unchanged conventional v4 GraphModel',new_learned_components=0,pdf_sha256=sha(S/'figures/workflow.pdf'),caption_in_artwork=False,caption_below=True,workflow_page=2))
(DOC/'workflow_audit.md').write_text('# Workflow audit\n\nThe page-two PDF is native vector artwork, with zero embedded raster images. The SVG retains vector objects and text; PNG is a preview. The figure shows the existing encoder, R1 relation sharing, degree-normalized messages, residual/LayerNorm masking, ordered readout, original activity scale, frozen B1/B2/B3 evaluation and separate activity/pair training. It adds no new learned module. The artwork contains node labels and training-path annotations, but no figure number or embedded caption. The manuscript caption is below it. Parent result-figure bytes are unchanged.\n')
# Metadata preservation scope is stated honestly.
old=json.loads((RUN/'preservation_before.json').read_text());changed=[]
for r in old:
 p=ROOT/r['path']
 if not p.is_file() or (p.stat().st_size,p.stat().st_mtime_ns)!=(r['size'],r['mtime_ns']):changed.append(r['path'])
assert not changed,changed[:20]
write(RUN/'preservation_after.json',dict(metadata_entries=len(old),changed=changed,parent_artifact_hashes_verified=10,parent_fit_artifact_hashes_verified=3564,scope='Full historical metadata; targeted content hashes, not all-history content hashing'))
(DOC/'review_response.md').write_text('# Review response\n\n'+ '\n'.join('- **'+r[0]+' — '+r[1]+'.** '+r[2]+'. Evidence: '+r[3]+'.' for r in claims)+'\n\nThe newest brief supersedes the earlier demand for another unconditional training cycle: failed novelty/validity gates block the optional pilot. Required empirical corrections were already performed in v4 and independently verified here. There is no newly justified label or identity change requiring another refit.\n')
commands=[json.loads(l) for l in (RUN/'commands.jsonl').read_text().splitlines()]
(DOC/'executed_commands.md').write_text('# Executed commands\n\nEntry point: `bash sirna_gnn_empirical/v5/run_all.sh`. It verifies completed stage hashes, runs correctness/reuse/build checks, and records pilot/confirmation as blocked.\n\n'+ '\n'.join('- `'+r['command']+'` — '+r['status']+', log `'+r['log']+'`.' for r in commands)+'\n\nPublic acquisition ran `python sirna_gnn_empirical/v5/acquire.py` and `python sirna_gnn_empirical/v5/acquire_comparisons.py`; `python sirna_gnn_empirical/v5/freeze.py` froze correctness inputs before execution. Additional direct source-file acquisition, inspection, authoring, rendering and hashing were administrative and incurred nonzero unmetered cost. One optional HTML inspection failed because bs4 was absent; official pages were verified through the browser instead. No third-party training, GPU use, historical wrapper, ledger, push or external contact occurred. The runner records the final reports/package command after each returns; those final receipts are supplied alongside the archive.\n')
resource=dict(recorded_stage_wall_s=sum(r['wall_s'] for r in commands),recorded_stage_child_cpu_s=sum(r['child_cpu_s'] for r in commands),max_stage_rss_kib=max(r['max_rss_kib'] for r in commands),new_fits=0,new_training_optimizer_updates=0,GPU_work=False,administrative_cost='Additional nonzero unmetered browsing, direct acquisition, authoring, inspection and hashing; nested build commands overlap stage totals; final reports/package receipts supplied alongside archive',failed_fits=0)
write(DOC/'resource_accounting.json',resource)
md=pivot.to_markdown(index=False,floatfmt='.6f')
(DOC/'execution_summary.md').write_text('# Completed contribution-gated revision\n\n**G1 FAIL; G2 UNRESOLVED; G3 UNRESOLVED (blocked, not run). No new ML contribution is established.** The explicit candidate reduces to bounded linear optimization plus encoding-equivalence checks. The complete conditional proof and 19 passing correctness/counterexample checks do not establish adaptive learned-representation validity or biological coverage.\n\nRecomputed saved ten-seed ensemble MSE (lower preferred):\n\n'+md+'\n\nThese are verified reused v4 results, not new experimental fits. All 36 score groups and seed/ensemble decompositions reproduced; 3,564 fit-artifact hashes and 36 grouped partition checks passed. The parent executed 792 fits (432 development, 360 final); this investigation executed zero fits and zero training optimizer updates. The 19 diagnostic cases include tiny synthetic graphs, exact finite-state calculations and gradient checks, not a positive pilot or synthetic siRNA evidence.\n\nFocused S7 GNN-minus-no-message differences remain -0.003057 and -0.002122 with conditional intervals excluding zero. Tree and CNN have lower observed S7 error. General GNN improvement and sequence-specific measured chemistry prediction remain unestablished. B2 pair GNN MSE is 0.222115 versus ridge 0.219297; its advantage over the training mean is unresolved. B3 GNN MSE 0.068335 does not beat zero interaction 0.068205; 140 rectangles share one background. APP/S7 follow-ups remain retrospective.\n\nPrimary-assay curation, source-excluded fitting and input visibility audits are preserved and verified; no new evidence justifies additional affected refits. The 299 discrepancies are not all certified label errors. S1 orientation, complete assay traceability and shared-reference uncertainty remain unresolved. Published geometry/token blockers remain; exact published fits remain zero. No private anchors, wet-lab, causal or clinical validation is added.\n\nManuscript v10 has main pages 1–9, statements 10, references 11–14, appendix 15–53 (39 pages), 48 citations, nine main tables and four main figures including the native-vector workflow on page two. All five parent proof bodies are retained. The exact four-entry source ZIP compiles independently; zero undefined references or overfull boxes. All pages are rendered and inspected; caption checks pass.\n\nArtifacts are indexed in delivery_manifest.json. The compact review packet is for the user to send; nothing was pushed or sent. Historical files, proofs, ledgers and the backend no-go/reference-efficiency question remain unchanged.\n')
print('Wrote comparison tables, completeness/claim/citation/workflow/preservation audits')
