from common import *
import shutil,re,subprocess,sys
old=ROOT/'papers/interaction_recoverability_iclr2027/v9';S=PAPER/'source';A=PAPER/'artifacts';A.mkdir(exist_ok=True)
if not S.exists():shutil.copytree(old/'source',S)
if not (PAPER/'typesetting').exists():shutil.copytree(old/'typesetting',PAPER/'typesetting')
subprocess.run([sys.executable,str(HERE/'workflow.py')],check=True)
shutil.copy2(A/'figures/workflow.pdf',S/'figures/workflow.pdf')
main=(old/'source/main.tex').read_text();app=(old/'source/appendix.tex').read_text();bib=(old/'source/references.bib').read_text()
main=main.replace('The artwork uses column vectors; Appendix~\\ref{app:architecture} gives equivalent row-vector equations.','The artwork and Appendix~\\ref{app:architecture} use row-vector messages.')
# Resolve the abstract's ambiguous scope: historical routing versus new sensitivity controls.
main=main.replace('Relation sharing gives a narrow conditional error reduction on Davis S7, without establishing a message-passing advantage.','The historical routing comparison gives a narrow conditional S7 error reduction, distinct from the focused message-passing comparisons.')
needle='The contribution is the controlled support comparison and the source, calibration and measured-effect diagnostics that qualify it.'
main=main.replace(needle,needle+' Predictor interaction attribution, compatible-model identification and support-dependent recoverability are distinct prior problems \\citep{chen2023harsanyi,xia2023ncm,tan2024consistency,lengerich2020pure,kuskova2026real}. Our finite input diagnostics do not establish a new learned interaction method.')
main=main.replace('AI tools assisted with code, analysis, source inspection, writing, figure construction and manuscript checks.','AI tools assisted with research framing, experiment design, source and chemistry reconciliation, mathematical formulations and proof writing, implementation, result interpretation, literature comparison, prose, figure construction and manuscript checks.')
# Full relevant mathematical content remains unchanged; research attempt is not grafted onto the empirical paper.
needle='\\section{'
# Add a scoped paragraph in the existing prior-work appendix before the next subsection.
anchor='\\label{app:published}'
if anchor not in app:
 anchor='\\label{app:priorwork}'
if anchor in app:
 at=app.index('\n',app.index(anchor))
 addition=r'''
\paragraph{Interaction attribution and identification.}
HarsanyiNet computes attributions of its specially structured predictor under its receptive-field and masking requirements \citep{chen2023harsanyi}; this network has neither those units nor their attribution theorem. Neural causal identification uses graph-consistent structural models and matched distributions \citep{xia2023ncm}; molecular adjacency supplies neither a causal diagram nor those distributional constraints. Continuous neural partial-identification results require bounded Lipschitz mechanisms, latent-support regularity and controlled approximation and matching schedules \citep{tan2024consistency}. Those assumptions are not verified here, and local neural extrema would not provide certified outer bounds. Functional-ANOVA purification fixes a weighting-dependent decomposition of a predictor \citep{lengerich2020pure}; it does not supply an unobserved joint assay mean. The G-NAVAR preprint's population result assumes fixed direct-sum edge spaces, positive product support and regularity, with additional HOFD and faithfulness conditions in its shared-modulator formulation \citep{kuskova2026real}. A covariance-rank diagnostic is not sufficient identification. None of these results is invoked as a guarantee for this siRNA learner. Four forward passes, encoding-equivalence checks and optimization over a compatible finite response class do not by themselves establish a new ML contribution.
'''
 app=app[:at]+addition+app[at:]
else:raise RuntimeError('Locate existing prior-work section before insertion')
newbib=r'''
@inproceedings{chen2023harsanyi, title={HarsanyiNet: Computing Accurate Shapley Values in a Single Forward Propagation}, author={Chen, Lu and Lou, Siyu and Zhang, Keyan and Huang, Jin and Zhang, Quanshi}, booktitle={International Conference on Machine Learning}, pages={4804--4825}, year={2023}, url={https://proceedings.mlr.press/v202/chen23s.html}}
@inproceedings{xia2023ncm, title={Neural Causal Models for Counterfactual Identification and Estimation}, author={Xia, Kevin and Pan, Yushu and Bareinboim, Elias}, booktitle={International Conference on Learning Representations}, year={2023}, url={https://arxiv.org/abs/2210.00035}}
@inproceedings{tan2024consistency, title={Consistency of Neural Causal Partial Identification}, author={Tan, Jiyuan and Blanchet, Jose and Syrgkanis, Vasilis}, booktitle={Advances in Neural Information Processing Systems}, year={2024}, note={Updated arXiv version 3, June 2025}, url={https://arxiv.org/abs/2405.15673}}
@inproceedings{lengerich2020pure, title={Purifying Interaction Effects with the Functional ANOVA: An Efficient Algorithm for Recovering Identifiable Additive Models}, author={Lengerich, Benjamin and Tan, Sarah and Chang, Chun-Hao and Hooker, Giles and Caruana, Rich}, booktitle={International Conference on Artificial Intelligence and Statistics}, pages={2402--2412}, year={2020}, url={https://proceedings.mlr.press/v108/lengerich20a.html}}
@article{kuskova2026real, title={When Are Neural Interaction Discoveries Real? Identifiability, Recoverability, and a Pre-Fit Diagnostic}, author={Kuskova, Valentina and Zaytsev, Dmitry and Coppedge, Michael}, journal={arXiv preprint arXiv:2606.08390}, year={2026}, url={https://arxiv.org/abs/2606.08390}}
'''
(S/'main.tex').write_text(main);(S/'appendix.tex').write_text(app);(S/'references.bib').write_text(bib+newbib)
proofs=re.findall(r'\\begin\{proof\}.*?\\end\{proof\}',(old/'source/appendix.tex').read_text(),re.S)
assert all(p in app for p in proofs)
write(A/'authoring_proof_preservation.json',dict(parent='v9',complete_proofs_preserved=len(proofs),verbatim=True))
write(A/'manuscript_changes.json',dict(parent='v9',new='v10',scientific_numbers_changed=False,changes=['native vector workflow','clarified historical-routing abstract sentence','scoped primary prior-work comparison; no new-method claim'],proofs_preserved=True))
print('Prepared factual v10 manuscript; all parent proof bodies retained')
