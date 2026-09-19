"""Claim-specific primary citation audit; no imported statistical guarantee."""
from common import *
import re
P=ROOT/'papers/interaction_recoverability_iclr2027/v6';S=RUN/'sources/citations';D=ROOT/'docs/sirna_gnn_empirical/v3';bib='\n'.join((P/n).read_text() for n in ['references.bib','references_new.bib']);keys=re.findall(r'@\w+\{([^,]+),',bib);tex='\n'.join(p.read_text() for p in P.rglob('*.tex') if '_renderdeps' not in p.parts);cited=set()
for q in re.findall(r'\\cite\w*\{([^}]+)\}',tex):cited.update(k.strip() for k in q.split(','))
assert cited<=set(keys),cited-set(keys)
scopes={
'huesken2005':'Early neural siRNA design and assay provenance; not a new campaign cohort.',
'huesken2005correction':'Primary correction distinguishes reported predictions from experimental supplement values; no data added here.',
'huesken2006correction':'Primary correction of filtering/assay description; no historical numerical values imported as new evidence.',
'gneiting2007':'Distribution scoring is a distinct target from point-prediction squared error; no GNN risk theorem invoked.',
'elbashir2001':'Mammalian-cell RNAi motivation only; not validation of this predictor.',
'reynolds2004':'Sequence-design criteria under the reported screen, not universal potency guarantees.',
'uitei2004':'Sequence preferences in specified mammalian/chick assays; no assumption of universal transport.',
'khvorova2003':'Strand bias and terminal duplex stability motivate strand-specific inputs.',
'schwarz2003':'Asymmetric RISC assembly is biological context, not a measured mechanism of our model.',
'jackson2003':'Efficacy and transcript off-target specificity differ; no off-target score inferred here.',
'birmingham2006':'Seed-match association concerns off-targets, distinct from activity regression.',
'allerson2005':'Fully modified motifs can alter potency/stability in their assays; does not isolate the current B2 bundle.',
'bramsen2009':'Primary Table1 and SupplementalTable1 supply exact chemical states and measured eGFP B3 panel; methods and dependence checked.',
'nair2014':'GalNAc hepatocyte delivery is contextual and cannot establish our cohort transport or clinical benefit.',
'gilmer2017':'MPNN framework is a graph-prediction precedent; no exact reproduction or new architecture claimed.',
'battaglia2018':'Relational inductive-bias framework is conceptual context.',
'xu2019':'Graph expressivity does not certify the ordered finite-data predictor; no expressivity theorem applied.',
'zaheer2017':'Permutation-invariant set pooling differs from the actual ordered readout; theorem not transferred.',
'cmsirna2026':'Database provenance and chemistry annotation; primary patent checking remains separate.',
'davis2025':'Primary full text/workbook for separate S1 quarantine and S7 endpoint evidence; units/orientation not guessed.',
'meg2026':'Primary abstract establishes chemistry-aware graph precedent; exact code claims refer to pinned inspected source and failed invocation.',
'ridge1970':'Quadratic regularization precedent; own weighted objective stated, no original iid risk result imported.',
'extratrees2006':'Randomized feature/cut-point tree ensemble; actual installed settings specified independently.',
'numpy2020':'Numerical software attribution, installed version separately recorded.',
'scipy2020':'Scientific software attribution, not numerical certification.',
'matplotlib2007':'Plotting implementation attribution; exact mark-level exports retained.',
'adam2015':'Adaptive-moment optimizer precedent, not convergence under current nonconvex grouped data.',
'adamw2019':'Decoupled weight decay distinguished from data gradient; exact update library settings recorded.',
'layernorm2016':'Per-state normalization precedent, actual coordinate axis and epsilon explicit.',
'rgcn2018':'Relation-specific transforms precedent, different task from biological outcome prediction.',
'gcn2017':'Different graph normalization baseline context; not reproduced here.',
'graphsage2017':'Neighbor aggregation/sampling precedent; actual current dense aggregation supplied.',
'kim2014':'Local token convolution precedent, not a reproduced sentence classifier.',
'ensemble2017':'Member prediction variation is distinct from biological/sample uncertainty; no coverage claim imported.',
'torch2019':'Executed model/autograd software attribution.',
'ensi2025':'Primary abstract and pinned official implementation establish geometric chemistry-aware precedent; missing PDB/Rosetta prevents exact matched fit.',
'rf2001':'Tree ensemble precedent, not identification of ExtraTrees with a reproduced original random forest.',
'leakage2023':'Leakage and dependence safeguards motivate explicit split roles; no claim all biological dependence is removed.',
'bootstrap1979':'Resampling precedent; current component scheme defined directly, no iid coverage theorem transferred.',
'dropout2014':'Stochastic masking precedent, actual .1 probability and evaluation convention stated.',
'sklearn2011':'Executed classical estimators/software attribution.',
'wilds2021':'Distribution-shift evaluation context; APP/S7 are not independent new studies.',
'oligogym2025':'Standardized oligonucleotide datasets/splits motivate explicit data roles; not a matched reproduced benchmark.'}
audit=[]
for key in keys:
 rec=json.loads((S/(key+'.json')).read_text());paths=[];ep=S/(key+'.europepmc.json')
 if ep.exists() and json.loads(ep.read_text()).get('resultList',{}).get('result'):paths.append(ep);access='Primary author abstract and deposited bibliographic identity via Europe PMC'
 elif key in ['gneiting2007','ridge1970','leakage2023']:paths=[S/(key+'.pdf'),S/(key+'.txt')];access='Primary paper PDF, relevant opening/model text inspected'
 elif key=='bootstrap1979':paths=[S/'bootstrap1979.primary.pdf',S/'bootstrap1979.primary.txt'];access='Original primary paper copy, opening resampling definition and bibliographic header inspected'
 else:paths=[S/(key+'.html'),S/(key+'.txt')];access='Primary publisher/author/conference text or official software citation page inspected'
 assert all(p.exists() for p in paths),(key,paths)
 assert key in scopes,key
 audit.append(dict(key=key,cited=key in cited,primary_url=rec['primary_url'],access_scope=access,current_claim=scopes[key],evidence=[dict(path=str(p.relative_to(RUN)),sha256=sha(p)) for p in paths],new_published_guarantee_invoked=False))
write(RUN/'citation_audit.json',audit)
text='# Current citation audit\n\n'+f'{len(keys)} bibliography entries: {len(keys)-2} research/method/software papers and two primary corrections. All are cited in the current empirical manuscript. No unrelated theory or uncited padding is included. Primary access failures (403, redirect or bot responses) were retained and resolved through the named primary-paper/author-abstract alternatives. HTTP 200 alone was not accepted as article verification.\n\nThe original HTML-to-text helper lacked an optional parser in the core environment; standard-library extraction completed the administrative read without changing the fitting environment. The preserved older citation audit was inspected but did not substitute for the current source/claim checks. Detailed model comparisons use the pinned code. Abstract-only evidence supports narrow contextual claims, not unseen published methods or numerical superiority. No published theorem is invoked for this GNN.\n\n| Key | Primary evidence | Verified use and limit |\n|---|---|---|\n'
for q in audit:text+=f"| {q['key']} | [{q['access_scope']}]({q['primary_url']}) | {q['current_claim']} |\n"
(D/'citation_audit.md').write_text(text);print(len(audit),'entries',len(cited),'cited')
