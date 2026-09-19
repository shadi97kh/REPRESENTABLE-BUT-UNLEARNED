"""Source-first independent chemical/assay reconciliation. No model predictions read."""
from common import *
import re,html,xml.etree.ElementTree as ET,collections,ast
V3=ROOT/'runs/sirna_gnn_empirical/v3-20260914T013314Z';out=RUN/'adjudication';out.mkdir(exist_ok=True)
xml=V3/'sources/PMC2685080.xml';workbook=V3/'sources/PMC2685080/gkp106_nar-02440-c-2008-File010.xls'
root=ET.parse(xml).getroot();table=root.find('.//table-wrap[@id="T1"]')
alias=json.loads((V3/'b3/mapping_rules.json').read_text())['aliases'];alias['AP']='2-Aminopropyl_sugar' # Table1 explicitly defines 2-prime aminopropyl.
# Independent parser: replace chemical tags by lexical tokens and then tokenize flat text.
states={};stateaudit=[]
for tr in table.findall('.//tr'):
 c=tr.findall('td')
 if len(c)!=6:continue
 name,role,kind,level,viability=[''.join(x.itertext()).strip() for x in c[:5]]
 markup=ET.tostring(c[5],encoding='unicode');marked=re.sub(r'<(?:sub|sup)>(.*?)</(?:sub|sup)>',lambda m:'['+re.sub('<[^>]+>','',m[1])+']',markup)
 plain=re.sub(r'\s+','',html.unescape(re.sub('<[^>]+>','',marked)));tokens=re.findall(r'([ACGUT])((?:\[[^\]]+\])*)',plain)
 residue=re.sub(r'[ACGUT](?:\[[^\]]+\])*','',plain);residue=re.sub(r"[\s′'53–-]",'',residue)
 mods=[re.findall(r'\[([^\]]+)\]',x[1]) for x in tokens];unknown=sorted({m for z in mods for m in z if m not in alias})
 seq=''.join(x[0] for x in tokens);chem=[sorted(alias[m] for m in z if m in alias) for z in mods]
 state=dict(name=name,role=role,sequence=seq,chemistry=chem,chemical_kind=kind,primary_mean=level,primary_viability=viability,markup=markup,unresolved=unknown+[residue] if residue else unknown,printed_orientation='Printed oligonucleotide order and AS/SS role retained; Table1 does not repeat terminal direction labels per row',length=len(seq),backbone='not individually position-specified; no resolved stereochemistry inferred',terminal='no separate terminal annotation in Table1')
 states[name]=state;stateaudit.append(state)
write(out/'primary_states.json',states)
legacy=json.loads((V3/'b3/primary_strand_states.json').read_text());checks=[]
for name,s in legacy.items():
 q=states[name];ok=q['sequence']==s['sequence'] and q['chemistry']==[sorted(n['mods']) for n in s['nodes']];assert ok,name;checks.append(dict(name=name,independent_parser_agrees=ok))
export(pd.DataFrame(checks),'adjudication/independent_parser_checks.csv')
def key(seq,mods):return (seq,tuple(tuple(sorted(z)) for z in mods))
lookup={role:collections.defaultdict(list) for role in ['AS','SS']}
for name,s in states.items():
 if not s['unresolved']:lookup[s['role']][key(s['sequence'],s['chemistry'])].append(name)
# Read the named assay fields; direct xlrd accesses verify original row and column cells.
import xlrd
book=xlrd.open_workbook(workbook);sheet=book.sheet_by_name('Sheet1');headers=sheet.row_values(0);ix={str(n).strip():i for i,n in enumerate(headers)}
primary=[];byid=collections.defaultdict(list)
for j in range(1,sheet.nrows):
 a,s=str(sheet.cell_value(j,0)).strip(),str(sheet.cell_value(j,1)).strip()
 if a not in states or s not in states:continue
 try:y=float(sheet.cell_value(j,ix['Relative eGFP levels']));sd=float(sheet.cell_value(j,ix['SD (eGFP levels)']));v=float(sheet.cell_value(j,ix['Relative Viability']))
 except ValueError:continue
 z=dict(AS=a,SS=s,primary_sheet='Sheet1',primary_row=j+1,primary_column='D',remaining=y,primary_activity=1-y,primary_sd=sd,primary_viability=v,pair_cell=sheet.cell_value(j,2));primary.append(z);byid[a,s].append(z)
pandas=pd.read_excel(workbook,sheet_name='Sheet1')
for z in primary:
 j=z['primary_row']-2;assert pandas.iloc[j,0]==z['AS'] and pandas.iloc[j,1]==z['SS'];assert float(pandas.iloc[j,3])==z['remaining']
# Primary marginal means in Table1 corroborate the named eGFP column, independent of release agreement.
marginal=[]
for n,s in states.items():
 pair=(n,'W207') if s['role']=='AS' else ('W053',n);z=byid.get(pair,[]);mt=re.match(r'([\d.]+)\s*±\s*([\d.]+)',s['primary_mean'])
 if len(z)==1 and mt:
  residual=z[0]['remaining']-float(mt[1]);marginal.append(dict(strand=n,role=s['role'],printed_mean=float(mt[1]),supplement_mean=z[0]['remaining'],residual=residual,within_printed_rounding=abs(residual)<=.005000001))
export(pd.DataFrame(marginal),'adjudication/primary_marginal_checks.csv')
# Aliases are copied from the original documented parser, then raw workbook annotations parsed independently.
aliasnode=next(n for n in ast.parse((ROOT/'sirna_gnn_empirical/v1/common.py').read_text()).body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='ALIASES' for t in n.targets));release_alias=ast.literal_eval(aliasnode.value)
def release_state(seq,mm,pp):
 seq=str(seq).strip();mods=[[] for _ in seq]
 if str(mm).strip() not in ['0','0.0']:
  types=str(mm).split('*');positions=str(pp).split('*');assert len(types)==len(positions)
  for m,positions in zip(types,positions):
   m=release_alias.get(m.strip(),m.strip())
   for token in positions.split(','):mods[int(float(token.strip()))-1].append(m)
 return seq,[sorted(set(x)) for x in mods]
admitted={x['record_id']:x for x in rows_all() if x['source_family']=='19282453'}
oldmap=pd.read_csv(V3/'b3/record_identity_mapping.csv').fillna('').set_index('record_id');rows=[]
for f in ['train_88_1.xlsx','valid_88_1.xlsx','test_88.xlsx','patent.xlsx']:
 p=OLD/'raw/ensi/ENsiRNA-mod/dataset'/f;df=pd.read_excel(p).fillna('')
 for j,x in df.iterrows():
  if str(x.get('source','')).strip()!='19282453':continue
  rid=f'ensi:{f}:{j+2}';a,am=release_state(x['anti raw seq'],x['anti mod'],x['anti pos']);s,sm=release_state(x['sense raw seq'],x['sense mod'],x['sense pos']);aa=lookup['AS'][key(a,am)];ss=lookup['SS'][key(s,sm)]
  if rid in admitted:
   frozen=admitted[rid]
   for st,seq,mods in [('guide',a,am),('passenger',s,sm)]:assert key(seq,mods)==key(frozen[st]['sequence'],[[m for m in n['mods'] if m!='unmodified_as_annotated'] for n in frozen[st]['nodes']]),rid
  primaryrows=byid.get((aa[0],ss[0]),[]) if len(aa)==len(ss)==1 else []
  category='unmatched_reported_chemical_identity' if not aa or not ss else 'ambiguous_chemical_identity' if len(aa)!=1 or len(ss)!=1 else 'ambiguous_or_missing_primary_measurement' if len(primaryrows)!=1 else 'dose_unresolved' if float(x['cc'])!=10 else 'unique_reported_identity_primary_assay_reanchored'
  z=dict(record_id=rid,admitted_historical=rid in admitted,release_workbook=str(p.relative_to(ROOT)),release_sheet='Sheet1',release_row=j+2,release_persistent_id=x['ID'],release_source=str(x['source']),release_raw_pct=float(x['PCT']),release_activity=float(x['PCT'])/100,release_concentration=float(x['cc']),release_concentration_units='unverified in released workbook; primary reference documents 10 nM',release_assay_metadata='not individually provided; no row-level assay provenance key',guide_sequence=a,passenger_sequence=s,guide_length=len(a),passenger_length=len(s),guide_positional_chemistry=json.dumps(am),passenger_positional_chemistry=json.dumps(sm),release_raw_guide_mod=x['anti mod'],release_raw_guide_positions=x['anti pos'],release_raw_passenger_mod=x['sense mod'],release_raw_passenger_positions=x['sense pos'],AS_candidates='|'.join(aa),SS_candidates='|'.join(ss),primary_pair_cardinality=len(primaryrows),category=category,primary_assay='HeLa-eGFP; 10 nM; 72 h; INTERFERin; confocal eGFP; mismatch-siRNA normalized',primary_replication='triplicate assays repeated twice; aggregate SD provided; covariance and run-level observations unavailable',outcome_convention='primary relative eGFP fraction -> inhibition=1-remaining; released PCT -> PCT/100; no clipping',identity_selection='sequence/positional chemistry/source/dose only; no outcome or prediction search',orientation='AS and SS both retain printed oligonucleotide order; no swap/reversal',unresolved_chemical_fields='per-bond stereochemistry and explicit full terminal/backbone specification absent from release; literal reported identity, not certified complete molecule',adjudication='A separately named primary-assay target is justified on uniquely reported identities; discrepancy alone does not identify which historical label/mapping is wrong.')
  if len(primaryrows)==1:
   z.update(primaryrows[0]);z['residual']=z['release_activity']-z['primary_activity'];z['absolute_difference']=abs(z['residual'])
  rows.append(z)
df=pd.DataFrame(rows)
unique=df.category.eq('unique_reported_identity_primary_assay_reanchored'); multiplicity=df[unique & df.admitted_historical].groupby(['AS','SS']).size();ambiguous=set(multiplicity[multiplicity>1].index); df['released_identity_cardinality']=[int(multiplicity.get((a,b),0)) for a,b in zip(df.AS,df.SS)]; df.loc[unique & df.released_identity_cardinality.gt(1),'category']='ambiguous_multiple_released_rows_per_primary_identity'
export(df,'adjudication/row_reconciliation.csv')
# Audit all above-threshold cases and deterministic SHA-ordered agreeing controls.
legacydiff=pd.read_csv(V3/'b3/primary_reconciliation.csv');bad=set(legacydiff.loc[legacydiff.absolute_difference>.05,'record_id']);assert len(bad)==299
verified=df[df.admitted_historical & df.category.eq('unique_reported_identity_primary_assay_reanchored')]
for rid in bad:
 q=df[df.record_id==rid].iloc[0];old=legacydiff[legacydiff.record_id==rid].iloc[0];assert q.AS==old.AS and q.SS==old.SS and abs(q.absolute_difference-old.absolute_difference)<1e-12
controls=verified[~verified.record_id.isin(bad)].copy();controls['order']=controls.record_id.map(lambda s:hashlib.sha256(s.encode()).hexdigest());sample=controls.sort_values('order').head(32).record_id
export(df[df.record_id.isin(bad|set(sample))],'adjudication/direct_cell_review.csv')
# Fresh primary-reanchored subset, with unresolved source identities quarantined. No original values overwritten.
changes=verified[['record_id','release_activity','primary_activity','AS','SS','primary_row','primary_column','category']].copy();export(changes,'adjudication/primary_reanchored_targets.csv')
summary=dict(raw_source_records=len(df),historically_admitted_source_records=int(df.admitted_historical.sum()),categories_all=df.category.value_counts().to_dict(),categories_admitted=df[df.admitted_historical].category.value_counts().to_dict(),legacy_matched_count=1604,legacy_above_005=299,legacy_discrepancies_independently_reproduced=True,direct_cell_checks=len(bad)+len(sample),independent_primary_state_checks=len(checks),primary_marginal_checks=len(marginal),primary_marginal_rounding_failures=sum(not x['within_printed_rounding'] for x in marginal),primary_reanchored_source_rows=len(verified),primary_reanchored_other_source_rows=1003,primary_reanchored_changes_above_005=int((verified.absolute_difference>.05).sum()),primary_reanchored_max_residual=float(verified.absolute_difference.max()),new_AP_alias='Primary AP explicitly means 2-prime aminopropyl; matches sugar alias only, never base alias',decision='Create primary-reanchored sensitivity target on uniquely reported identities, quarantine remaining Bramsen rows. This is not a certified repair of ENsiRNA labels or evidence of incorrect biological measurements; missing full assay/chemical traceability remains explicit.',forbidden_searches='No transformations optimized against values, no model predictions read, no clipping, no row-shift/permutation fitted to agreement',raw_sha256={str(p.relative_to(ROOT)):sha(p) for p in [xml,workbook]},shared_reference='No endpoint covariance inferred from SD. Distinct rectangles are dependent.')
write(out/'summary.json',summary);print(json.dumps(summary,indent=2))
