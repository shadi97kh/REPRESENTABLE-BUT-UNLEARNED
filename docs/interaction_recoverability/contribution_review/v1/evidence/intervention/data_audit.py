"""Chemical identity and row provenance audit; labels are never fitted here."""
from __future__ import annotations
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
import math
from .accounting import snapshot, elapsed, sha, write_json


def main():
    start = snapshot()
    from rdkit import Chem, rdBase
    from rdkit.Chem.Scaffolds import MurckoScaffold
    root = Path('data/intervention_sources')
    manifest = json.loads((root / 'manifest-v1.json').read_text())
    cache = {}
    def identity(smiles):
        if smiles in cache:
            return cache[smiles]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            result = dict(valid=False, canonical_isomeric=None, canonical_connectivity=None,
                          scaffold=None, specified_stereocenters=None, unspecified_stereocenters=None)
        else:
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
            centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            result = dict(valid=True,
                          canonical_isomeric=Chem.MolToSmiles(mol, isomericSmiles=True),
                          canonical_connectivity=Chem.MolToSmiles(mol, isomericSmiles=False),
                          scaffold=MurckoScaffold.MurckoScaffoldSmiles(mol=mol),
                          specified_stereocenters=sum(c != '?' for _, c in centers),
                          unspecified_stereocenters=sum(c == '?' for _, c in centers))
        cache[smiles] = result
        return result
    allrows, summaries, rowsets = [], {}, {}
    for file in manifest['files']:
        assert sha(file['path']) == file['sha256']
        name = file['name']
        rows = list(csv.DictReader(Path(file['path']).open()))
        output = []
        for i, raw in enumerate(rows):
            smiles = raw.get('combined_mol_SMILES', raw.get('SMILES'))
            value = raw.get('expt_Hela', raw.get('Target', raw.get('TARGET')))
            row = dict(dataset=name, source_row_number=i+2, source_url=file['url'],
                       source_sha256=file['sha256'], raw_fields=raw,
                       raw_smiles=smiles, **identity(smiles),
                       measurement_kind='published_aggregated_delivery_endpoint',
                       raw_hela_label=value, hela_label_numeric=float(value),
                       raw_raw2647_label=raw.get('expt_Raw'),
                       components={k: {'raw': raw[k], **identity(raw[k])}
                                   for k in ('A_smiles','B_smiles','C_smiles') if k in raw},
                       assay_cell_line='HeLa', assay_source='paper/repository dataset identity; not a row field',
                       formulation=None, cargo=None, dose=None, timing=None, endpoint_scale=None,
                       individual_replicates=None, replicate_aggregation=None,
                       missing_assay_fields='not represented in these CSVs; no row-specific values imputed',
                       split_origin=('published_scaffold_'+name.split('_')[-1]) if name.startswith('transma') else None,
                       action_pair_group=None)
            if not math.isfinite(row['hela_label_numeric']):
                raise ValueError('nonfinite source label')
            output.append(row)
        allrows.extend(output); rowsets[name] = output
        groups = defaultdict(list)
        for r in output:
            if r['valid']: groups[r['canonical_isomeric']].append(r)
        conflicts = []
        for k, rr in groups.items():
            labels = [r['hela_label_numeric'] for r in rr]
            conflict = max(labels) - min(labels) > 1e-6
            for r in rr:
                r['canonical_duplicate_count'] = len(rr)
                r['conflicting_structure_labels'] = conflict
                r['future_pair_eligibility'] = 'unresolved_assay_metadata' if not conflict else 'exclude_unresolved_conflicting_labels'
            if conflict:
                conflicts.append(dict(canonical_isomeric=k, rows=[r['source_row_number'] for r in rr], labels=labels))
        summaries[name] = dict(rows=len(rows), raw_unique=len({r['raw_smiles'] for r in output}),
                               valid_rows=sum(r['valid'] for r in output),
                               canonical_isomeric_unique=len(groups),
                               canonical_connectivity_unique=len({r['canonical_connectivity'] for r in output if r['valid']}),
                               duplicate_groups=sum(len(rr)>1 for rr in groups.values()),
                               duplicate_excess_rows=sum(len(rr)-1 for rr in groups.values()),
                               conflicting_groups=conflicts,
                               rows_with_unspecified_stereocenters=sum(bool(r['unspecified_stereocenters']) for r in output),
                               component_counts={k:dict(raw=len({r['components'][k]['raw'] for r in output}),
                                                        canonical=len({r['components'][k]['canonical_isomeric'] for r in output}))
                                                 for k in output[0]['components']})
    overlaps = []
    names = list(rowsets)
    for i, a in enumerate(names):
        for b in names[i+1:]:
            ga, gb = defaultdict(list), defaultdict(list)
            for r in rowsets[a]: ga[r['canonical_isomeric']].append(r)
            for r in rowsets[b]: gb[r['canonical_isomeric']].append(r)
            common = (set(ga) & set(gb)) - {None}
            comparisons = [(k,ra,rb) for k in common for ra in ga[k] for rb in gb[k]]
            conflicts = [dict(canonical_isomeric=k, left_row=ra['source_row_number'],
                              right_row=rb['source_row_number'], left_label=ra['raw_hela_label'],
                              right_label=rb['raw_hela_label']) for k,ra,rb in comparisons
                         if abs(ra['hela_label_numeric']-rb['hela_label_numeric'])>1e-6]
            overlaps.append(dict(left=a,right=b,unique_canonical_overlap=len(common),
                                 compared_row_pairs=len(comparisons), conflicting_row_pairs=len(conflicts),
                                 conflicting_unique_structures=len({r['canonical_isomeric'] for r in conflicts}),
                                 conflicts=conflicts))
    tr,te = rowsets['transma_train'], rowsets['transma_test']
    scaffold_overlap = set(r['scaffold'] for r in tr)&set(r['scaffold'] for r in te)
    audit = dict(rdkit_version=rdBase.rdkitVersion, canonicalization='sanitized canonical isomeric SMILES; no salt stripping, tautomer/protonation normalization or stereochemical imputation',
                 datasets=summaries, cross_dataset_overlap=overlaps,
                 transma_scaffold_overlap=dict(count=len(scaffold_overlap), scaffolds=sorted(scaffold_overlap),
                                              includes_empty_acyclic_scaffold=('' in scaffold_overlap)),
                 source_hashes_verified=True, training_steps=0, measured_pairs_exported=0,
                 pair_policy='group by verified assay and chemical identity before pairing; conflicts quarantined; repeated controls share one cluster',
                 resources=elapsed(start))
    out = Path('runs/intervention/audit-v1')
    out.mkdir(parents=True, exist_ok=True)
    rowpath = out/'rows.jsonl'
    with rowpath.open('x') as f:
        for row in allrows: f.write(json.dumps(row,allow_nan=False)+'\n')
    Path(str(rowpath)+'.sha256').write_text(f'{sha(rowpath)}  {rowpath.name}\n')
    write_json(out/'audit.json',audit)
    print(json.dumps({k:{kk:vv for kk,vv in v.items() if kk!='conflicting_groups'}|{'conflicting_groups':len(v['conflicting_groups'])} for k,v in summaries.items()},indent=2))
    print(json.dumps([{k:v for k,v in x.items() if k!='conflicts'} for x in overlaps],indent=2))
    print(json.dumps(audit['transma_scaffold_overlap']))
    print(json.dumps(audit['resources']))


if __name__ == '__main__':
    argparse.ArgumentParser(description=__doc__).parse_args()
    main()
