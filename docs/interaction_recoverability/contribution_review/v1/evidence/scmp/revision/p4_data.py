"""Row-level P3 reconciliation; no fitting and no edits to historical records."""
from __future__ import annotations

from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

from .p4_resources import checked_record, sha


def norm(s):
    return s.upper().replace('T', 'U')


def rc(s):
    return norm(s).translate(str.maketrans('ACGU', 'UGCA'))[::-1]


def transcripts(path, wanted):
    """Keep the actual transcript ID, sequence and coordinate system."""
    out = defaultdict(list)
    cur, buf = None, []
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if cur and cur[1] in wanted:
                    out[cur[1]].append({'id': cur[0], 'seq': norm(''.join(buf))})
                fields = line[1:].strip().split('|')
                cur, buf = (fields[0], fields[5].upper()), []
            elif cur and cur[1] in wanted:
                buf.append(line.strip())
        if cur and cur[1] in wanted:
            out[cur[1]].append({'id': cur[0], 'seq': norm(''.join(buf))})
    return dict(out)


def all_positions(s, sub):
    p = s.find(sub)
    while p >= 0:
        yield p
        p = s.find(sub, p + 1)


def components(seqs, k=13):
    """Transitive global shared-k-mer components; no gene-qualified loophole."""
    parent = list(range(len(seqs)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    seen = {}
    for i, s in enumerate(seqs):
        for token in {s[j:j+k] for j in range(max(1, len(s)-k+1))}:
            if token in seen:
                a, b = root(i), root(seen[token])
                parent[max(a, b)] = min(a, b)
            else:
                seen[token] = i
    return [root(i) for i in range(len(seqs))]


def partition(groups, seed, train_fraction, dev_fraction):
    keys = np.array(sorted(set(groups)))
    np.random.RandomState(seed).shuffle(keys)
    a, b = int(len(keys)*train_fraction), int(len(keys)*(train_fraction+dev_fraction))
    sets = [set(keys[:a]), set(keys[a:b]), set(keys[b:])]
    return [np.array([i for i, g in enumerate(groups) if g in s], dtype=int) for s in sets]


def cross_partition(rows, train, test):
    if not train or not test:
        return {'status': 'empty_partition'}
    seqs = [r['guide'] for r in rows]
    lookup = defaultdict(set)
    for i in train:
        for j in range(len(seqs[i])-12):
            lookup[seqs[i][j:j+13]].add(i)
    pairs = set()
    for i in test:
        for j in range(len(seqs[i])-12):
            for other in lookup[seqs[i][j:j+13]]:
                pairs.add((other, i))
    encoded = np.array([[ord(c) for c in s] for s in seqs], dtype=np.uint8)
    nearest, examples = [], []
    for j in test:
        dist = (encoded[train] != encoded[j]).sum(axis=1)
        i = train[int(np.argmin(dist))]
        nearest.append(int(dist.min()))
        if dist.min() <= 2 and len(examples) < 12:
            examples.append({'train_id': rows[i]['row_id'], 'test_id': rows[j]['row_id'],
                             'mismatches': int(dist.min()), 'different_gene': rows[i]['gene'] != rows[j]['gene']})
    overlaps = []
    for i in train:
        a = rows[i]
        for j in test:
            b = rows[j]
            if (a['gene'], a['legacy_transcript_id']) != (b['gene'], b['legacy_transcript_id']):
                continue
            n = min(a['legacy_pos']+a['legacy_match_nt'], b['legacy_pos']+b['legacy_match_nt']) - max(a['legacy_pos'], b['legacy_pos'])
            if n > 0:
                overlaps.append((i, j, n))
    return {'n_train': len(train), 'n_test': len(test),
            'exact_duplicate_pairs': sum(seqs[i] == seqs[j] for i, j in pairs),
            'shared_any_13mer_pairs': len(pairs),
            'shared_13mer_test_rows': len({j for _, j in pairs}),
            'cross_gene_shared_13mer_pairs': sum(rows[i]['gene'] != rows[j]['gene'] for i, j in pairs),
            'nearest_hamming_mismatch_histogram': dict(sorted(Counter(nearest).items())),
            'maximum_aligned_identity': 1-min(nearest)/21,
            'same_transcript_overlapping_site_pairs': len(overlaps),
            'overlap_test_rows': len({j for _, j, _ in overlaps}),
            'near_duplicate_examples': examples,
            'replicate_leakage_status': 'raw assay replicate IDs unavailable; duplicate sequences and overlapping sites audited separately'}


def audit_rows(cfg):
    from experiments._huesken import load_verified
    recs, _, _ = load_verified()
    ann = list(csv.DictReader(open('data/Huesken_2431_annotated.tsv'), delimiter='\t'))
    by_guide = defaultdict(list)
    for i, a in enumerate(ann):
        by_guide[norm(a['Antisense_21mer'])].append((i, a))
    tx = transcripts('data/gencode/gencode.v47.pc_transcripts.fa.gz',
                     {a['Gene'].upper() for a in ann})
    rows, legacy_hits = [], []
    for idx, r in enumerate(recs):
        guide = norm(r['sequence'])
        candidates = by_guide[guide]
        # Match the old last-assignment gene map for RECONSTRUCTION ONLY.
        a = candidates[-1][1]
        gene = a['Gene']
        row = {'row_id': f'dsir:{idx}', 'dsir_index': idx, 'guide': guide,
               'gene': gene, 'accession': a['Accession_number'], 'efficacy': r['efficacy'],
               'annotation_rows': [i for i, _ in candidates],
               'annotation_gene_choices': sorted({v['Gene'] for _, v in candidates}),
               'site19': norm(a['Sense_19mer']), 'replicate_id': None,
               'published_partition': 'train' if idx < 2182 else 'test',
               'context_status': 'unresolved', 'verified_construct': False,
               'legacy_status': 'gene_absent_in_exact_case_sensitive_lookup', 'mappings': []}
        row['sense19_matches_guide'] = row['site19'] == rc(guide[:19])
        old = []
        for ti, t in enumerate(tx.get(gene, [])):
            p, mlen = t['seq'].find(rc(guide)), 21
            if p < 0:
                p, mlen = t['seq'].find(rc(guide)[2:]), 19
            if p >= 0:
                old.append((ti, p, mlen, t))
            # The supported matched site is the annotated 19-mer. Coordinates
            # are NEVER padded to 21 bases when only 19 were matched.
            for start in all_positions(t['seq'], row['site19']):
                center = start + len(row['site19'])//2
                lo, hi = center-25, center+25
                row['mappings'].append({'transcript_id': t['id'], 'iso_index': ti,
                                        'start': start, 'end': start+19,
                                        'window_start': lo, 'window_end': hi,
                                        'window': t['seq'][lo:hi] if lo >= 0 else '',
                                        'transcript_len': len(t['seq'])})
        if gene in tx:
            row['legacy_status'] = 'target_not_found'
        if old:
            ti, p, mlen, t = old[0]
            row.update(legacy_status='located', legacy_pos=p, legacy_match_nt=mlen,
                       legacy_transcript_id=t['id'], legacy_isoform_hits=len(old),
                       legacy_window=t['seq'][p-15:p+35], legacy_transcript_len=len(t['seq']))
            row['legacy_local_window_count'] = len({v['seq'][q-15:q+35] for _, q, _, v in old})
            legacy_hits.append(row)
        rows.append(row)
    hull = {}
    for r in legacy_hits:
        lo, hi = hull.setdefault(r['gene'], [10**10, -1])
        hull[r['gene']] = [min(lo, r['legacy_pos']), max(hi, r['legacy_pos']+21)]
    legacy_pilot, short = [], []
    for r in legacy_hits:
        lo, hi = hull[r['gene']]
        p = r['legacy_pos']
        inside = p-15 >= lo and p+35 <= hi
        if inside and len(r['legacy_window']) == 50:
            r['legacy_status'] = 'pilot_included'
            r['cluster'] = min(r['guide'][i:i+13] for i in range(9))
            legacy_pilot.append(r)
        elif inside:
            r['legacy_status'] = 'audit_interior_but_truncated_transcript_window'
            short.append({k: r[k] for k in ['row_id', 'gene', 'guide', 'legacy_transcript_id',
                                          'legacy_pos', 'legacy_transcript_len']})
            short[-1]['actual_window_nt'] = len(r['legacy_window'])
        else:
            r['legacy_status'] = 'legacy_gene_hull_edge'
    # Hulls can only be formed in a single transcript's coordinates and within
    # one accession. The shared reporter insert remains an assumption.
    coord = defaultdict(list)
    for r in rows:
        for m in r['mappings']:
            coord[(r['gene'], r['accession'], m['transcript_id'])].append((r, m))
    coord_summary, best = [], {}
    for key, hits in coord.items():
        starts, ends = [m['start'] for _, m in hits], [m['end'] for _, m in hits]
        count = len({r['guide'] for r, _ in hits})
        coord_summary.append({'gene': key[0], 'accession': key[1], 'transcript_id': key[2],
                              'start': min(starts), 'end': max(ends), 'distinct_guides': count})
        score = (-count, key[2])
        if key[:2] not in best or score < best[key[:2]]:
            best[key[:2]] = score
        for _, m in hits:
            m['conditional_interior'] = (len(m['window']) == 50 and
                                        m['window_start'] >= min(starts) and m['window_end'] <= max(ends))
    for r in rows:
        good = [m for m in r['mappings'] if m['conditional_interior']]
        unique = sorted({m['window'] for m in good})
        r['distinct_local_windows'] = len({m['window'] for m in r['mappings'] if len(m['window']) == 50})
        r['conditional_distinct_local_windows'] = len(unique)
        r['isoform_class'] = ('no_usable_local_window' if not unique else
                             'identical_local_windows' if len(unique) == 1 else 'different_local_windows')
        selected_id = best.get((r['gene'], r['accession']), (None, None))[1]
        selected = [m for m in good if m['transcript_id'] == selected_id]
        r['selection_transcript_id'] = selected_id
        if selected and len(r['annotation_gene_choices']) == 1 and r['sense19_matches_guide']:
            m = min(selected, key=lambda m: m['start'])
            r.update(context_status='inferred_contiguous_native_insert', window=m['window'],
                     window_start=m['window_start'], window_end=m['window_end'],
                     alternative_window=unique[-1],
                     assumption='one_shared_contiguous_native_insert_per_annotation_accession_unverified')
    stats = []
    for gene in sorted({r['gene'] for r in rows}):
        rr = [r for r in rows if r['gene'] == gene]
        ll = [r for r in rr if 'legacy_pos' in r]
        stats.append({'gene': gene, 'rows': len(rr), 'unique_guides': len({r['guide'] for r in rr}),
                      'annotation_accessions': sorted({r['accession'] for r in rr}),
                      'legacy_located_rows': len(ll),
                      'legacy_pilot_rows': sum(r['legacy_status'] == 'pilot_included' for r in rr),
                      'legacy_selected_transcript_ids': sorted({r['legacy_transcript_id'] for r in ll}),
                      'mixed_transcript_coordinates_in_legacy_hull': len({r['legacy_transcript_id'] for r in ll}) > 1,
                      'corrected_conditional_rows': sum(r['context_status'] != 'unresolved' for r in rr)})
    from .rna_pilot import novel_sequence_split
    leaks = []
    for seed in range(5):
        tr, te, split = novel_sequence_split(legacy_pilot, .25, seed)
        leaks.append({'seed': seed, 'legacy_group_check': split,
                      **cross_partition(legacy_pilot, tr, te)})
    subset = [r for r in legacy_pilot if r['context_status'] != 'unresolved']
    seq_counts = Counter(r['guide'] for r in rows)
    return {'rows': rows, 'summary': {
                'full_rows': len(rows), 'unique_guides': len(seq_counts),
                'duplicate_guide_rows_excess': sum(v-1 for v in seq_counts.values()),
                'targets': len(stats), 'legacy_located_rows': len(legacy_hits),
                'legacy_located_targets': len({r['gene'] for r in legacy_hits}),
                'legacy_audit_interior': len(legacy_pilot)+len(short), 'legacy_pilot_rows': len(legacy_pilot),
                'legacy_pilot_targets': len({r['gene'] for r in legacy_pilot}),
                'legacy_pilot_unique_guides': len({r['guide'] for r in legacy_pilot}),
                'legacy_pilot_isoform_ambiguous_rows': sum(r['legacy_isoform_hits'] > 1 for r in legacy_pilot),
                'legacy_ambiguous_identical_local_windows': sum(r['legacy_isoform_hits'] > 1 and r['legacy_local_window_count'] == 1 for r in legacy_pilot),
                'legacy_ambiguous_different_local_windows': sum(r['legacy_isoform_hits'] > 1 and r['legacy_local_window_count'] > 1 for r in legacy_pilot),
                'corrected_conditional_intersection_with_p3': len(subset),
                'corrected_context_classes_in_p3': dict(Counter(r['isoform_class'] for r in legacy_pilot)),
                'legacy_exclusions': dict(Counter(r['legacy_status'] for r in rows)),
                'assay_replicate_count': None, 'replicate_reason': 'distributed rows have no raw replicate IDs',
                'label_above_one': sum(r['efficacy'] > 1 for r in rows),
                'label_max': max(r['efficacy'] for r in rows),
                'verified_construct_rows': 0,
                'site_definition': '19-mer sequence is identifiable; physical construct locus and raw assay replicates are not'},
            'per_gene': stats, 'three_short_exclusions': short,
            'transcript_accession_hulls': coord_summary, 'legacy_cross_partition': leaks,
            'published_global_duplicates': len(set(r['guide'] for r in rows[:2182]) & set(r['guide'] for r in rows[2182:])),
            'policy': cfg['rna'], 'raw_labels_clipped': False}


def historical_evidence(cfg):
    records = {k: checked_record(p) for k, p in cfg['historical_records'].items()}
    hashes = {p: sha(p) for p in cfg['historical_records'].values()}
    p3 = records['p3_pilot']
    ref = p3['summary']['matched_sequence_chemistry']['spearman_per_seed']
    contrasts = {arm: {'paired_deltas': [a-b for a, b in zip(v['spearman_per_seed'], ref)],
                       'mean_delta': float(np.mean(np.array(v['spearman_per_seed'])-ref))}
                 for arm, v in p3['summary'].items()}
    aliases = [
        {'names': ['matched_sequence_chemistry', 'target_sequence'], 'actual_class': 'numpy ridge solve',
         'feature_dim': 85, 'fitted_coefficients': 85, 'distinct_function_fits_per_split': 1,
         'estimator': 'intercept + guide onehot; target arm permutes those features'},
        {'names': ['mfe_gnn', 'anchor_A'], 'actual_class': 'numpy ridge solve',
         'feature_dim': 89, 'fitted_coefficients': 89, 'distinct_function_fits_per_split': 1,
         'estimator': 'guide onehot + four exact-prior/global window summaries; no MFE computation'},
        {'names': ['sampled_ensemble_gnn', 'repaired_A_plus_B'], 'actual_class': 'numpy ridge solve',
         'feature_dim': 90, 'fitted_coefficients': 90, 'distinct_function_fits_per_split': 1,
         'estimator': 'previous four summaries + mean edge count over 64 valid draws; not a learned residual; cached sample shared'}]
    for a in aliases:
        a['prediction_hashes'] = None
        a['prediction_hash_status'] = 'unavailable: historical record stores metrics, not predictions or coefficients; no unbudgeted refit'
        a['score_vector_sha256'] = {arm: hashlib.sha256(np.asarray(p3['summary'][arm]['spearman_per_seed'], dtype='<f8').tobytes()).hexdigest() for arm in a['names']}
    p2 = records['p2']
    by_rung = []
    timing = defaultdict(list)
    for rung in range(4):
        cs = [c for c in p2['cases'] if c['case']['rung'] == rung]
        by_rung.append({'rung': rung, 'n': len(cs), 'coverage': {
            arm: sum(c['arms'][arm]['decision_certified'] for c in cs)/len(cs)
            for arm in cs[0]['arms']}})
    for case, raw in zip(p2['cases'], p2['raw']):
        for arm, r in raw.items():
            if 'profile' in r:
                timing[arm].append((r['profile']['wall_s'], case['arms'][arm]['wall_s']))
    fixed_timing = {arm: {'median_immediate_own_wall_s': float(np.median([x for x, _ in v])),
                          'median_reported_later_wall_s': float(np.median([y for _, y in v]))}
                    for arm, v in timing.items()}
    return {'hashes': hashes, 'p1_stages': records['p1']['stages'],
            'p2_by_rung': by_rung, 'p2_policy_gate': p2['policy_gate'],
            'p2_timing_reconciliation': fixed_timing,
            'p2_predictor': 'random CertifiedModel plus random parameter perturbation; no checkpoint loaded',
            'p2_lp_invalid_for_dag': 'degree<=1 excludes internal path vertices of degree 2; do not reuse as a sound DAG relaxation',
            'p2_timing_defect': 'Meter.summary reads a live wall clock again during score_case after all arms; charges later arms to earlier arms',
            'p2_no_learned_policy_test': True,
            'p3_aliases': aliases, 'p3_full_precision_contrasts': contrasts,
            'p3_seed_contract': 'split changes for each seed; ridge has no initialization; first ensemble arm seed 0 fills per-row sample cache reused by later seeds and aliases',
            'p3_mc_contract': 'SE of mean |S| in base-pair counts, not efficacy predictions; 64 draws in fitted arm despite config 256; independent backtraces with one pseudorandom stream, not MCMC',
            'p3_prior_contract': 'pilot theta=0; sensitivity theta=0.5*reported_temperature on first 24 windows. This is inverse-temperature/edge reward, not temperature. No old prediction sensitivity was measured.',
            'p3_eligible_increase': 4, 'p3_distinct_fitted_functions': 3,
            'p3_claim_limit': 'summary ridge comparison only; five inspected partitions are not independent biological datasets'}
