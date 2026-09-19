"""Administrative preservation, Markdown-link, manifest and ZIP checks only."""
from pathlib import Path, PurePosixPath
import hashlib
import json
import re
import resource
import subprocess
import time
from urllib.parse import unquote
import zipfile

ROOT = Path('/home/shadi/iclr2027')
OUT = ROOT / 'docs/interaction_recoverability/contribution_review/v1'
start_wall, start_cpu = time.perf_counter(), time.process_time()

def digest(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for part in iter(lambda: f.read(1024 * 1024), b''):
            h.update(part)
    return h.hexdigest()

def git(*args):
    return subprocess.check_output(['git', '--no-optional-locks', *args], cwd=ROOT)

def clean_status(data):
    return [line for line in data.decode().splitlines()
            if not line[3:].startswith('docs/interaction_recoverability/contribution_review/')]

baseline = json.loads((OUT / 'administrative/preservation_before.json').read_text())
changed, missing = [], []
for rel, expected in baseline['files'].items():
    p = ROOT / rel
    if not p.is_file():
        missing.append(rel)
    elif p.stat().st_size != expected['bytes'] or digest(p) != expected['sha256']:
        changed.append(rel)
status = git('status', '--short')
diff = git('diff', '--binary')
cached = git('diff', '--cached', '--binary')
status_same = clean_status(status) == clean_status((OUT / 'administrative/git_status_before.txt').read_bytes())
diff_same = diff == (OUT / 'administrative/git_diff_before.patch').read_bytes()
cached_same = cached == (OUT / 'administrative/git_diff_cached_before.patch').read_bytes()
new_paths = [p.decode() for p in git('ls-files', '-c', '-o', '--exclude-standard', '-z').split(b'\0') if p]
new_outside = [p for p in new_paths if p not in baseline['files']
               and not (ROOT / p).is_relative_to(OUT)]
assert not changed and not missing and not new_outside and status_same and diff_same and cached_same
preservation = {'scope': baseline['scope'], 'baseline_files_checked': len(baseline['files']),
                'changed': changed, 'missing': missing, 'new_git_visible_files_outside_review': new_outside,
                'git_status_outside_review_equal': status_same, 'tracked_diff_equal': diff_same,
                'cached_diff_equal': cached_same,
                'process_cpu_seconds_through_preservation_check': time.process_time() - start_cpu,
                'wall_seconds_through_preservation_check': time.perf_counter() - start_wall}
(OUT / 'administrative/preservation_after.json').write_text(json.dumps(preservation, indent=2) + '\n')
(OUT / 'administrative/git_status_after.txt').write_bytes(status)
(OUT / 'administrative/git_diff_after.patch').write_bytes(diff)
(OUT / 'administrative/git_diff_cached_after.patch').write_bytes(cached)

inv = json.loads((OUT / 'administrative/evidence_inventory.json').read_text())
for entry in inv['evidence']:
    assert digest(ROOT / entry['original_path']) == entry['sha256']
    assert digest(OUT / entry['packet_path']) == entry['sha256']

excluded = {'review_packet.zip', 'review_packet.zip.sha256', 'archive_checks.json', 'manifest.json', 'integrity_checks.json'}
payload = sorted(p for p in OUT.rglob('*') if p.is_file() and str(p.relative_to(OUT)) not in excluded)
names = {str(p.relative_to(OUT)) for p in payload} | {'manifest.json', 'integrity_checks.json'}
checked, external, fragments, broken = [], set(), [], []

def destinations(p):
    text = p.read_text()
    text = re.sub(r'^```[^\n]*\n.*?^```\s*$', '', text, flags=re.M | re.S)
    targets = re.findall(r'!?\[[^\]\n]+\]\(([^\)\n]+)\)', text)
    targets += re.findall(r'^\s*\[[^\]]+\]:\s*(\S+)', text, re.M)
    targets += re.findall(r'<a\s[^>]*href=[\"\']([^\"\']+)', text)
    return targets

for p in payload:
    if p.suffix != '.md':
        continue
    for raw in destinations(p):
        target = raw.split(' "')[0].strip('<>')
        if re.match(r'[a-zA-Z][\w+.-]*:', target):
            external.add(target)
            continue
        filepart, _, fragment = unquote(target).partition('#')
        dest = (p.parent / filepart).resolve() if filepart else p
        info = {'from': str(p.relative_to(OUT)), 'target': target}
        if not dest.is_relative_to(OUT) or str(dest.relative_to(OUT)) not in names:
            broken.append(info)
        else:
            checked.append(info)
        if fragment:
            # Preserved review evidence uses GitHub-style source line locators.
            match = re.fullmatch(r'L([0-9]+)(?:-L([0-9]+))?', fragment)
            assert match is not None, ('Unsupported local fragment', info)
            line_count = len(dest.read_text().splitlines())
            first = int(match.group(1))
            last = int(match.group(2) or match.group(1))
            assert 1 <= first <= last <= line_count, info
            fragments.append({**info, 'valid_source_line_locator': True, 'file_lines': line_count})
assert not broken, broken
records = []
for p in payload:
    records.append({'path': str(p.relative_to(OUT)), 'bytes': p.stat().st_size, 'sha256': digest(p)})
manifest = {'format': 'sha256 manifest for static review payload',
            'evidence_paths_are_original_repository_relative_under': 'evidence/',
            'files': records,
            'excluded_from_manifest': sorted(excluded),
            'closure': 'All active local Markdown file links remain inside the packet. Code imports are not an executable dependency closure.'}
(OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
integrity = {'status': 'PASS', 'evidence_files_identical': len(inv['evidence']),
             'historical_commitments_matching': len(inv['historical_commitments']),
             'historical_preservation': preservation, 'manifest_payload_files': len(records),
             'markdown_files_checked': sum(p.suffix == '.md' for p in payload),
             'active_local_destinations_checked': len(checked), 'broken_local_links': broken,
             'local_fragments': fragments, 'external_urls_not_mass_probed': sorted(external),
             'archive_planned_entries': len(names),
             'archive_scope': 'Payload plus manifest and this report. The outside archive_checks.json checks the final ZIP bytes.',
             'administrative_only': True, 'project_modules_imported': False,
             'earlier_administrative_attempts': 'administrative/operation_notes.json',
             'process_cpu_seconds_before_integrity_report': time.process_time() - start_cpu,
             'wall_seconds_before_integrity_report': time.perf_counter() - start_wall}
(OUT / 'integrity_checks.json').write_text(json.dumps(integrity, indent=2) + '\n')
zip_path = OUT / 'review_packet.zip'
with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for name in sorted(names):
        z.write(OUT / name, name)
with zipfile.ZipFile(zip_path) as z:
    members = z.namelist()
    assert len(members) == len(set(members))
    assert set(members) == names
    for name in members:
        q = PurePosixPath(name)
        assert not q.is_absolute() and '..' not in q.parts and '\\' not in name
        assert hashlib.sha256(z.read(name)).hexdigest() == digest(OUT / name), name
    assert z.testzip() is None
    archived_manifest = json.loads(z.read('manifest.json'))
    for entry in archived_manifest['files']:
        data = z.read(entry['path'])
        assert len(data) == entry['bytes'] and hashlib.sha256(data).hexdigest() == entry['sha256']
sha = digest(zip_path)
(OUT / 'review_packet.zip.sha256').write_text(f'{sha}  review_packet.zip\n')
# Recheck historical bytes after packaging; only new-version files were written.
assert all((ROOT / rel).is_file() and digest(ROOT / rel) == val['sha256']
           for rel, val in baseline['files'].items())
report = {'status': 'PASS', 'zip_file': 'review_packet.zip', 'zip_sha256': sha,
          'zip_bytes': zip_path.stat().st_size, 'archive_entries': len(names),
          'every_entry_matches_local_bytes': True, 'all_manifest_hashes_match': True,
          'unique_relative_paths_no_traversal': True, 'zip_crc_pass': True,
          'archive_local_links_resolve': True, 'active_local_destinations_checked': len(checked),
          'historical_files_rechecked_after_packaging': len(baseline['files']),
          'historical_changes': [], 'administrative_only': True,
          'helper_invocation': 'finalize_packet.py',
          'cost_scope': 'This helper process and waited Git children only; excludes prior review, authoring, retrieval and prepare_packet.py.',
          'process_cpu_seconds_before_report_write': time.process_time() - start_cpu,
          'waited_child_cpu_seconds': resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime + resource.getrusage(resource.RUSAGE_CHILDREN).ru_stime,
          'wall_seconds_before_report_write': time.perf_counter() - start_wall,
          'parent_max_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(OUT / 'archive_checks.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
