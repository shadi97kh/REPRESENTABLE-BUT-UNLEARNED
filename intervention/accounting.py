"""Measured process/child CPU accounting; no historical budget is inferred."""
from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path
import resource
import time


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def snapshot():
    me = resource.getrusage(resource.RUSAGE_SELF)
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    return dict(wall=time.monotonic(), cpu=me.ru_utime + me.ru_stime,
                child_cpu=children.ru_utime + children.ru_stime,
                peak_rss_kib=me.ru_maxrss)


def elapsed(start):
    end = snapshot()
    cpu = end['cpu'] - start['cpu']
    child = end['child_cpu'] - start['child_cpu']
    return dict(wall_s=end['wall'] - start['wall'], process_cpu_s=cpu,
                waited_child_cpu_s=child, cpu_core_hours=(cpu + child) / 3600,
                gpu_hours=0, peak_rss_kib=end['peak_rss_kib'],
                thread_limits={k: os.environ.get(k) for k in
                               ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')},
                historical_remaining_allowance=None,
                scope='instrumented region; excludes earlier tool/editor work and imports before snapshot')


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as out:
        json.dump(data, out, indent=2, allow_nan=False)
        out.write('\n')
    Path(str(path) + '.sha256').write_text(f'{sha(path)}  {path.name}\n')


def write_bytes(path, content):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != content:
            raise FileExistsError(f'refusing changed content: {path}')
    else:
        with path.open('xb') as out:
            out.write(content)
    return dict(path=str(path), sha256=sha(path), bytes=path.stat().st_size)
