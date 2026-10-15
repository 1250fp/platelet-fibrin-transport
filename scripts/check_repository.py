"""Read-only integrity and outcome checks for the distributed research snapshot."""
import csv
import hashlib
import json
import math
from pathlib import Path
import re
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]


def check_distribution_documents():
    roots = ['docs', 'validation', 'data', 'archive', 'configs', 'scripts', 'src', 'tests', 'paper']
    paths = list(ROOT.glob('*.md'))
    for directory in roots:
        paths.extend(p for p in (ROOT/directory).rglob('*')
                     if p.is_file() and p.suffix in {'.md', '.json', '.csv', '.py', '.cpp', '.tex', '.bib'}
                     and not {'build', '__pycache__', 'raw'}.intersection(p.relative_to(ROOT).parts))
    personal_root = re.compile(r'(?:[A-Za-z]:[\\/]+Users[\\/]+|[A-Za-z]:[\\/]+tesisprog[\\/]+)', re.I)
    for path in paths:
        contents = path.read_text(encoding='utf-8')
        assert not personal_root.search(contents), f'Personal filesystem root in {path}'
        if path.suffix != '.md':
            continue
        for target in re.findall(r'\]\(([^)]+)\)', contents):
            target = unquote(target.strip('<>').split('#', 1)[0])
            if not target or re.match(r'^[A-Za-z][A-Za-z0-9+.-]*:', target):
                continue
            assert (path.parent/target).exists(), f'Broken relative link in {path}: {target}'


def main():
    manifest = json.loads((ROOT/'data/file_manifest.json').read_text())
    assert manifest['schema_version'] == 2
    assert len({r['file'] for r in manifest['records']}) == len(manifest['records'])
    for record in manifest['records']:
        p = ROOT/record['file']
        assert p.is_file(), f'Missing {p}'
        actual = hashlib.sha256(p.read_bytes()).hexdigest()
        assert actual == record['sha256'], f'Snapshot changed: {p}'
        assert p.stat().st_size == record['bytes'], f'Size changed: {p}'
        assert re.fullmatch(r'[0-9a-f]{64}', record['source_sha256'])
        for alias in record.get('aliases', []):
            assert not (ROOT/alias).exists(), f'Redundant alias restored: {alias}'
    folder = ROOT/'data/processed/homogeneous_pilot'
    with (folder/'state_metrics.csv').open() as f: metrics = list(csv.DictReader(f))
    canonical = {r['canonical_geometry'] for r in metrics}
    assert len(metrics) == 144 and len(canonical) == 117
    total = 0
    for sid in canonical:
        with (folder/'states'/sid/'transport_merged.csv').open() as f: rows=list(csv.DictReader(f))
        assert len(rows)==256 and {int(r['particle']) for r in rows}==set(range(256))
        t=[float(r['time_s']) for r in rows]
        assert all(math.isfinite(v) and 0 < v <= 100 for v in t)
        censored=sum(int(r['outcome'])==0 for r in rows)
        assert all(int(r['outcome']) in (0,1,2) for r in rows)
        assert all(float(r['time_s'])==100 for r in rows if int(r['outcome'])==0)
        for m in [m for m in metrics if m['canonical_geometry']==sid]:
            assert math.isclose(float(m['S100']),censored/256,abs_tol=1e-12)
            assert math.isclose(float(m['RMST100_s']),sum(t)/256,abs_tol=1e-9)
        total+=len(rows)
    assert total==29952
    temporal=list((ROOT/'validation/temporal_20261009/events').glob('*.csv'))
    assert len(temporal)==9
    count=0
    for path in temporal:
        with path.open() as f: rows=list(csv.DictReader(f))
        assert len(rows)==8192
        count+=len(rows)
    assert count==73728
    for filename in ['main2vPT.tex','validation_supplement.tex']:
        assert r'\bibliography{refs}' in (ROOT/'paper'/filename).read_text()
    check_distribution_documents()
    print(f'PASS: {len(manifest["records"])} distributed hashes with original provenance; 117 unique production states / {total} events; {count} temporal events; shared bibliography; Markdown links and personal-path checks')

if __name__ == '__main__':
    main()
