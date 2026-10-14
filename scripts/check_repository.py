"""Read-only integrity and outcome checks for the distributed research snapshot."""
import csv
import hashlib
import json
import math
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]

def main():
    manifest = json.loads((ROOT/'data/file_manifest.json').read_text())
    for record in manifest['records']:
        p = ROOT/record['file']
        assert p.is_file(), f'Missing {p}'
        actual = hashlib.sha256(p.read_bytes()).hexdigest()
        assert actual == record['sha256'], f'Snapshot changed: {p}'
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
    print(f'PASS: {len(manifest["records"])} source hashes; 117 unique production states / {total} events; {count} temporal events; shared bibliography')

if __name__ == '__main__':
    main()
