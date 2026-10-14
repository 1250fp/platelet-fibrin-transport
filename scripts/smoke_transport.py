"""Small synthetic installation/regression check, not a clot validation study."""
import csv
import math
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]

def main():
    exe = ROOT/'build/bin'/('transport_reflection_long.exe' if os.name == 'nt' else 'transport_reflection_long')
    if not exe.is_file():
        raise RuntimeError('First run scripts/build_transport.py')
    with tempfile.TemporaryDirectory() as directory:
        p = Path(directory)
        scene = p/'scene.tsv'; scene.write_text('0 0 0 0 1 0 5\n')
        domain = p/'domain.txt'; domain.write_text('2 2 -5 -5 5 5\n1 1\n1 1\n')
        queries = p/'steps.tsv'
        queries.write_text('-2 0 2.5 4 0 0 0.25\n-2 3 2.5 4 0 0 0.25\n-2 0 2.5 0 0 4 0.25\n')
        result = p/'result.tsv'
        subprocess.run([str(exe),'steps',str(scene),str(domain),str(queries),str(result)],check=True)
        rows = [[float(x) for x in row.split()] for row in result.read_text().splitlines()]
        assert len(rows) == 3
        assert math.isclose(rows[0][0], -4.5, abs_tol=1e-7) and rows[0][4] >= 1
        assert math.isclose(rows[1][0], 2, abs_tol=1e-8) and rows[1][4] == 0
        assert rows[2][2] == 5 and rows[2][3] == 2
        initial = p/'initial.tsv'; initial.write_text('-2 0 2.5\n'*16)
        events = p/'events.csv'
        subprocess.run([str(exe),'run',str(scene),str(domain),str(initial),str(events),
                        '16','0.01','0.000015625','1234','0.25','0.4','0'],check=True)
        with events.open() as f: result = list(csv.DictReader(f))
        assert len(result) == 16
        assert all(0 < float(r['time_s']) <= 0.01 for r in result)
        assert all(float(r['min_clearance_um']) >= 0.25-2e-10 for r in result)
    print('PASS: solid reflection, clear motion, axial exit, and 16 synthetic Brownian histories')

if __name__ == '__main__':
    main()
