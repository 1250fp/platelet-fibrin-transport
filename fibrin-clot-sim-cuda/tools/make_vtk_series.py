#!/usr/bin/env python3
"""Index a completed native run for ParaView. Python standard library only.

The solver saves Params and the three VTK families in the same output event.
Pair ordered snapshot times with contiguous frame indices, never guessed times.
Existing identical indexes are preserved; conflicting indexes are not replaced.
"""
import argparse
import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile

FAMILIES = {'FibrinNetwork': 'FibrinNetwork_', 'Platelet': 'Platelet', 'PlateletConn': 'PlateletConn'}


def single_dir(run, pattern):
    found = [p for p in run.glob(pattern) if p.is_dir()]
    if len(found) != 1:
        raise ValueError(f'Expected one {pattern} directory, found {len(found)}')
    return found[0]


def build_series(run, expected_final_time=None):
    run = Path(run).resolve(strict=True)
    if (run / 'exit_code.txt').read_text().strip() != '0':
        raise ValueError('Solver did not finish successfully; no indexes written')
    animation = single_dir(run, 'Animation*')
    params = single_dir(run, 'Params*')
    times = []
    for path in params.glob('Param_*.sta'):
        filename_time = float(path.stem.removeprefix('Param_'))
        with path.open() as stream:
            header = stream.readline().split()
        if len(header) != 2 or header[0] != 'time':
            raise ValueError(f'Missing time header: {path.name}')
        time = float(header[1])
        if not math.isfinite(time) or time < 0 or abs(time - filename_time) > 1e-5:
            raise ValueError(f'Invalid or inconsistent time: {path.name}')
        times.append(time)
    times.sort()
    if not times or any(a >= b for a, b in zip(times, times[1:])):
        raise ValueError('Snapshot times must be nonempty and strictly increasing')
    if expected_final_time is not None:
        if not math.isfinite(expected_final_time) or abs(times[-1] - expected_final_time) > 1e-5:
            raise ValueError(f'Last snapshot is {times[-1]} s, expected {expected_final_time} s')

    pending = []
    for family, prefix in FAMILIES.items():
        pattern = re.compile(re.escape(prefix) + r'(\d+)\.vtk')
        frames = sorted((int(m[1]), p) for p in animation.iterdir()
                        if (m := pattern.fullmatch(p.name)) and p.is_file())
        if [i for i, p in frames] != list(range(1, len(times) + 1)):
            raise ValueError(f'{family}: frame indices must be 1..{len(times)}, matching Params')
        for _, path in frames:
            with path.open('rb') as stream:
                if not stream.readline().startswith(b'# vtk DataFile Version '):
                    raise ValueError(f'Empty or invalid VTK header: {path.name}')
        data = {'file-series-version': '1.0',
                'files': [{'name': p.name, 'time': t} for (_, p), t in zip(frames, times)]}
        target = animation / (family + '.vtk.series')
        if target.exists():
            if json.loads(target.read_text(encoding='utf-8-sig')) != data:
                raise ValueError(f'Conflicting existing index; refusing to overwrite: {target.name}')
        else:
            pending.append((target, data))

    # Validate all three families before creating any index; publish each atomically.
    for target, data in pending:
        fd, temporary = tempfile.mkstemp(prefix='.' + target.name, suffix='.tmp', dir=animation)
        try:
            with os.fdopen(fd, 'w', encoding='utf-8', newline='\n') as stream:
                json.dump(data, stream, indent=2, allow_nan=False)
                stream.write('\n')
            os.replace(temporary, target)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    print(f'Validated 3 series, {len(times)} frames each, {times[0]:g} to {times[-1]:g} seconds; '
          f'created {len(pending)} indexes in {animation}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    parser.add_argument('--expected-final-time', type=float)
    args = parser.parse_args()
    try:
        build_series(args.run, args.expected_final_time)
    except (ValueError, OSError) as error:
        print(f'Series generation failed: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
