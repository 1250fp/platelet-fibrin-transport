"""Build unchanged archived C++17 engines; generated binaries are not tracked."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--compiler', default=os.environ.get('CXX', 'g++'))
    args = p.parse_args()
    compiler = shutil.which(args.compiler)
    if not compiler:
        p.error('C++17 compiler not found; use --compiler with its executable path')
    out = ROOT/'build/bin'; out.mkdir(parents=True, exist_ok=True)
    names = ['transport_reflection_long', 'transport_reflection_temporal_20261009',
             'common_domain_connectivity', 'homogeneous_clearance']
    for name in names:
        target = out/(name + ('.exe' if os.name == 'nt' else ''))
        subprocess.run([compiler, '-std=c++17', '-O3', '-DNDEBUG',
                        str(ROOT/'src/continuous_transport'/f'{name}.cpp'), '-o', str(target)], check=True)
        print(target)

if __name__ == '__main__':
    main()
