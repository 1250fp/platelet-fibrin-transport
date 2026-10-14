"""Compile both manuscripts in isolation, then refresh their versioned PDFs."""
import argparse
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--latexmk', default='latexmk')
    args = parser.parse_args()
    exe = shutil.which(args.latexmk)
    if not exe:
        parser.error('Install TeX with latexmk, or supply --latexmk PATH')
    paper = ROOT/'paper'
    build = paper/'build'; build.mkdir(exist_ok=True)
    for name in ['main2vPT', 'validation_supplement']:
        subprocess.run([exe, '-pdf', '-interaction=nonstopmode', '-halt-on-error',
                        '-outdir=build', name+'.tex'], cwd=paper, check=True)
        log = (build/(name+'.log')).read_text(errors='replace')
        if any(x in log for x in ['There were undefined references', 'Citation `', 'LaTeX Error']):
            raise RuntimeError(f'Unresolved citations or build errors in {name}')
        shutil.copy2(build/(name+'.pdf'), paper/(name+'.pdf'))
        print('Refreshed:', paper/(name+'.pdf'))

if __name__ == '__main__':
    main()
