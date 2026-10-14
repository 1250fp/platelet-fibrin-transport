"""Regenerate outcome plots from included events, without rerunning simulations."""
from pathlib import Path
import importlib.util
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

def main():
    path = ROOT/'archive/research_drivers/update_nine_replica_15s_assets_20261001.py'
    spec = importlib.util.spec_from_file_location('archived_plotting', path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    folder = ROOT/'data/processed/homogeneous_pilot'
    rows = module.read_csv(folder/'state_metrics.csv')
    pairs = module.read_csv(folder/'contrastes_vs_control.csv')
    bykey = {(int(r['replicate']), r['f0_nN'], int(r['mechanical_time_s'])):r for r in rows}
    data = {sid:np.genfromtxt(folder/'states'/sid/'transport_merged.csv', delimiter=',', names=True, ndmin=1)
            for sid in sorted({r['canonical_geometry'] for r in rows})}
    if len(data) != 117 or sum(len(d) for d in data.values()) != 29952:
        raise RuntimeError('Unexpected production cohort size')
    out = ROOT/'build/reproduced'; out.mkdir(parents=True, exist_ok=True)
    module.make_assets(rows, pairs, data, bykey, out)
    print('Reproduced existing-data figures in', out)

if __name__ == '__main__':
    main()
