"""End-of-experiment analysis only; no post-hoc margin or sample-size changes."""
from pathlib import Path
import json
import numpy as np
from scipy.stats import beta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'outputs/temporal_confirmation_20261009'


def verdict(low, high, margin):
    if low > -margin and high < margin:
        return 'within_margin'
    if high < -margin or low > margin:
        return 'outside_margin'
    return 'inconclusive'


def survival_interval(a, b):
    difference = a.astype(int)-b.astype(int)
    n = len(difference)
    plus, minus = int((difference == 1).sum()), int((difference == -1).sum())
    tail = .025/4
    def bounds(k):
        return (0. if k == 0 else float(beta.ppf(tail, k, n-k+1)),
                1. if k == n else float(beta.ppf(1-tail, k+1, n-k)))
    p, m = bounds(plus), bounds(minus)
    low, high = p[0]-m[1], p[1]-m[0]
    return dict(difference=float(difference.mean()), lower=low, upper=high,
                discordant_plus=plus, discordant_minus=minus,
                verdict=verdict(low, high, .02))


def rmst_interval(a, b, seed):
    rng = np.random.default_rng(seed)
    draws = []
    for offset in range(0, 20000, 100):
        idx = rng.integers(0, len(a), (min(100, 20000-offset), len(a)))
        af, bf = a[idx].mean(axis=1), b[idx].mean(axis=1)
        draws.append((af-bf)/bf)
    low, high = map(float, np.quantile(np.concatenate(draws), [.0125, .9875]))
    return dict(difference=float((a.mean()-b.mean())/b.mean()), lower=low, upper=high,
                verdict=verdict(low, high, .03), bootstrap_seed=seed, bootstrap_resamples=20000)


def survival(data, t):
    exits = np.sort(data['time_s'][data['outcome'] != 0])
    return 1-np.searchsorted(exits, t, side='right')/len(data)


def main():
    from temporal_confirmation_20261009 import load, save, sha, validate_events, verify_hashes
    execution = load(OUT/'execution_summary.json')
    p = load(OUT/'protocol.json')
    assert execution['completed'] and len(execution['jobs']) == len(p['jobs'])
    assert len({j['job']['name'] for j in execution['jobs']}) == len(p['jobs'])
    verify_hashes(p['input_sha256'])
    records = {}
    for result in execution['jobs']:
        j = result['job']
        assert result['exit_code'] == 0 and result['sha256'] == sha(j['output'])
        data = validate_events(j['output'], j)
        records.setdefault((j['state_id'], j['method']), []).append(data)
    combined, metrics, comparisons = {}, [], []
    for si, state in enumerate(p['states']):
        sid = state['state_id']
        reference_xyz = np.loadtxt(state['initial'])
        for method in p['steps_s']:
            data = np.concatenate(records[sid, method])
            data = data[np.argsort(data['particle'])]
            assert np.array_equal(data['particle'], np.arange(p['N_per_geometry']))
            assert np.array_equal(np.c_[data['x0'], data['y0'], data['z0']], reference_xyz)
            combined[sid, method] = data
            np.savetxt(OUT/f'{sid}_{method}.csv', np.column_stack([data[n] for n in data.dtype.names]),
                       delimiter=',', header=','.join(data.dtype.names), comments='', fmt='%.17g')
            metrics.append(dict(state_id=sid, method=method, N=len(data),
                 escaped=int((data['outcome'] != 0).sum()), censored=int((data['outcome'] == 0).sum()),
                 S100=float((data['outcome'] == 0).mean()), RMST100_s=float(data['time_s'].mean()),
                 minimum_checked_clearance_um=float(data['min_clearance_um'].min())))
        for ci, (coarse, fine) in enumerate([('dt15', 'dt7'), ('dt31', 'dt15')]):
            a, b = combined[sid, coarse], combined[sid, fine]
            S = survival_interval(a['outcome'] == 0, b['outcome'] == 0)
            R = rmst_interval(a['time_s'], b['time_s'], 2026100917+10*si+ci)
            comparisons.append(dict(state_id=sid, coarse=coarse, fine=fine,
                 role='primary' if ci == 0 else 'secondary', S=S, RMST=R,
                 joint_within_margin=S['verdict'] == R['verdict'] == 'within_margin'))
    primary = [q for q in comparisons if q['role'] == 'primary']
    overall = all(q['joint_within_margin'] for q in primary)
    results = dict(completed=True, independent_geometries=len(p['states']),
         N_per_geometry=p['N_per_geometry'], metrics=metrics, comparisons=comparisons,
         all_primary_geometries_within_margin=overall, numerical_failures=0,
         scope='Three selected geometries; frozen transport; specified release, radius, diffusion, boundaries and100s. Not all117states; not mechanical or stereological validation.',
         uniform_curve_equivalence_tested=False, manuscript_modified=False,
         source_protocol_sha256=sha(OUT/'protocol.json'))
    save(OUT/'results.json', results)
    plt.rcParams.update({'font.family':'DejaVu Serif', 'mathtext.fontset':'stix', 'font.size':10})
    t = np.linspace(0, 100, 1001)
    fig, axes = plt.subplots(3, 2, figsize=(10, 10), layout='constrained')
    colors = {'dt31':'#7b8c99', 'dt15':'#277a9c', 'dt7':'#c35e35'}
    labels = {'dt31':'0.03125 ms', 'dt15':'0.015625 ms', 'dt7':'0.0078125 ms'}
    for row, state in enumerate(p['states']):
        sid = state['state_id']
        for method in p['steps_s']:
            axes[row,0].step(t, survival(combined[sid,method], t), where='post',
                              color=colors[method], label=labels[method])
        axes[row,0].set(title=f"R{state['replicate']}, $f_0$={state['f0_nN']:g} nN, $t_m$={state['mechanical_time_s']:g} s",
                         ylabel='$S(t)$', ylim=(0, 1.02))
        axes[row,0].legend(fontsize=8)
        difference = 100*(survival(combined[sid,'dt15'],t)-survival(combined[sid,'dt7'],t))
        axes[row,1].axhspan(-2, 2, color='#e6eedf')
        axes[row,1].axhline(0, color='gray', lw=.7)
        axes[row,1].step(t, difference, where='post', color='#277a9c')
        axes[row,1].set(ylabel='Diferencia de supervivencia (pp)',
                         title='0.015625 − 0.0078125 ms · descriptivo')
        for ax in axes[row]:
            ax.set_xlabel('Tiempo de transporte (s)')
    fig.savefig(OUT/'supervivencia.png', dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(1,2,figsize=(10,4.2), layout='constrained')
    for ax, key, margin, label in [(axes[0],'S',2,'Diferencia S(100 s) (pp)'),
                                 (axes[1],'RMST',3,'Diferencia relativa RMST(100 s) (%)')]:
        ax.axvspan(-margin,margin,color='#e6eedf')
        ax.axvline(0,color='gray',lw=.7)
        for i, q in enumerate(primary):
            v=q[key]
            color='#277a9c' if v['verdict']=='within_margin' else '#b46624'
            ax.plot([100*v['lower'],100*v['upper']],[i,i],color=color,lw=2)
            ax.scatter(100*v['difference'],i,color=color)
        ax.set(yticks=range(3),yticklabels=[s['state_id'] for s in p['states']],xlabel=label)
        ax.invert_yaxis()
    fig.suptitle('Contraste primario · IC marginales97.5% · 8192 pares por geometría')
    fig.savefig(OUT/'equivalencia.png',dpi=180)
    plt.close(fig)
    translation={'within_margin':'Cumple','outside_margin':'Fuera del margen','inconclusive':'Inconcluso'}
    lines=['# Prueba temporal confirmatoria del transporte', '',
           '## Conclusión', '',
           ('**Los dos criterios primarios cumplen en las tres geometrías comprobadas.** El paso0.015625ms queda respaldado, dentro de los márgenes y alcance prefijados, frente a0.0078125ms.' if overall else
            '**No se demuestra equivalencia conjunta en todas las geometrías comprobadas.** No corresponde marcar la convergencia temporal completa como aprobada ni ampliar los márgenes para conseguirlo.'), '',
           'Se completaron73728 trayectorias:8192 posiciones por geometría,3 geometrías y3 pasos. No se repitió la mecánica, no se modificó la producción y no hubo análisis de resultados intermedios.', '',
           '## Resultados', '',
           '| Estado | Paso (ms) | Escapan | Censurados | S(100 s) | RMST(100 s) |',
           '|---|---:|---:|---:|---:|---:|']
    for m in metrics:
        lines.append(f"| {m['state_id']} | {1000*p['steps_s'][m['method']]:g} | {m['escaped']} | {m['censored']} | {100*m['S100']:.2f}% | {m['RMST100_s']:.3f} s |")
    lines += ['', '## Contraste primario:0.015625 menos0.0078125ms', '',
              '| Estado | ΔS e IC (pp) | Decisión S | ΔRMST relativo e IC (%) | Decisión RMST |',
              '|---|---:|---|---:|---|']
    for q in primary:
        s,r=q['S'],q['RMST']
        lines.append(f"| {q['state_id']} | {100*s['difference']:.3f} [{100*s['lower']:.3f}, {100*s['upper']:.3f}] | {translation[s['verdict']]} | {100*r['difference']:.3f} [{100*r['lower']:.3f}, {100*r['upper']:.3f}] | {translation[r['verdict']]} |")
    lines += ['', '![Criterios primarios](equivalencia.png)', '',
              'Los márgenes permanecen ±2pp y ±3%. Se conservó el procedimiento conservador del piloto anterior: cotas exactas para discordancias y bootstrap pareado de20000 remuestreos. Los IC son condicionales a cada geometría; no miden la variación entre redes. La cobertura bootstrap es aproximada.', '',
              '## Contraste diagnóstico:0.03125 menos0.015625ms', '',
              '| Estado | ΔS e IC (pp) | ΔRMST relativo e IC (%) | Conjunto |', '|---|---:|---:|---|']
    for q in comparisons:
        if q['role'] != 'secondary':
            continue
        s,r=q['S'],q['RMST']
        lines.append(f"| {q['state_id']} | {100*s['difference']:.3f} [{100*s['lower']:.3f}, {100*s['upper']:.3f}] | {100*r['difference']:.3f} [{100*r['lower']:.3f}, {100*r['upper']:.3f}] | {'Cumple' if q['joint_within_margin'] else 'No establecido'} |")
    lines += ['', '## Supervivencia', '', '![Supervivencia y diferencias](supervivencia.png)', '',
              'Las curvas y las bandas±2pp son descriptivas. La prueba fijada corresponde a los indicadores a100s, no a equivalencia uniforme en todo el intervalo.', '',
              '## Qué no concluye esta prueba', '',
              '- No valida el número de plaquetas, la estereología ni la calibración del tiempo mecánico.',
              '- No demuestra convergencia espacial o conectividad completa.',
              '- No certifica las117 geometrías de producción ni otros radios, D, fronteras o redes móviles.',
              '- Si un intervalo cruza el margen, el resultado es inconcluso, no automáticamente incorrecto. Si queda completamente fuera, es evidencia de una diferencia relevante bajo el criterio fijado.',
              '- Si el paso actual no cumple, el paso más fino no queda por ello validado; necesitaría una referencia adicional.', '',
              '## Integridad', '',
              'Hashes de fuentes y geometrías verificados; replays del ejecutor histórico idénticos en todos los campos; acoplamiento browniano contrastado en espacio libre; errores numéricos nunca convertidos en censura. Geometrías de producción y minipaper intactos.', '',
              'Protocolo completo: [protocolo.md](protocolo.md). Resultados y contadores: `results.json`; registros por trayectoria y bloques originales conservados.', '']
    (OUT/'informe.md').write_text('\n'.join(lines),encoding='utf-8')
    print(json.dumps({'status':'analysis_complete','all_primary_geometries_within_margin':overall}),flush=True)


if __name__ == '__main__':
    main()
