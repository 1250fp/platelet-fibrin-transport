"""Cylinder clot XML for Michael's bend-model. Usage:
python gen_cylinder_min.py template.xml out.xml [seed] [R] [H]"""
import re, sys
import numpy as np

tpl_path, out_path = sys.argv[1], sys.argv[2]
seed, R, H = (int(sys.argv[3]) if len(sys.argv) > 3 else 1,
              float(sys.argv[4]) if len(sys.argv) > 4 else 24.0,
              float(sys.argv[5]) if len(sys.argv) > 5 else 32.0)
Rs = 24.0                                   # template sphere radius (um)
rng = np.random.default_rng(seed)

# 1. read the template (Michael's sphere XML)
txt = open(tpl_path, encoding="utf-8").read()
grab = lambda tag, f: np.array([list(map(f, s.split())) for s in re.findall(rf"<{tag}>(.*?)</{tag}>", txt)])
tn, tp, tl = grab("node", float), grab("plt", float), grab("link", int)

def cyl(n):                                 # uniform points in the cylinder
    r, th, z = R * np.sqrt(rng.random(n)), 2 * np.pi * rng.random(n), H * (rng.random(n) - 0.5)
    return np.c_[r * np.cos(th), r * np.sin(th), z]
dist = lambda p: np.sqrt(((p[:, None] - p[None]) ** 2).sum(-1))

# 2. same densities as the template
scale = np.pi * R**2 * H / (4 / 3 * np.pi * Rs**3)
n_nodes, n_plts = round(len(tn) * scale), round(len(tp) * scale)

# 3. nodes + target links per node, from the template's mix
deg = np.bincount(tl.ravel(), minlength=len(tn))
vals, cnt = np.unique(deg[deg > 0], return_counts=True)
nodes = cyl(n_nodes)
stubs = rng.choice(vals, size=n_nodes, p=cnt / cnt.sum()).astype(float)
if stubs.sum() % 2: stubs[rng.integers(n_nodes)] += 1

# 4. length preference w(d) = template link lengths / template pair distances
bins = np.arange(0, 2 * max(Rs, np.hypot(2 * R, H)) + 1, 1.0)
hl = np.histogram(np.linalg.norm(tn[tl[:, 0]] - tn[tl[:, 1]], axis=1), bins)[0]
hd = np.histogram(dist(tn)[np.triu_indices(len(tn), 1)], bins)[0]
w = np.convolve(np.where(hd > 0, hl / np.maximum(hd, 1), 0.0), [.25, .5, .25], "same")
W = (w / w.max())[np.clip(np.digitize(dist(nodes), bins) - 1, 0, len(w) - 1)]
np.fill_diagonal(W, 0)

# 5. join nodes: busiest node picks a partner with prob ~ w(d) * open slots
links = []
while (stubs > 0).sum() >= 2:
    o = np.flatnonzero(stubs > 0)
    i = rng.choice(o[stubs[o] == stubs[o].max()])
    p = W[i] * stubs; p[i] = 0
    if p.sum() <= 0: stubs[i] = 0; continue
    j = rng.choice(len(p), p=p / p.sum())
    links.append((min(i, j), max(i, j))); stubs[i] -= 1; stubs[j] -= 1
    W[i, j] = W[j, i] = 0                   # no duplicate links

# 6. platelets: uniform, centres at least 2.26 um apart
plts = []
while len(plts) < n_plts:
    c = cyl(1)[0]
    if all(np.linalg.norm(c - q) >= 2.26 for q in plts): plts.append(c)

# 7. write: template <settings> copied verbatim
f = lambda v: f"{v:.5f}".rstrip("0").rstrip(".")
xyz = lambda tag, P: [f"        <{tag}>{f(x)} {f(y)} {f(z)}</{tag}>" for x, y, z in P]
out = ['<?xml version="1.0" encoding="UTF-8"?>', "<data>",
       re.search(r"<settings>.*?</settings>", txt, re.S).group(0),
       "    " + re.search(r"<nodes\b[^>]*>", txt).group(0), *xyz("node", nodes), "    </nodes>",
       "    " + re.search(r"<plts\b[^>]*>(?:<!--.*?-->)?", txt).group(0), *xyz("plt", plts), "    </plts>",
       "    <links>", *[f"        <link>{a} {b}</link>" for a, b in links], "    </links>",
       "    <fixed>", "    </fixed></data>", ""]
open(out_path, "w", encoding="utf-8", newline="\n").write("\n".join(out))
print(f"nodes={len(nodes)} links={len(links)} plts={len(plts)} -> {out_path}")
