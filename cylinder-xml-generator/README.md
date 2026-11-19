# cylinder-xml-generator

Makes cylinder-shaped clot input XMLs for the 3D fibrin/platelet CUDA model (`fibrin-clot-sim-cuda`, a modified version of the model by C. Michael et al., *Commun. Biol.* 6, 869 (2023), https://doi.org/10.1038/s42003-023-05240-z). The new network is statistically matched to a template XML from the original model.

## Usage

```bash
pip install numpy
python gen_cylinder_min.py template_sample1.xml out.xml [seed] [R] [H]
```

The defaults are seed 1, radius R = 24 µm and height H = 32 µm, with the cylinder centred at the origin along z. Examples:

```bash
python gen_cylinder_min.py template_sample1.xml cylinder_R24_H32_seed1.xml 1
python gen_cylinder_min.py template_sample1.xml cylinder_R24_H32_seed2.xml 2
```

The same seed always gives the same file. Seeds 1–3 reproduce the example inputs in `fibrin-clot-sim-cuda/examples/` byte for byte. The output XML can be passed straight to `bend-model`.

## Method

The template is a sphere of radius 24 µm. The steps below match the code comments:

1. **Read the template**: the node, platelet and link lists.
2. **Same densities**: the node and platelet counts are scaled by the volume ratio, cylinder / template sphere.
3. **Nodes**: placed uniformly in the cylinder. Each node gets a target number of links (its degree), drawn from the template's degree distribution.
4. **Link-length preference**: w(d) = histogram of template link lengths ÷ histogram of all template node-pair distances (1 µm bins), smoothed with [0.25, 0.5, 0.25]. This makes link lengths follow the template's.
5. **Linking**: greedy stub matching. The node with the most open slots picks a partner with probability ∝ w(distance) × the partner's open slots. Self-links and duplicate links are not allowed.
6. **Platelets**: placed uniformly in the cylinder, centres at least 2.26 µm apart.
7. **Write**: the template's `<settings>` block (all model parameters) is copied verbatim, followed by the new nodes, platelets and links.

To change model parameters, edit the `<settings>` in the template, or use another template from the original code.

## Files

- `gen_cylinder_min.py`: the generator (Python 3 + NumPy).
- `template_sample1.xml`: template, an unmodified copy of `simulations/StandardSimulations/Code_Input/data_Cir_0.003_0.0008_23.811_23.811_23.811_0.5_set5_sample1strain_switch0.1tnd13_F4.xml` from the original model code (Code Ocean capsule 2690247), Copyright (c) 2022 Zhiliang Xu, MIT License. It uses the "Standard" parameter set.

## License

MIT, see `LICENSE`. The generator is Copyright (c) 2026 Paolo Teppa. The template is covered by the original code's MIT license, which is reproduced in `LICENSE`.
