# dft-loc

**Localized-orbital bond analysis from DFT output.**

`dft-loc` parses XYZ geometry files and quantum-chemistry log files (Boys-localized MO coefficients) to infer bond types (single / double / triple) and per-atom valency summaries.

---

## Installation

```bash
pip install dft-loc
```

Or install from source:

```bash
git clone https://github.com/sdebnath26/DFT-LOC_project.git
cd DFT-LOC_project
pip install -e .
```

---

## CLI usage

```bash
# Show help
dft-loc --help

# Run bond analysis
dft-loc analyze --xyz examples/pyrrole.xyz --log examples/fb_pyrrole.log

# Limit to the last 8 occupied MOs
dft-loc analyze --xyz examples/pyrrole.xyz --log examples/fb_pyrrole.log --max-occupied-mos 8
```

Example output:
```
Loaded: examples/pyrrole.xyz (9 atoms)
Loaded: examples/fb_pyrrole.log (243 lines)

Bond summary
- 0-N:1-C        double   shared_mos=[3, 6]  d=1.371 Å
- 0-N:4-C        double   shared_mos=[3, 7]  d=1.371 Å
...

Valency summary
- 0-N   valency=3 (single=1, double=2, triple=0)
- 1-C   valency=3 (single=1, double=2, triple=0)
...
```

---

## Python API

```python
from dft_loc import (
    open_xyz,
    open_log,
    parse_mo_coefficients,
    build_atom_mo_coefficients,
    get_top_atom_coeff,
    build_connectivity,
    bond_distances,
    infer_bond_pair_indices,
    classify_bonds,
    valency_atom,
    atom_labels,
)

mol = open_xyz("examples/pyrrole.xyz")
log_lines = open_log("examples/fb_pyrrole.log")

mo_dict = parse_mo_coefficients(log_lines)
atom_coeff = build_atom_mo_coefficients(mo_dict)
top_coeff = get_top_atom_coeff(atom_coeff)

connectivity = build_connectivity(mol)
distances = bond_distances(mol, connectivity)
bond_indices = infer_bond_pair_indices(top_coeff, connectivity)
bond_types = classify_bonds(bond_indices)

for label in atom_labels(mol):
    v = valency_atom(label, bond_types)
    print(f"{label}: valency={v}")
```

---

## Build from source

```bash
pip install build
python -m build
```

This creates `dist/dft_loc-*.whl` and `dist/dft_loc-*.tar.gz`.

---

## License

MIT — see [LICENSE](LICENSE).