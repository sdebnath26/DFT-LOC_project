# LOC-code

Refactored localized-orbital analysis pipeline.

## What it does
- Parses XYZ geometry files.
- Parses occupied MO coefficients from log files.
- Builds geometry-based connectivity.
- Infers bond multiplicity from shared occupied MO indices.
- Prints per-atom valency-like summaries.

## Run
```bash
cd LOC-code
python main.py --xyz data/pyrrole.xyz --log data/fb_pyrrole.log
```

Optional:
```bash
python main.py --max-occupied-mos 8
```

## Legacy copy
The pre-refactor code has been preserved under `LOC-code/legacy_backup/` with the original module layout.
