# Mutation Impact Analysis

`scripts/mutation_impact_analysis.py` runs a Modeller-based mutation scan and
uses ICed-ENM eigenvalues to calculate mutation impact scores.

## Requirements

The ICed-ENM executable is built with `make`.

The mutation pipeline additionally requires a Python environment with Modeller:

```bash
conda env create -f envs/mutation_impact_environment.yml
conda activate iced-enm-mutation
```

After installing Modeller, set the Modeller license key in:

```bash
$CONDA_PREFIX/lib/modeller-10.8/modlib/modeller/config.py
```

## Usage

```bash
python3 scripts/mutation_impact_analysis.py input.pdb chain \
  --iced-enm-bin bin/ICed_ENM_NMA
```

## Main Options

```text
input_pdb                  Input PDB file.
chain                      Chain ID(s), for example A or AB.
--entropy-percent X        Percentage of CA modes used for entropy. Default: 5.
--entropy-mode N           Absolute number of modes used for entropy.
--jobs N                   Number of parallel mutation/NMA jobs. Default: 4.
--nma-core N               OpenMP cores per ICed_ENM_NMA run. Default: 1.
--max-attempts N           Mutation-generation attempts per mutant. Default: 3.
--seed N                   Base seed for deterministic mutant generation.
--output-dir DIR           Output directory. Default: <pdb_id>_<chain>.
--iced-enm-bin PATH        Path to ICed_ENM_NMA.
--modeller-path PATH       Optional Modeller package/modlib/install path.
```

`--entropy-percent` and `--entropy-mode` are mutually exclusive.

During mutant NMA, total CPU use is approximately `--jobs` x `--nma-core`.
The default keeps `--nma-core 1` because multiple mutant structures are already
processed in parallel by `--jobs`.

## Output Layout

For `1CRN.pdb` chain `A`, the default output directory is:

```text
1CRN_A/
```

Main outputs:

```text
1CRN_A/mutant_pdbs/              Mutant PDB files
1CRN_A/evals/                    WT and mutant eigenvalue files
1CRN_A/logs/                     Mutation and ICed-ENM logs
1CRN_A/mutation_impact_scores.tsv
1CRN_A/residue_impact_scores.tsv
1CRN_A/residue_impact_bfactor.pdb
```

`mutation_impact_scores.tsv` contains one row per mutation:

```text
chain  resid  wt_residue  mut_residue  mutation_impact_score(|ΔS|)
```

`residue_impact_scores.tsv` contains one row per residue:

```text
chain  resid  wt_residue  residue_impact_score(Ω)
```

`residue_impact_bfactor.pdb` contains the original input structure with
visualization-scaled residue impact scores written into the B-factor column.
Values above the 90th percentile are clipped to the maximum color, and the
result is mapped to 0-100. To color it in PyMOL:

```pymol
load output_dir/residue_impact_bfactor.pdb
hide everything
show cartoon
spectrum b, blue_white_red, minimum=0, maximum=100
```
