# Examples

This directory contains a small example structure, `1CRN.pdb`, and two runnable
walkthrough scripts.

Run all commands in this README from the repository root unless noted.

## 1. Normal Mode Analysis

The NMA example does not require the Python mutation-analysis environment. It
only requires the Fortran build dependencies described in the main README:
`gfortran`, OpenMP support, and OpenBLAS.

```bash
cd examples
./run_nma_example.sh
```

This script builds the executable with:

```bash
make -C ..
```

and then runs:

```bash
../bin/ICed_ENM_NMA 1CRN.pdb A --mode 10 --out-prefix 1CRN_A
```

The generated NMA example outputs are:

```text
1CRN_A_eval.txt
1CRN_A_evec_CA.txt
1CRN_A_RMSF_CA.txt
1CRN_A_ref_CA.pdb
```

## 2. Mutation Impact Analysis

The mutation-impact example requires Modeller. Create and activate the Conda
environment from the repository root:

```bash
conda env create -f envs/mutation_impact_environment.yml
conda activate iced-enm-mutation
```

After installing the environment, set your Modeller license key in:

```bash
$CONDA_PREFIX/lib/modeller-10.8/modlib/modeller/config.py
```

Replace:

```python
license = 'XXXX'
```

with your Modeller license key, then verify the installation:

```bash
python3 -c "import modeller; print(modeller.__version__)"
```

Then run the example:

```bash
cd examples
./run_mutation_impact_example.sh
```

This script builds the NMA executable and runs:

```bash
python3 ../scripts/mutation_impact_analysis.py 1CRN.pdb A \
  --iced-enm-bin ../bin/ICed_ENM_NMA \
  --output-dir 1CRN_A
```

The mutation-impact example outputs are written to:

```text
examples/1CRN_A/
```

Key output files are:

```text
mutation_impact_scores.tsv
residue_impact_scores.tsv
residue_impact_bfactor.pdb
mutation_impact_scores.failed.txt
```

`mutation_impact_scores.tsv` contains mutation-level scores for individual
single-residue substitutions.

`residue_impact_scores.tsv` contains residue-level impact scores after
post-processing and smoothing.

`residue_impact_bfactor.pdb` stores the residue-level impact score in the PDB
B-factor column for visualization.

## Visualizing The B-Factor PDB

In PyMOL:

```pymol
load 1CRN_A/residue_impact_bfactor.pdb
hide everything
show cartoon
spectrum b, blue_red, minimum=0, maximum=100
```

The B-factor values are visualization-scaled to 0-100. Values above the 90th
percentile are clipped to the maximum color.

## Notes

`1CRN.pdb` is included only as a lightweight usage example. For scientific use
or redistribution, cite the original structure source.

The example output files are intentionally small so that users can quickly check
whether the code builds, Modeller is configured correctly, and the pipeline
finishes end to end.
