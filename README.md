# ICed-ENM Mutation Impact

Internal-coordinate elastic network model tools for normal mode analysis and
mutation impact scoring.

## Contents

```text
src/                         Fortran ICed-ENM NMA engine
scripts/                     Python mutation impact analysis pipeline
envs/                        Conda environment for Modeller-based analysis
examples/                    NMA and mutation-impact walkthrough
docs/                        Command-line and output documentation
Makefile                     Build recipe for ICed_ENM_NMA
```

## Build

The NMA engine requires `gfortran`, OpenMP support, and OpenBLAS.

```bash
make
```

By default, the Makefile looks for OpenBLAS at:

```text
/opt/homebrew/opt/openblas
```

If OpenBLAS is installed elsewhere:

```bash
make OPENBLAS_PREFIX=/path/to/openblas
```

The executable is written to:

```text
bin/ICed_ENM_NMA
```

## Quick NMA Example

```bash
cd examples/
./run_nma_example.sh
```

This compiles the code and runs ICed-ENM NMA on chain A of `1CRN.pdb`.

## Mutation Impact Analysis

Create the Python environment:

```bash
conda env create -f envs/mutation_impact_environment.yml
conda activate iced-enm-mutation
```

Set your Modeller license key in the environment's Modeller config file before
running mutation analysis:

```bash
nano "$CONDA_PREFIX/lib/modeller-10.8/modlib/modeller/config.py"
```

Replace:

```python
license = 'XXXX'
```

with your Modeller license key, then verify:

```bash
python3 -c "import modeller; print(modeller.__version__)"
```

Example command:

```bash
python3 scripts/mutation_impact_analysis.py examples/1CRN.pdb A \
  --iced-enm-bin bin/ICed_ENM_NMA
```

Or run the example script:

```bash
cd examples
./run_mutation_impact_example.sh
```

The full mutation-impact scan performs 19 substitutions for each residue in the
selected chain, so larger proteins can take time.

The pipeline parallelizes mutant generation and mutant NMA with `--jobs`.
Each ICed-ENM run uses `--nma-core 1` by default to avoid multiplying Python
process parallelism by OpenMP parallelism. Increase `--nma-core` only when you
also reduce `--jobs` or have enough CPU cores.

## Outputs

NMA outputs are described in [docs/ICed_ENM_NMA.md](docs/ICed_ENM_NMA.md).

Mutation impact outputs are described in
[docs/mutation_impact_analysis.md](docs/mutation_impact_analysis.md).

The mutation pipeline also writes `residue_impact_bfactor.pdb`, where a
visualization-scaled residue impact score is stored in the PDB B-factor column.
Values above the 90th percentile are clipped to the maximum color and mapped to
0-100. In PyMOL:

```pymol
load residue_impact_bfactor.pdb
hide everything
show cartoon
spectrum b, blue_red, minimum=0, maximum=100
```

## Example Structure

The included example structures are provided only as usage examples. For
publication or redistribution, cite the original structure source.

## Citation

If you use this software, please cite the associated bioRxiv preprint:

```text
Lee, B. H., Scaramozzino, D., Piticchio, S., and Orellana, L.
Mutation-induced reshaping of protein conformational dynamics predicted by a coarse-grained modeling framework.
https://doi.org/10.64898/2026.03.29.715126
```

See [CITATION.cff](CITATION.cff) for citation metadata.

## License

This project is distributed under the BSD 3-Clause License. See
[LICENSE](LICENSE).
