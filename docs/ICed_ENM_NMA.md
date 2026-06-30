# ICed_ENM_NMA

`ICed_ENM_NMA` is the Fortran ICed-ENM normal mode analysis executable.

## Usage

```bash
bin/ICed_ENM_NMA input.pdb chain [options]
```

Example:

```bash
bin/ICed_ENM_NMA examples/1CRN.pdb A --mode 20 --core 4 --cutoff 8.0
```

## Main Options

```text
--core N              Number of CPU cores used by OpenMP/OpenBLAS. Default: 4.
--mode N              Number of modes to output. Default: 3. Use 0 for full DOF.
--cutoff X            Cartesian cutoff distance in Angstrom. Default: 8.0.
--out-prefix PREFIX   Prefix added to output file names.
```

## Optional Outputs

```text
--write-IC            Write evec_IC.txt.
--write-CC            Write evec_CC_raw.txt, evec_CC.txt, and ref_CC.pdb.
--write-CA-raw        Write evec_CA_raw.txt.
--write-variance      Write variance.txt and variance_cumulative.txt.
--write-all           Write all optional outputs.
```

## Movie Outputs

```text
--movie WEIGHT NFRAME Generate mode trajectory PDB files.
--movie-mode N        Generate movies for modes 1..N.
--movie-only N        Generate only mode N.
```

`--movie-mode` and `--movie-only` are mutually exclusive.

## Default Output Files

Without optional output flags, the main outputs are:

```text
eval.txt
evec_CA.txt
RMSF_CA.txt
ref_CA.pdb
```

If `--out-prefix 1CRN_A` is used, the files become:

```text
1CRN_A_eval.txt
1CRN_A_evec_CA.txt
1CRN_A_RMSF_CA.txt
1CRN_A_ref_CA.pdb
```
