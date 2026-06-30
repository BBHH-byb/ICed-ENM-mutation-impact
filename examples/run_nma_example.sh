#!/usr/bin/env bash
set -euo pipefail

# Run this script from examples.

make -C ..

../bin/ICed_ENM_NMA 1CRN.pdb A --mode 10 --out-prefix 1CRN_A

echo "NMA example complete."
echo "Generated: 1CRN_A_eval.txt, 1CRN_A_evec_CA.txt, 1CRN_A_RMSF_CA.txt, 1CRN_A_ref_CA.pdb"
