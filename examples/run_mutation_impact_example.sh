#!/usr/bin/env bash
set -euo pipefail

# Run this script from examples after activating the Modeller environment.

make -C ..

python3 ../scripts/mutation_impact_analysis.py 1CRN.pdb A \
  --iced-enm-bin ../bin/ICed_ENM_NMA \
  --output-dir 1CRN_A

echo "Mutation-impact example complete."
echo "Generated: 1CRN_A/mutation_impact_scores.tsv, 1CRN_A/residue_impact_scores.tsv, and 1CRN_A/residue_impact_bfactor.pdb"
