"""Generate a mutation FASTA file from MSGF+ output tables.
"""

import os
import re

import pandas as pd

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'
PEPTIDE_COLUMN_INDEX = 9


def collect_peptides(input_directory='msgf_out/', experiment='01'):
    """Collect unique peptides from MSGF+ result tables."""
    peptides = set()
    file_prefix = f'FDR_TAR_{experiment}'

    for filename in os.listdir(input_directory):
        if not filename.startswith(file_prefix):
            continue

        filepath = os.path.join(input_directory, filename)
        df = pd.read_csv(filepath, sep='\t', header=None)

        for peptide in df.iloc[:, PEPTIDE_COLUMN_INDEX]:
            peptides.add(re.sub(r'[0-9+.]', '', str(peptide)))

    return sorted(peptides)


def write_mutation_fasta(peptides, output_path):
    """Write all single amino acid substitution variants to FASTA."""
    with open(output_path, 'w') as fasta_file:
        for peptide in peptides:
            for position, original_aa in enumerate(peptide):
                for mutation_index, aa in enumerate(AMINO_ACIDS, start=1):
                    if aa == original_aa:
                        continue

                    mutated_peptide = peptide[:position] + aa + peptide[position + 1:]
                    fasta_file.write(f'>mutate_{peptide}_{position + 1}_{mutation_index}\n')
                    fasta_file.write(f'{mutated_peptide}\n')


def generate_mutate_fasta(experiment='01', input_directory='msgf_out/', output_directory='.'):
    """Generate the mutation FASTA for a given experiment prefix."""
    peptides = collect_peptides(input_directory=input_directory, experiment=experiment)
    output_path = os.path.join(output_directory, f'mutate_peptides_{experiment}.fasta')
    write_mutation_fasta(peptides, output_path)


if __name__ == '__main__':
    generate_mutate_fasta()
