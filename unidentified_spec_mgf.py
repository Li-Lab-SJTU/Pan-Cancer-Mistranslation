"""Remove identified spectra from an MGF file.
"""

from pathlib import Path
import os
import tempfile

import pandas as pd


FRACTION = '01CPTAC_BCprospective_W_BI_20160911_BL_f01'
MSGF_RESULT_FILE = Path('msgf_out') / f'FDR_TAR_{FRACTION}.tsv'
INPUT_MGF = Path('mgf') / f'{FRACTION}.mgf'
OUTPUT_MGF = Path('mgf') / f'{FRACTION}_filtered.mgf'


def load_spectra_to_remove(result_file):
    """Load spectrum titles that should be removed from the MGF file."""
    result_table = pd.read_csv(result_file, sep='\t', header=None)
    spectra = result_table.iloc[:, 2].dropna().astype(str).unique()
    return set(spectra)


def extract_title(title_line):
    """Extract the spectrum title from an MGF TITLE line."""
    return title_line.strip().split('.')[-2]


def remove_spectra_from_mgf(input_file, output_file, spectra_to_remove):
    """Write an MGF file with selected spectra removed.

    If the input and output paths are the same, the function writes to a
    temporary file first and then replaces the original file safely.
    """
    input_path = Path(input_file)
    output_path = Path(output_file)

    same_path = input_path.resolve() == output_path.resolve()
    target_path = output_path
    temp_path = None

    if same_path:
        temp_fd, temp_name = tempfile.mkstemp(suffix='.mgf', dir=str(output_path.parent))
        os.close(temp_fd)
        temp_path = Path(temp_name)
        target_path = temp_path

    with input_path.open('r') as infile, target_path.open('w') as outfile:
        spectrum_lines = []
        write_spectrum = True

        for line in infile:
            if line.startswith('BEGIN IONS'):
                spectrum_lines = [line]
                write_spectrum = True
            elif line.startswith('TITLE='):
                spectrum_lines.append(line)
                if extract_title(line) in spectra_to_remove:
                    write_spectrum = False
            elif line.startswith('END IONS'):
                spectrum_lines.append(line)
                if write_spectrum:
                    outfile.writelines(spectrum_lines)
            else:
                spectrum_lines.append(line)

    if temp_path is not None:
        temp_path.replace(output_path)


if __name__ == '__main__':
    spectra_to_remove = load_spectra_to_remove(MSGF_RESULT_FILE)
    remove_spectra_from_mgf(INPUT_MGF, OUTPUT_MGF, spectra_to_remove)
