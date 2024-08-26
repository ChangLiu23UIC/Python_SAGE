from cleavage_methods import peptide_cleavage
from multiprocessing import Pool
from b_y_ion_spectrum import *


class Protein:
    def __init__(self, uniprot, peptides, rev, ions, ms1):
        self.uniprot = uniprot
        self.rev = rev
        self.peptides = peptides
        self.ions = ions
        self.ms1 = ms1


def process_entry(entry, rule):
    """
    Process a single FASTA entry.
    :param entry:
    :param rule:
    :return:
    """
    lines = entry.strip().split('\n')
    if len(lines) < 2 or '|' not in lines[0]:
        return None, None
    header = lines[0]
    sequence = ''.join(lines[1:])
    if "X" in sequence or "U" in sequence:
        return Protein("None", ["NONE"], False, {}, [])
    else:
        uniprot_id = header.split('|')[1]
        rev = header.split('|')[0].startswith("rev_")
        peptides = peptide_cleavage(rule, sequence)
        spectrum = {}
        ms1 = sorted([protein_weight(peptide) for peptide in peptides])
        #     spectrum[peptide] = cal_b_y_ion_mass(peptide)
        protein = Protein(uniprot_id, peptides, rev, spectrum, ms1)
        return protein


def run_fasta_to_class_parallel(rule, filename):
    """
    A function to read the fasta file into a {Uniprot_id: cleaved Sequence} dictionary format
    and process in parallel.
    :param filename:
    :return:
    """
    with open(filename, 'r') as file_read:
        file_content = file_read.read()
        entries = [entry for entry in file_content.split('\n>') if entry.strip()]

    with Pool() as pool:
        results = pool.starmap(process_entry, [(entry, rule) for entry in entries])

    # Convert the results into a dictionary, filtering out None entries
    protein_results = {protein.uniprot: protein for protein in results}

    return protein_results


if __name__ == '__main__':
    import time

    start_time = time.time()

    # Example usage
    rule = 'trypsin'
    filename = 'human_1_rev.fasta'
    uniprot_sequences = run_fasta_to_class_parallel(rule, filename)

    end_time = time.time()


    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.2f} seconds")