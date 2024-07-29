from expasy import expasy_rules
import re
from concurrent.futures import ProcessPoolExecutor
import time

def peptide_cleavage(rule, sequence):
    """
    Method for peptide cleavage with the expasy rule selected with regular expression
    :param rule:
    :param sequence:
    :return:
    """
    exp_rule = expasy_rules[rule]
    cut_end_pos_list = [i.start() for i in re.finditer(exp_rule, sequence)] + [(len(sequence) - 1)]
    cut_start_pos_list = [0] + [i + 1 for i in cut_end_pos_list]
    cut_start_pos_list.pop()

    pieces_list = {"start_pos": cut_start_pos_list, "end_pos": cut_end_pos_list}

    result = [sequence[pieces_list["start_pos"][i]:pieces_list["end_pos"][i] + 1]
              for i in range(0, len(pieces_list["start_pos"]) - 1)]
    return result

def process_entry(args):
    """
    Process a single FASTA entry.
    :param args: tuple containing entry and rule
    :return:
    """
    entry, rule = args
    lines = entry.strip().split('\n')
    if len(lines) < 2 or '|' not in lines[0]:
        return None, None
    header = lines[0]
    sequence = ''.join(lines[1:])
    uniprot_id = header.split('|')[1]
    return uniprot_id, peptide_cleavage(rule, sequence)

def cleave_fasta_with_enzyme_parallel(rule, filename):
    """
    A function to read the fasta file into a {Uniprot_id: cleaved Sequence} dictionary format
    and process in parallel.
    :param filename:
    :return:
    """
    with open(filename, 'r') as file_read:
        file_content = file_read.read()
        entries = [entry for entry in file_content.split('>') if entry.strip()]

    with ProcessPoolExecutor() as executor:
        results = list(executor.map(process_entry, [(entry, rule) for entry in entries]))

    # Convert the results into a dictionary
    uniprot_dict = {uniprot_id: cleaved_sequences for uniprot_id, cleaved_sequences in results}
    return uniprot_dict

if __name__ == '__main__':
    start_time = time.time()

    # Example usage
    rule = 'trypsin'
    filename = 'human.fasta'
    uniprot_sequences = cleave_fasta_with_enzyme_parallel(rule, filename)

    end_time = time.time()

    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.2f} seconds")

    # Print results
    for uniprot_id, sequences in uniprot_sequences.items():
        print(f"UniProt ID: {uniprot_id}")
        print(f"Cleaved Sequences: {sequences}\n")
