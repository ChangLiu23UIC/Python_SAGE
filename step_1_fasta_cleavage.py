from expasy import expasy_rules
import re
from multiprocessing import Pool


def peptide_cleavage(rule:str, sequence:str, mis_cleavage:int = 2, min_len:int = 7, max_len:int = 50):
    """
    Method for peptide cleavage with the expasy rule selected with regular expression. Default min len is 7 and max len
    is 50 with miscleavage of 2.
    :param rule:
    :param sequence:
    :return:
    """
    exp_rule = expasy_rules[rule]
    cut_end_pos_list = [i.start() for i in re.finditer(exp_rule, sequence)] + [(len(sequence) - 1)]
    cut_start_pos_list = [0] + [i + 1 for i in cut_end_pos_list[:-1]]

    pieces_list = {"start_pos": cut_start_pos_list, "end_pos": cut_end_pos_list}

    result = [
        sequence[cut_start_pos_list[i]:cut_end_pos_list[j]]
        for i in range(len(cut_start_pos_list))
        for j in range(i, min(i + mis_cleavage + 1, len(cut_start_pos_list)))
        if min_len <= len(sequence[cut_start_pos_list[i]:cut_end_pos_list[j]]) <= max_len
    ]

    return result


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

    with Pool() as pool:
        results = pool.starmap(process_entry, [(entry, rule) for entry in entries])

    # Convert the results into a dictionary, filtering out None entries
    uniprot_dict = {uniprot_id: cleaved_sequences for uniprot_id, cleaved_sequences in results if
                    uniprot_id is not None}

    return uniprot_dict


if __name__ == '__main__':
    import time

    start_time = time.time()

    # Example usage
    rule = 'trypsin'
    filename = 'human.fasta'
    uniprot_sequences = cleave_fasta_with_enzyme_parallel(rule, filename)

    end_time = time.time()


    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.2f} seconds")