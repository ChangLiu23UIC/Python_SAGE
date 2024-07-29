
def read_fasta(filename):
    """
    A function to read the fasta file into a {Uniprot_id: Sequence} dictionary format
    :param filename:
    :return:
    """
    with open(filename, 'r') as file_read:
        uniprot_dict = {}
        file_content = file_read.read()
        entries = file_content.split('>')
        for entry in entries:
            lines = entry.strip().split('\n')
            header = lines[0]
            sequence = ''.join(lines[1:])
            uniprot_id = header.split('|')[1]
            uniprot_dict[uniprot_id] = sequence
        return uniprot_dict
